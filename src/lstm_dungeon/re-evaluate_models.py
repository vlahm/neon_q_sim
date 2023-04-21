
# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-20

from pathlib import Path
from neuralhydrology.evaluation import metrics
from neuralhydrology.nh_run import start_run, continue_run, eval_run, finetune
import os
import pickle
import pandas as pd
import numpy as np
from datetime import datetime
import re

wd = Path(r.r2pyenv['wdir'])
os.chdir(wd)
config_dirs = r.r2pyenv['runset']
runs_to_reeval = r.r2pyenv['runrange']
phase = r.r2pyenv['phase']

for run in runs_to_reeval:

    runid = str(run)
    config_dir = Path(config_dirs, 'run' + runid)

    runs = os.listdir(Path(wd, 'out', 'lstm_runs'))
    try:
        rundir = [x for x in runs if re.search(phase + runid, x)][0]
    except:
        print('no run data stored for run ' + runid)
        continue
    
    rundir = wd / 'out' / 'lstm_runs' / rundir

    try:
        epoch_list = os.listdir(rundir / 'test')
    except:
        print('cant re-evaluate ' + str(run))
        continue

    epoch_list.sort()
    final_epoch = epoch_list[-1]

    old_eval = rundir / 'test' / final_epoch / 'test_results.p'
    old_metrics = rundir / 'test' / final_epoch / 'test_metrics.csv'
    os.rename(old_eval, str(old_eval) + '.bak')
    os.rename(old_metrics, str(old_metrics) + '.bak')

    try:
        eval_run(run_dir=rundir, period='test')
    except:
        os.rename(str(old_eval) + '.bak', old_eval)
        os.rename(str(old_metrics) + '.bak', old_metrics)
    
    os.rename(old_eval, rundir / 'test' / final_epoch / 'test_results_re-eval.p')
    os.rename(old_metrics, rundir / 'test' / final_epoch / 'test_metrics_re-eval.csv')
    os.rename(str(old_eval) + '.bak', old_eval)
    os.rename(str(old_metrics) + '.bak', old_metrics)

    # with open(rundir / 'test' / final_epoch / 'test_results.p', 'rb') as fp:
    #     results = pickle.load(fp)
    


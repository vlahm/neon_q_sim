
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

    eval_run(run_dir=rundir, period='test')

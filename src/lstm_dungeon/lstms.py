#!/usr/bin/env python
# coding: utf-8

# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-01

import pickle
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from neuralhydrology.evaluation import metrics
from neuralhydrology.nh_run import start_run, continue_run, eval_run, finetune
import os
import pandas as pd
import numpy as np
import yaml
from datetime import datetime
import re
from operator import itemgetter
import torch
from itertools import groupby
import functools

#wd = r.r2pyenv['wdir']
#confdir = Path(r.r2pyenv['confdir'])
#rundir = Path(r.r2pyenv['rundir'])
#strtgy = r.r2pyenv['strategy']
#config_rel = r.r2pyenv['runset']
wd = '/hpc/home/mjv22/q_sim'
confdir = Path('/hpc/home/mjv22/q_sim/lstm_configs')
rundir = Path('/hpc/home/mjv22/q_sim/lstm_runs')
#config_rel = 'runs_3003-3032' #gen search
config_rel = 'runs_3033-3812' #spec search
IS_GENERALIST = False
ENSEMBLE = False

taskID=int(os.environ['SLURM_ARRAY_TASK_ID'])

os.chdir(wd)

def format_eval_metrics(results_obj, target, prefix):

    eval_metrics = pd.DataFrame()
    for key in results_obj:
        #           [basin][res][xarray][each metric specified in config]
        qobs = results_obj[key]['1D']['xr'][target + '_obs']
        qsim = results_obj[key]['1D']['xr'][target + '_sim']

        em = metrics.calculate_all_metrics(qobs.isel(time_step=-1), qsim.isel(time_step=-1))
        em = pd.DataFrame.from_dict(em, orient='index')#.reset_index()
        em.iloc[:, 0] = em.iloc[:, 0].round(3)
        eval_metrics = pd.concat([eval_metrics, em], axis=1)

    eval_metrics.columns = list(results_obj.keys())

    eval_metrics = eval_metrics.transpose()
    eval_metrics.columns = [prefix + '_' + x for x in eval_metrics.columns]
    eval_metrics = eval_metrics.reset_index()
    eval_metrics.rename(columns={'index':'site_code'}, inplace=True)

    return eval_metrics

all_runs = os.listdir(confdir)
config_dir_or_dirs = Path(confdir, config_rel)

id_range = re.match('^runs_([0-9]+)-([0-9]+)$', config_rel).groups()
id_range = list(range(int(id_range[0]), int(id_range[1]) + 1))

run = id_range[taskID - 1]
runid = str(run)

config_dir = Path(confdir, config_dir_or_dirs.stem, 'run' + runid)

results_ft1 = None
if IS_GENERALIST or ENSEMBLE:

    config_file1 = Path(config_dir, 'continue' + runid + '.yml')
    with open(config_file1, 'r') as fp:
        config_deets1 = yaml.safe_load(fp)

    finetune(config_file=config_file1)
    results = os.listdir(rundir)
    finetune1 = [x for x in results if re.match('run' + runid, x)]
    if len(finetune1) != 1:
        raise ValueError('need exactly one output directory for ' + str(rundir) + '/run' + runid)
        sys.exit(1)

    eval_run(run_dir=finetune1, period='test')
    epoch_list = os.listdir(Path(finetune1, 'test'))
    epoch_list.sort()
    final_epoch_ft = epoch_list[-1]

    with open(finetune1 / 'test' / final_epoch_ft / 'test_results.p', 'rb') as fp:
        results_ft1 = pickle.load(fp)

results_ft2 = None
if not IS_GENERALIST:

    config_file2 = Path(config_dir, 'finetune' + runid + '.yml')
    with open(config_file2, 'r') as fp:
        config_deets2 = yaml.safe_load(fp)

    finetune(config_file2)
    results = os.listdir(rundir)
    finetune2 = [x for x in results if re.match('finetune' + runid, x)]
    if len(finetune2) != 1:
        raise ValueError('need exactly one output directory for ' + str(rundir) + '/finetune' + runid)
        sys.exit(1)

    eval_run(run_dir=finetune2, period='test')
    epoch_list = os.listdir(Path(finetune2, 'test'))
    epoch_list.sort()
    final_epoch_ft = epoch_list[-1]

    with open(finetune2 / 'test' / final_epoch_ft / 'test_results.p', 'rb') as fp:
        results_ft2 = pickle.load(fp)

## compile specs

if results_ft2 is not None:
    results_ft = results_ft2
    #config_deets = config_deets2
    target = config_deets2['target_variables'][0]
else:
    results_ft = results_ft1
    #config_deets = config_deets1
    target = config_deets1['target_variables'][0]

eval_metrics_fine = format_eval_metrics(results_ft, target=target, prefix='fine')
#subset discharge_obs from results_ft
discharge_obs = results_ft[list(results_ft.keys())[0]]['1D']['xr'][target + '_obs'].to_dataframe().reset_index('date')
#retrieve the max date from discharge_obs
max_date = discharge_obs['date'].max()
#filter NaNs from discharge_obs
discharge_obs = discharge_obs[discharge_obs[target + '_obs'].notna()]

## record config deets, eval metrics, and notes for this run

with open(Path(confdir, config_dir_or_dirs.stem, 'addtl_notes.yml'), 'r') as f:
    addtl_notes = yaml.safe_load(f)
        
#addtl_notes.pop('run_dir')

results_file = Path('accumulated_results.csv')
if os.path.exists(results_file):
    accum_deets = pd.read_csv(results_file)
else:
    accum_deets = pd.DataFrame()

newrow = pd.DataFrame(data = {'datetime': [datetime.now().strftime('%Y-%m-%d %H:%M:%S')]})

if results_ft1 is not None:

    config_deets1.pop('train_basin_file')
    config_deets1.pop('validation_basin_file')
    config_deets1.pop('test_basin_file')
    config_deets1.pop('per_basin_train_periods_file')
    config_deets1.pop('per_basin_validation_periods_file')
    config_deets1.pop('per_basin_test_periods_file')
    config_deets1.pop('device')

    for i, key in enumerate(config_deets1, start=1):
        newrow.insert(i, key, str(config_deets1[key]))

    for i, key in enumerate(addtl_notes, start=2):
        newrow.insert(i, key, str(addtl_notes[key]))

if results_ft2 is not None:

    config_deets2.pop('train_basin_file')
    config_deets2.pop('validation_basin_file')
    config_deets2.pop('test_basin_file')
    config_deets2.pop('per_basin_train_periods_file')
    config_deets2.pop('per_basin_validation_periods_file')
    config_deets2.pop('per_basin_test_periods_file')
    config_deets2['learning_rateFINE'] = config_deets2.pop('learning_rate')
    config_deets2['epochsFINE'] = config_deets2.pop('epochs')
    config_deets2['experiment_nameFINE'] = config_deets2.pop('experiment_name')

    for i, key in enumerate(config_deets2, start=1):
        if key in newrow.columns:
            kk = key + '2'
        else:
            kk = key
        newrow.insert(i, kk, str(config_deets2[key]))

    eval_metrics = eval_metrics_fine
    newrow = pd.concat([newrow, eval_metrics], axis=1)

accum_deets = pd.concat([accum_deets, newrow])
accum_deets.to_csv(results_file, index=False, encoding='utf-8', header=True)

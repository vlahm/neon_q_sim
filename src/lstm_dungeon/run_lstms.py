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

## setup
wd = r.r2pyenv['wdir']
confdir = Path(r.r2pyenv['confdir'])
rundir = Path(r.r2pyenv['rundir'])
strtgy = r.r2pyenv['strategy']
config_rel = r.r2pyenv['runset']

os.chdir(wd)

IS_GENERALIST = True if strtgy == 'generalist' else False

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

    # eval_metrics.to_csv(Path(most_recent_rund, 'eval_metrics.csv'),
    #                     index=False,
    #                     encoding='utf-8',
    #                     header=False)

    eval_metrics = eval_metrics.transpose()
    eval_metrics.columns = [prefix + '_' + x for x in eval_metrics.columns]
    eval_metrics = eval_metrics.reset_index()
    eval_metrics.rename(columns={'index':'site_code'}, inplace=True)

    return eval_metrics

all_runs = os.listdir(confdir)

combined_run_bool = [not bool(re.match('run([0-9]+)', x)) for x in all_runs]
combined_run_dirs = [x for (x, y) in zip(all_runs, combined_run_bool) if y]
combined_run_dirs = list(zip(combined_run_dirs,
    [int(re.match('runs_[0-9]+-([0-9]+)$', x).group(1)) for x in combined_run_dirs]))
normal_run_dirs = [x for (x, y) in zip(all_runs, combined_run_bool) if not y]
normal_run_dirs = list(zip(normal_run_dirs,
    [int(re.match('run([0-9]+)$', x).group(1)) for x in normal_run_dirs]))

# config_rel = max(combined_run_dirs + normal_run_dirs, key=itemgetter(1))[0]
config_dir_or_dirs = Path(confdir, config_rel)

runs_in_dir = [x for x in os.listdir(config_dir_or_dirs) if re.match('^run[0-9]+$', x)]
multirun_dir = bool(len(runs_in_dir))
cfg_fls = os.listdir(config_dir_or_dirs)
turboroutine = any([re.search('(pretrain.yml|pretrained_model_loc.txt)$', x) for x in cfg_fls])

if multirun_dir:
    id_range = re.match('^runs_([0-9]+)-([0-9]+)$', config_rel).groups()
    id_range = list(range(int(id_range[0]), int(id_range[1]) + 1))
    
    cfg_fls2 = os.listdir(Path(config_dir_or_dirs, 'run' + str(id_range[0])))
else:
    id_range = [int(re.match('run([0-9]+)', config_rel).group(1))]

with open(Path(rundir, 'failed_runs.txt'), 'a') as f:
    f.write('---\n')

transfer_learning = os.path.exists(os.path.join(config_dir_or_dirs, 'pretrained_model_loc.txt'))

if turboroutine and not transfer_learning:

    pre_config = Path(config_dir_or_dirs, 'pretrain.yml')
    start_run(config_file=pre_config)
    turborund = sorted(rundir.iterdir(), key=os.path.getmtime)[-1]

if transfer_learning:

    with open(os.path.join(config_dir_or_dirs, 'pretrained_model_loc.txt'), 'r') as pm:
        pretrained_model_loc = pm.read().splitlines()[0]

for run in id_range:

    runid = str(run)
    path_extra = config_dir_or_dirs.stem if multirun_dir else ''
    config_dir = Path(confdir, path_extra, 'run' + runid)
    if turboroutine:
        config_file = Path(config_dir, 'continue' + runid + '.yml')
    else:
        config_file = Path(config_dir, 'run' + runid + '.yml')
    nhm_config_file = Path(config_dir, 'run' + runid + 'nhm.yml')
    finetune_config_file = Path(config_dir, 'finetune' + runid + '.yml')
    is_finetune_run = os.path.exists(finetune_config_file)

    try:

        ## real data initial train OR continue train

        with open(config_file, 'r') as fp:
            config_deets = yaml.safe_load(fp)

        if turboroutine and not transfer_learning:
            with open(config_file, 'a') as fp:
                fp.write(f'base_run_dir: {turborund}')
            finetune(config_file=config_file)
        elif turboroutine and transfer_learning:
            with open(config_file, 'a') as fp:
                fp.write(f'base_run_dir: {pretrained_model_loc}')
            finetune(config_file=config_file)
        elif not transfer_learning:
            #train base model and save output directory
            start_run(config_file=config_file)

        most_recent_rund = sorted(rundir.iterdir(), key=os.path.getmtime)[-1]

        ## finetune model

        if is_finetune_run and not IS_GENERALIST:

            with open(Path(config_dir, 'finetune' + runid + '.yml'), 'a') as fp:
                fp.write(f'base_run_dir: {most_recent_rund}')

            finetune(Path(config_dir, 'finetune' + runid + '.yml'))
            finetune_dir = sorted(rundir.iterdir(), key=os.path.getmtime)[-1]

    except Exception as e:
        print(e)
        with open(Path(wd, 'failed_runs.txt'), 'a') as f:
            f.write(str(run) + '\n')

    ## evaluate fine-tuned run on test set; load predictions

    results_ft = None
    if is_finetune_run and not IS_GENERALIST:

        eval_run(run_dir=finetune_dir, period='test')

        epoch_list = os.listdir(Path(finetune_dir, 'test'))
        epoch_list.sort()
        final_epoch_ft = epoch_list[-1]

        with open(finetune_dir / 'test' / final_epoch_ft / 'test_results.p', 'rb') as fp:
            results_ft = pickle.load(fp)

    if (turboroutine and not is_finetune_run) or (turboroutine and IS_GENERALIST):

        eval_run(run_dir=most_recent_rund, period='test')

        epoch_list = os.listdir(Path(most_recent_rund, 'test'))
        epoch_list.sort()
        final_epoch_ft = epoch_list[-1]

        with open(most_recent_rund / 'test' / final_epoch_ft / 'test_results.p', 'rb') as fp:
            results_ft = pickle.load(fp)

    ## compute all metrics implemented in the neuralHydrology package.
    ## (additional hydrological signatures implemented in neuralhydrology.evaluation.signatures)

    target = config_deets['target_variables'][0]

    if results_ft is not None:
        eval_metrics_fine = format_eval_metrics(results_ft, target=target, prefix='fine')

        #subset discharge_obs from results_ft
        discharge_obs = results_ft[list(results_ft.keys())[0]]['1D']['xr'][target + '_obs'].to_dataframe().reset_index('date')
        #retrieve the max date from discharge_obs
        max_date = discharge_obs['date'].max()
        #filter NaNs from discharge_obs
        discharge_obs = discharge_obs[discharge_obs[target + '_obs'].notna()]

    ## record config deets, eval metrics, and notes for this run
    
    if multirun_dir:
        with open(Path(confdir, path_extra, 'addtl_notes.yml'), 'r') as f:
            addtl_notes = yaml.safe_load(f)
    else:
        with open(Path(confdir, 'run' + runid, 'addtl_notes.yml'), 'r') as f:
            addtl_notes = yaml.safe_load(f)
            
    addtl_notes.pop('run_dir')

    results_file = Path(rundir, 'accumulated_results.csv')
    if os.path.exists(results_file):
        accum_deets = pd.read_csv(results_file)
    else:
        accum_deets = pd.DataFrame()

    newrow = pd.DataFrame(data = {'datetime': [datetime.now().strftime('%Y-%m-%d %H:%M:%S')]})

    config_deets.pop('train_basin_file')
    config_deets.pop('validation_basin_file')
    config_deets.pop('test_basin_file')
    config_deets.pop('per_basin_train_periods_file')
    config_deets.pop('per_basin_validation_periods_file')
    config_deets.pop('per_basin_test_periods_file')
    config_deets.pop('device')

    for i, key in enumerate(config_deets, start=1):
        newrow.insert(i, key, str(config_deets[key]))

    for i, key in enumerate(addtl_notes, start=2):
        newrow.insert(i, key, str(addtl_notes[key]))

    if results_ft is not None:

        if turboroutine and not is_finetune_run:
            with open(Path(config_dir, 'continue' + runid + '.yml'), 'r') as fp:
                finetune_deets = yaml.safe_load(fp)
        else:
            with open(Path(config_dir, 'finetune' + runid + '.yml'), 'r') as fp:
                finetune_deets = yaml.safe_load(fp)

        finetune_deets.pop('train_basin_file')
        finetune_deets.pop('validation_basin_file')
        finetune_deets.pop('test_basin_file')
        finetune_deets.pop('per_basin_train_periods_file')
        finetune_deets.pop('per_basin_validation_periods_file')
        finetune_deets.pop('per_basin_test_periods_file')
        finetune_deets['learning_rateFINE'] = finetune_deets.pop('learning_rate')
        finetune_deets['epochsFINE'] = finetune_deets.pop('epochs')
        finetune_deets['experiment_nameFINE'] = finetune_deets.pop('experiment_name')

        for i, key in enumerate(finetune_deets, start=1):
            if key in newrow.columns:
                kk = key + '2'
            else:
                kk = key
            newrow.insert(i, kk, str(finetune_deets[key]))

        eval_metrics = eval_metrics_fine
        newrow = pd.concat([newrow, eval_metrics], axis=1)
    else:
        raise ValueError('nothing to add to eval metrics')

    accum_deets = pd.concat([accum_deets, newrow])
    accum_deets.to_csv(results_file, index=False, encoding='utf-8', header=True)

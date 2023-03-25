#!/usr/bin/env python
# coding: utf-8

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
confdir = r.r2pyenv['confdir']
rundir = r.r2pyenv['rundir']
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

all_runs = os.listdir(Path(confdir))

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
    semiturboroutine = any([re.search('pretrain.yml$', x) for x in cfg_fls2])
else:
    id_range = [int(re.match('run([0-9]+)', config_rel).group(1))]

with open(Path(rundir, 'failed_runs.txt'), 'a') as f:
    f.write('---\n')

transfer_learning = os.path.exists(os.path.join(config_dir_or_dirs, 'pretrained_model_loc.txt'))

if turboroutine and not transfer_learning:

    pre_config = Path(config_dir_or_dirs, 'pretrain.yml')
    start_run(config_file=pre_config)
    turborund = sorted(Path(wd, 'runs').iterdir(), key=os.path.getmtime)[-1]
    # turborund = wd / Path('runs', 'runs_1260-1328_pretrain_1307_171204')

if transfer_learning:

    with open(os.path.join(config_dir_or_dirs, 'pretrained_model_loc.txt'), 'r') as pm:
        pretrained_model_loc = pm.read().splitlines()[0]

skip_continue = False

for run in id_range:

    runid = str(run)
    path_extra = config_dir_or_dirs.stem if multirun_dir else ''
    config_dir = Path(wd, 'run_configs', path_extra, 'run' + runid)
    if turboroutine:
        config_file = Path(config_dir, 'continue' + runid + '.yml')
    elif semiturboroutine:
        config_file = Path(config_dir, 'pretrain' + '.yml')
    else:
        config_file = Path(config_dir, 'run' + runid + '.yml')
    evaluating_initial_run = os.path.exists(Path(config_dir, '..', 'test_ranges.pkl')) #misnomer. nhm run might precede
    nhm_config_file = Path(config_dir, 'run' + runid + 'nhm.yml')
    finetune_config_file = Path(config_dir, 'finetune' + runid + '.yml')
    is_nhm_run = os.path.exists(nhm_config_file)
    is_finetune_run = os.path.exists(finetune_config_file)

    try:

        ## [OPTIONAL] NHM pre-train

        if is_nhm_run:

            with open(nhm_config_file, 'r') as fp:
                nhm_config_deets = yaml.safe_load(fp)

            start_run(config_file=nhm_config_file)
            nhm_rund = sorted(Path(wd, 'runs').iterdir(), key=os.path.getmtime)[-1]

        ## real data initial train OR continue train

        # config_file = Path('/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs/runs_424-453/run424/run424.yml')
        with open(config_file, 'r') as fp:
            config_deets = yaml.safe_load(fp)

        if not skip_continue:
            if is_nhm_run:
                nhm_rund = sorted(Path(wd, 'runs').iterdir(), key=os.path.getmtime)[-1]
                continue_run(run_dir=nhm_rund, config_file=config_file)
            elif turboroutine and not transfer_learning:
                with open(config_file, 'a') as fp:
                    fp.write(f'base_run_dir: {turborund}')
                # continue_run(run_dir=turborund, config_file=config_file)
                finetune(config_file=config_file)
            elif turboroutine and transfer_learning:
                with open(config_file, 'a') as fp:
                    fp.write(f'base_run_dir: {pretrained_model_loc}')
                finetune(config_file=config_file)
            elif not transfer_learning:
            ## train base model and save output directory
                start_run(config_file=config_file)

        if not skip_continue:
            most_recent_rund = sorted(Path(wd, 'runs').iterdir(), key=os.path.getmtime)[-1]
        else:
            most_recent_rund = Path('/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/runs/run1553_2309_164410')

        ## [OPTIONAL] finetune model

        if is_finetune_run and not IS_GENERALIST:

            with open(Path(config_dir, 'finetune' + runid + '.yml'), 'a') as fp:
                fp.write(f'base_run_dir: {most_recent_rund}')

            finetune(Path(config_dir, 'finetune' + runid + '.yml'))
            finetune_dir = sorted(Path(wd, 'runs').iterdir(), key=os.path.getmtime)[-1]

        if semiturboroutine:

            with open(Path(config_dir, 'continue' + runid + '.yml'), 'a') as fp:
                fp.write(f'base_run_dir: {most_recent_rund}')

            finetune(Path(config_dir, 'continue' + runid + '.yml'))
            finetune_dir = sorted(Path(wd, 'runs').iterdir(), key=os.path.getmtime)[-1]

    except Exception as e:
        print(e)
        with open(Path(wd, 'failed_runs.txt'), 'a') as f:
            f.write(str(run) + '\n')

    # base_model_nse = base_model_results[basin]['1D']['NSE']
    # finetune_nse = finetuned_results[basin]['1D']['NSE']
    # print(f'Basin {basin} base model performance: {base_model_nse:.3f}')
    # print(f'Performance after finetuning: {finetune_nse:.3f}')

    ## evaluate fine-tuned run on test set; load predictions

    results_ft = None
    # if re.search('finetune', str(finetune_dir)):
    if (is_finetune_run and not IS_GENERALIST) or semiturboroutine:

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

    ## evaluate base run on test set; load predictions

    if evaluating_initial_run:

        eval_run(run_dir=most_recent_rund, period='test')

        epoch_list = os.listdir(Path(most_recent_rund, 'test'))
        epoch_list.sort()
        final_epoch = epoch_list[-1]

        with open(most_recent_rund / 'test' / final_epoch / 'test_results.p', 'rb') as fp:
            results = pickle.load(fp)


    ## evaluate nhm run on test set; load predictions

    if is_nhm_run:

        eval_run(run_dir=nhm_rund, period='test')

        epoch_list_nhm = os.listdir(Path(nhm_rund, 'test'))
        epoch_list_nhm.sort()
        final_epoch_nhm = epoch_list_nhm[-1]

        with open(nhm_rund / 'test' / final_epoch_nhm / 'test_results.p', 'rb') as fp:
            results_nhm = pickle.load(fp)

    ## plot Q obs vs. fit
    target = config_deets['target_variables'][0]

    # all(np.isnan(results['FLNT']['1D']['xr']['discharge_sim']))
    # all(np.isnan(results_ft['FLNT']['1D']['xr']['discharge_sim']))
    # for xx in results:
    #     print(np.isnan(results[xx]['1D']['xr']['discharge_obs']).to_pandas().sum().loc[0])

    # if manual_rerun:
    #     pdf_name = 'Qobs_Qpred2.pdf'
    # else:
    pdf_name = 'Qobs_Qpred.pdf'

    pdf = matplotlib.backends.backend_pdf.PdfPages(Path(most_recent_rund, pdf_name))

    qsim_lst = []
    if evaluating_initial_run:

        for key in results:

            metr = 'NSE' if any([x == 'NSE' for x in list(results.values())[0]['1D'].keys()]) else 'discharge_NSE'

            if len(results[key]['1D']['xr'][target + '_sim'].dims) == 3:
                results[key]['1D']['xr'] = results[key]['1D']['xr'].mean(dim='samples', skipna=True)

            #              [basin][res][xarray][each metric specified in config]
            qsim = results[key]['1D']['xr'][target + '_sim']
            qsim.values[qsim.values < 0] = 0
            qobs = results[key]['1D']['xr'][target + '_obs']

            qsim_pd = qsim.to_dataframe().droplevel(1).reset_index('date').rename(columns={'discharge_sim': key})
            qsim_lst.append(qsim_pd)

            # qobs = results[key]['1D']['xr']['QObs(mm/d)_obs']
            # qsim = results[key]['1D']['xr']['QObs(mm/d)_sim']
            # if len(qobs[0]) > 1:
            #     [x[0] for x in np.array(qobs)]

            fig, ax = plt.subplots(figsize=(16, 10))
            ax.plot(qobs['date'], qobs)
            ax.plot(qsim['date'], qsim)

            if results_ft is not None and len(results_ft[key]) > 0:
                qsim_ft = results_ft[key]['1D']['xr'][target + '_sim']
                qsim_ft.values[qsim_ft.values < 0] = 0

                ax.plot(qsim_ft['date'], qsim_ft)
                ax.legend(['obs', 'pred_base', 'pred_fine'])
                ax.set_title('Site ' + key +\
                             f" test period\nbase NSE = {results[key]['1D'][metr]:.3f}" +\
                             f"; finetuned NSE = {results_ft[key]['1D'][metr]:.3f}")
            else:
                ax.legend(['obs', 'pred'])
                ax.set_title('Site ' + key + f" test period; NSE = {results[key]['1D']['NSE']:.3f}")

            # ax.set_ylabel('Discharge (L/s)')
            ax.set_ylabel('Discharge (mm/d)')
            fig = fig.get_figure()
            pdf.savefig(fig)
    else:

        # metr = 'NSE' if any([x == 'NSE' for x in list(results_ft.values())[0]['1D'].keys()]) else 'discharge_NSE'

        for key in results_ft:

            #           [basin][res][xarray][each metric specified in config]
            if len(results_ft[key]['1D']['xr'][target + '_sim'].dims) == 3: #average mc_dropout samples
                results_ft[key]['1D']['xr'] = results_ft[key]['1D']['xr'].mean(dim='samples', skipna=True)
                # qsim = results_ft[key]['1D']['xr'].mean(dim='samples', skipna=True)[target + '_sim']
            qsim = results_ft[key]['1D']['xr'][target + '_sim']
            qsim.values[qsim.values < 0] = 0
            qobs = results_ft[key]['1D']['xr'][target + '_obs']

            qsim_pd = qsim.to_dataframe().droplevel(1).reset_index('date').rename(columns={'discharge_sim': key})
            qsim_lst.append(qsim_pd)

            fig, ax = plt.subplots(figsize=(16, 10))
            ax.plot(qobs['date'], qobs)
            ax.plot(qsim['date'], qsim)
            ax.legend(['obs', 'pred'])
            skill = 'NA'
            try:
                skill = str(round(results_ft[key]['1D']['NSE'], 3))
            except KeyError:
                try:
                    skill = str(round(results_ft[key]['1D']['discharge_NSE'], 3))
                except KeyError:
                    pass
            ax.set_title('Site ' + key + " finetune test period; NSE = " + skill)
            ax.set_ylabel('Discharge (mm/d)')
            fig = fig.get_figure()
            pdf.savefig(fig)

    pdf.close()

    # functools.reduce(lambda x, y: pd.merge(x, y, 'left', on='date'), qsim_lst)\
    # .to_csv('/home/mike/Desktop/nhc_unhc_qsim/' + str(run) + '.csv')

    ## compute all metrics implemented in the neuralHydrology package.
    ## (additional hydrological signatures implemented in neuralhydrology.evaluation.signatures)
    if evaluating_initial_run:
        eval_metrics_base = format_eval_metrics(results, target=target, prefix='base')

    if results_ft is not None:
        eval_metrics_fine = format_eval_metrics(results_ft, target=target, prefix='fine')

        #subset discharge_obs from results_ft
        discharge_obs = results_ft[list(results_ft.keys())[0]]['1D']['xr'][target + '_obs'].to_dataframe().reset_index('date')
        #retrieve the max date from discharge_obs
        max_date = discharge_obs['date'].max()
        #filter NaNs from discharge_obs
        discharge_obs = discharge_obs[discharge_obs[target + '_obs'].notna()]




    if is_nhm_run:
        eval_metrics_nhm = format_eval_metrics(results_nhm, target=target, prefix='NHM')

    ## record config deets, eval metrics, and notes for this run
    if multirun_dir:
        with open(Path(wd, 'run_configs', path_extra, 'addtl_notes.yml'), 'r') as f:
            addtl_notes = yaml.safe_load(f)
    else:
        with open(Path(wd, 'run_configs', 'run' + runid, 'addtl_notes.yml'), 'r') as f:
            addtl_notes = yaml.safe_load(f)

    if os.path.exists('../accumulated_results.csv'):
        accum_deets = pd.read_csv('../accumulated_results.csv')
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

        if (turboroutine and not is_finetune_run) or semiturboroutine:
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

        if evaluating_initial_run:
            eval_metrics = eval_metrics_base.merge(eval_metrics_fine, how='left', on='site_code')
        else:
            eval_metrics = eval_metrics_fine
        newrow = pd.concat([newrow, eval_metrics], axis=1)
    else:
        if evaluating_initial_run:
            newrow = pd.concat([newrow, eval_metrics_base], axis=1)
        else:
            raise ValueError('nothing to add to eval metrics')

    if is_nhm_run:

        nhm_config_deets = {'learning_rateNHM': config_deets['learning_rate'],
                            'epochsFINE': config_deets['epochs']}

        eval_metrics = eval_metrics_base.merge(eval_metrics_fine, how='left', on='site_code')

    accum_deets = pd.concat([accum_deets, newrow])
    accum_deets.to_csv(Path(wd, '../accumulated_results.csv'),
        index=False, encoding='utf-8', header=True)

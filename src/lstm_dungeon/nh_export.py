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

## setup
wd = Path('/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods')
os.chdir(wd)

with open('../../data/nh_methods/run_name_courier.txt', 'r') as f:
    collection_name = f.read().split('\n')[0]
bestruns = pd.read_csv(wd / '../../out/best_run_for_each_site' / Path(collection_name + '.csv'))
# bestruns = bestruns[bestruns['fine_NSE'] >= 0.25]
# bestruns.reset_index(inplace = True, drop = True)
# bestruns.site_code = bestruns.site_code.apply(lambda x: x + '_extrapolate')

all_runs = os.listdir(wd / 'runs/')

for i in bestruns.index:

    site, runid, NSE = bestruns.iloc[i, :]
    if not bool(re.match(r'.*?_GAPPED$', site)):
        site = site + '_extrapolate'

    rundir = [x for x in all_runs if re.match('finetune' + str(runid) + '_', x)]
    # if len(rundir) > 1:
    #     raise ValueError('two runs with the same id')
    if not len(rundir):
        rundir = [x for x in all_runs if re.match('run' + str(runid) + '_', x)]
    rundir.sort()
    rundir = [rundir[len(rundir) - 1]]
    # if len(rundir) > 1:
    #     raise ValueError('two runs with the same id')
    testdir = wd / 'runs' / rundir[0] / 'test'
    epoch_list = os.listdir(testdir)
    epoch_list.sort()
    final_epoch = epoch_list[-1]

    try:
        with open(testdir / final_epoch / 'test_results.p', 'rb') as fp:
            results = pickle.load(fp)
        results
    except:
        print('failed at least one iteration. need to repeat this with old version of xarray')
        continue

    if len(results[site]['1D']['xr']['discharge_sim'].dims) == 3:
        results[site]['1D']['xr'] = results[site]['1D']['xr'].mean(dim='samples', skipna=True)
    qsim = results[site]['1D']['xr']['discharge_sim']
    qsim.values[qsim.values < 0] = 0
    qdate = results[site]['1D']['xr']['date']

    df = results[site]['1D']['xr'].to_dataframe()
    df.reset_index(inplace=True)

    if bool(re.match(r'.*?_GAPPED$', site)):
        outdir = wd / '../../out/predicted_q' / Path(str(collection_name))
    else:
        outdir = wd / '../../out/predicted_q' / Path(str(collection_name) + '_extrapolate')
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    df.loc[:, ['date', 'discharge_sim', 'discharge_obs']].to_csv(
        outdir / Path(site + '_' + str(runid) + '.csv'),
        index=False,
        encoding='utf-8')

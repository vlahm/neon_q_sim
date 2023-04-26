#!/usr/bin/env python
# coding: utf-8

# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-25

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

wd = r.r2pyenv['wdir']
confdir = Path(r.r2pyenv['confdir'])
rundir = Path(r.r2pyenv['rundir'])
IS_GENERALIST = True if r.r2pyenv['strategy'] == 'generalist' else Dalge
ENSEMBLE = r.r2pyenv['ensemble']
config_rel = r.r2pyenv['runset']

os.chdir(wd)

all_runs = os.listdir(confdir)
config_dir_or_dirs = Path(confdir, config_rel)

id_range = re.match('^runs_([0-9]+)-([0-9]+)$', config_rel).groups()
id_range = list(range(int(id_range[0]), int(id_range[1]) + 1))

for run in id_range:

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

    results_ft2 = None
    if not IS_GENERALIST:

        config_file2 = Path(config_dir, 'finetune' + runid + '.yml')
        with open(config_file2, 'r') as fp:
            config_deets2 = yaml.safe_load(fp)

        finetune(config_file=config_file2)
        results = os.listdir(rundir)
        finetune2 = [x for x in results if re.match('finetune' + runid, x)]
        if len(finetune2) != 1:
            raise ValueError('need exactly one output directory for ' + str(rundir) + '/finetune' + runid)


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

wd = '/hpc/path/to/working/directory'
confdir = Path('/hpc/path/to/lstm_configs')
rundir = Path('/hpc/path/to/lstm_runs')
config_rel = 'runs_3843-4112' #e.g.
IS_GENERALIST = False #e.g.
ENSEMBLE = False #e.g.

taskID=int(os.environ['SLURM_ARRAY_TASK_ID'])

os.chdir(wd)

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


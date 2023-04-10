import pickle
from pathlib import Path
import os
import pandas as pd
import re
import numpy as np
from operator import itemgetter
# from datetime import datetime

try:
    del ft_test_ranges
except:
    pass

try:
    del test_ranges
except:
    pass

# wd = r.r2pyenv['wdir']
config_dir = Path(r.r2pyenv['confdir'])
# strtgy = r.r2pyenv['strategy']
curr_config_rel = r.r2pyenv['runset']

all_runs = os.listdir(config_dir)

curr_config_dir = Path(config_dir, curr_config_rel)
finetuning_special_ranges = os.path.exists(Path(curr_config_dir, 'finetune_train_ranges.csv'))

cfg_fls = os.listdir(curr_config_dir)
finetune_dirs = [x for x in cfg_fls if re.match('finetune[0-9]+_[0-9]{3}', x)]
extra_finetunes = False
if len(finetune_dirs):
    extra_finetunes = True
    latest_finetune = sorted(finetune_dirs)[-1]
turboroutine = any([re.search('(pretrain.yml|pretrained_model_loc.txt)$', x) for x in cfg_fls])

a_rundir = [x for x in cfg_fls if re.match('^run[0-9]+$', x)][0]
cfg_fls2 = os.listdir(Path(curr_config_dir, str(a_rundir)))
semiturboroutine = any([re.search('pretrain.yml$', x) for x in cfg_fls2])

nhm_run = any([x for x in cfg_fls if re.match('^nhm_', x)])

## load preformatted CSVs for train/validate/test date ranges
train_ranges = pd.read_csv(Path(curr_config_dir, 'train_ranges.csv'),
                           dtype={'basin_id': str},
                           parse_dates=[1,2])

validation_ranges = pd.read_csv(Path(curr_config_dir, 'validation_ranges.csv'),
                                dtype={'basin_id': str},
                                parse_dates=[1,2])

try:
    test_ranges = pd.read_csv(Path(curr_config_dir, 'test_ranges.csv'),
                              dtype={'basin_id': str},
                              parse_dates=[1,2])
except:
    pass

if finetuning_special_ranges:
    ft_train_ranges = pd.read_csv(Path(curr_config_dir, 'finetune_train_ranges.csv'),
                                  dtype={'basin_id': str},
                                  parse_dates=[1,2])

    ft_validation_ranges = pd.read_csv(Path(curr_config_dir, 'finetune_validation_ranges.csv'),
                                       dtype={'basin_id': str},
                                       parse_dates=[1,2])

    try:
        ft_test_ranges = pd.read_csv(Path(curr_config_dir, 'finetune_test_ranges.csv'),
                                     dtype={'basin_id': str},
                                     parse_dates=[1,2])
    except:
        pass

if extra_finetunes:
    ftx_train_ranges = pd.read_csv(Path(curr_config_dir, latest_finetune, 'train_ranges.csv'),
                                   dtype={'basin_id': str},
                                   parse_dates=[1,2])

    ftx_validation_ranges = pd.read_csv(Path(curr_config_dir, latest_finetune, 'validation_ranges.csv'),
                                        dtype={'basin_id': str},
                                        parse_dates=[1,2])

    try:
        ftx_test_ranges = pd.read_csv(Path(curr_config_dir, latest_finetune, 'test_ranges.csv'),
                                      dtype={'basin_id': str},
                                      parse_dates=[1,2])
    except:
        pass

if nhm_run:
    nhm_train_ranges = pd.read_csv(Path(curr_config_dir, 'nhm_train_ranges.csv'),
                                   dtype={'basin_id': str},
                                   parse_dates=[1,2])

    nhm_validation_ranges = pd.read_csv(Path(curr_config_dir, 'nhm_validation_ranges.csv'),
                                        dtype={'basin_id': str},
                                        parse_dates=[1,2])

    try:
        nhm_test_ranges = pd.read_csv(Path(curr_config_dir, 'nhm_test_ranges.csv'),
                                      dtype={'basin_id': str},
                                      parse_dates=[1,2])
    except:
        pass

## convert to dicts and pickle
train_ranges = train_ranges.set_index('basin_id').to_dict(orient='index')
validation_ranges = validation_ranges.set_index('basin_id').to_dict(orient='index')
try:
    test_ranges = test_ranges.set_index('basin_id').to_dict(orient='index')
except:
    pass

if finetuning_special_ranges:
    ft_train_ranges = ft_train_ranges.set_index('basin_id').to_dict(orient='index')
    ft_validation_ranges = ft_validation_ranges.set_index('basin_id').to_dict(orient='index')
    try:
        ft_test_ranges = ft_test_ranges.set_index('basin_id').to_dict(orient='index')
    except:
        pass

if extra_finetunes:
    ftx_train_ranges = ftx_train_ranges.set_index('basin_id').to_dict(orient='index')
    ftx_validation_ranges = ftx_validation_ranges.set_index('basin_id').to_dict(orient='index')
    try:
        ftx_test_ranges = ftx_test_ranges.set_index('basin_id').to_dict(orient='index')
    except:
        pass

if nhm_run:
    nhm_train_ranges = nhm_train_ranges.set_index('basin_id').to_dict(orient='index')
    nhm_validation_ranges = nhm_validation_ranges.set_index('basin_id').to_dict(orient='index')
    try:
        nhm_test_ranges = nhm_test_ranges.set_index('basin_id').to_dict(orient='index')
    except:
        pass

def listify_range_elems(range_dict):

    for key in range_dict:
        for subkey in range_dict[key]:
            range_dict[key][subkey] = [range_dict[key][subkey]]

    return range_dict

train_ranges = listify_range_elems(train_ranges)
validation_ranges= listify_range_elems(validation_ranges)
try:
    test_ranges = listify_range_elems(test_ranges)
except:
    pass

if finetuning_special_ranges:
    ft_train_ranges = listify_range_elems(ft_train_ranges)
    ft_validation_ranges= listify_range_elems(ft_validation_ranges)
    try:
        ft_test_ranges = listify_range_elems(ft_test_ranges)
    except:
        pass

if extra_finetunes:
    ftx_train_ranges = listify_range_elems(ftx_train_ranges)
    ftx_validation_ranges= listify_range_elems(ftx_validation_ranges)
    try:
        ftx_test_ranges = listify_range_elems(ftx_test_ranges)
    except:
        pass

if nhm_run:
    nhm_train_ranges = listify_range_elems(nhm_train_ranges)
    nhm_validation_ranges= listify_range_elems(nhm_validation_ranges)
    try:
        nhm_test_ranges = listify_range_elems(nhm_test_ranges)
    except:
        pass

with open(Path(curr_config_dir, 'train_ranges.pkl'), 'wb') as pf:
    pickle.dump(train_ranges, pf)
with open(Path(curr_config_dir, 'validation_ranges.pkl'), 'wb') as pf:
    pickle.dump(validation_ranges, pf)
try:
    test_ranges
except NameError:
    pass
else:
    with open(Path(curr_config_dir, 'test_ranges.pkl'), 'wb') as pf:
        pickle.dump(test_ranges, pf)

if finetuning_special_ranges:
    with open(Path(curr_config_dir, 'finetune_train_ranges.pkl'), 'wb') as pf:
        pickle.dump(ft_train_ranges, pf)
    with open(Path(curr_config_dir, 'finetune_validation_ranges.pkl'), 'wb') as pf:
        pickle.dump(ft_validation_ranges, pf)
    try:
        ft_test_ranges
    except NameError:
        pass
    else:
        with open(Path(curr_config_dir, 'finetune_test_ranges.pkl'), 'wb') as pf:
            pickle.dump(ft_test_ranges, pf)

if extra_finetunes:
    with open(Path(curr_config_dir, latest_finetune, 'train_ranges.pkl'), 'wb') as pf:
        pickle.dump(ftx_train_ranges, pf)
    with open(Path(curr_config_dir, latest_finetune, 'validation_ranges.pkl'), 'wb') as pf:
        pickle.dump(ftx_validation_ranges, pf)
    try:
        ftx_test_ranges
    except NameError:
        pass
    else:
        with open(Path(curr_config_dir, latest_finetune, 'test_ranges.pkl'), 'wb') as pf:
            pickle.dump(ftx_test_ranges, pf)

if nhm_run:
    with open(Path(curr_config_dir, 'nhm_train_ranges.pkl'), 'wb') as pf:
        pickle.dump(nhm_train_ranges, pf)
    with open(Path(curr_config_dir, 'nhm_validation_ranges.pkl'), 'wb') as pf:
        pickle.dump(nhm_validation_ranges, pf)
    try:
        nhm_test_ranges
    except NameError:
        pass
    else:
        with open(Path(curr_config_dir, 'nhm_test_ranges.pkl'), 'wb') as pf:
            pickle.dump(nhm_test_ranges, pf)

# with open(Path(curr_config_dir, 'train_ranges.pkl'), 'rb') as pf:
#     zz = pickle.load(pf)
# with open(Path(curr_config_dir, 'test_ranges.pkl'), 'rb') as pf:
#     zz = pickle.load(pf)
# with open('/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs/run19/train_ranges.pkl', 'rb') as pf:
#     zz = pickle.load(pf)
# with open('/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs/run19/validation_ranges.pkl', 'rb') as pf:
#     zz = pickle.load(pf)
# with open('/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs/run51/test_ranges.pkl', 'rb') as pf:
#     zz = pickle.load(pf)
# with open('/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs/runs_582-582/train_ranges.pkl', 'rb') as pf:
#     zz = pickle.load(pf)

if turboroutine or semiturboroutine:
    turbofiles = os.listdir(curr_config_dir)
    turbodirs = [x for x in turbofiles if re.match('^run[0-9]', x)]
    for i in range(len(turbodirs)):
        turbofiles2 = os.listdir(Path(curr_config_dir, turbodirs[i]))
        to_pickle = [x for x in turbofiles2 if re.search('\\.csv$', x)]
        for j in range(len(to_pickle)):
            xx = pd.read_csv(Path(curr_config_dir, turbodirs[i], to_pickle[j]),
                                       dtype={'basin_id': str},
                                       parse_dates=[1,2])
            xx = xx.set_index('basin_id').to_dict(orient='index')
            xx = listify_range_elems(xx)
            with open(Path(curr_config_dir, turbodirs[i],
                           to_pickle[j].replace('.csv', '.pkl')), 'wb') as pf:
                pickle.dump(xx, pf)

## remove CSVs
# os.remove(Path(curr_config_dir, 'train_ranges.csv'))
# os.remove(Path(curr_config_dir, 'validation_ranges.csv'))
# os.remove(Path(curr_config_dir, 'test_ranges.csv'))

def one_off_pickler(x, runlist):

    for r in runlist:
        ranges = pd.read_csv(Path(curr_config_dir, 'run' + r, x + '.csv'),
                                    dtype={'basin_id': str},
                                    parse_dates=[1,2])
        ranges = ranges.set_index('basin_id').to_dict(orient='index')
        ranges = listify_range_elems(ranges)
        with open(Path(curr_config_dir, 'run' + r, x + '.pkl'), 'wb') as pf:
            pickle.dump(ranges, pf)


# runids = [str(x) for x in list(range(2393, 2398))]
# one_off_pickler('finetune_train_ranges', runids)
# one_off_pickler('finetune_validation_ranges', runids)
# one_off_pickler('test_ranges', runids)

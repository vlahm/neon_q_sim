
# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-20

import pickle
from pathlib import Path
import os
import pandas as pd
import re
import numpy as np
from sys import argv

curr_config_dir = Path(argv[1])
#config_dir = Path('in/lstm_configs', r.r2pyenv['confdir'])
#curr_config_rel = r.r2pyenv['runset']

def listify_range_elems(range_dict):
    for key in range_dict:
        for subkey in range_dict[key]:
            range_dict[key][subkey] = [range_dict[key][subkey]]
    return range_dict

#all_runs = os.listdir(config_dir)
#curr_config_dir = Path(config_dir, curr_config_rel)
cfg_fls = os.listdir(curr_config_dir)

#cfg_files2 = os.listdir(Path(curr_config_dir, rundirs[i]))
#to_pickle = [x for x in cfg_files2 if re.search('\\.csv$', x)]
#for j in range(len(to_pickle)):
xx = pd.read_csv(Path(curr_config_dir, 'train_ranges.csv'),
                 dtype={'basin_id': str},
                 parse_dates=[1,2])
xx = xx.set_index('basin_id').to_dict(orient='index')
xx = listify_range_elems(xx)
with open(Path(curr_config_dir, 'train_ranges.pkl'), 'wb') as pf:
    pickle.dump(xx, pf)

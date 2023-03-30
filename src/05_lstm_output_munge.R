# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-03-29

library(readr)
library(stringr)
library(reticulate)

#see step 2 in src/lstm_dungeon/README.txt for installing conda environment
use_condaenv('nh2')
xr <- import("xarray")
pd <- import("pandas")
np <- import("numpy")

source('src/00_helpers.R')

#source the following script (or step through it) if you're not using our pre-bundled data at
# []
# source('src/01_data_retrieval.R')

#source the following scripts (or step through them) if you're not using our pre-bundled data at
# []
# source('src/03_organize_camels_macrosheds_nhm.R')
# source('src/04_run_lstms.R')

#establish i/o directory locations
confdir <- file.path(getwd(), 'in/lstm_configs')
datadir <- file.path(getwd(), 'in/lstm_data')
rundir <- file.path(getwd(), 'out/lstm_runs')

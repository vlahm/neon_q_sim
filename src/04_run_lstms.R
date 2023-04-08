# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-04

library(readr)
library(stringr)
library(reticulate)

#see step 2 in src/lstm_dungeon/README.txt for installing conda environment
use_condaenv('nh2')
xr <- import("xarray")
pd <- import("pandas")
np <- import("numpy")

#pre-bundled in/out data available at: [**]
if(! exists('ts_plot')) source('src/00_helpers.R')
if(! exists('ms_areas')) source('src/01_data_retrieval.R')
if(! dir.exists('in/lstm_data')) source('src/03_organize_camels_macrosheds_nhm.R', local = new.env())

#establish i/o directory locations
confdir <- file.path(getwd(), 'in/lstm_configs')
datadir <- file.path(getwd(), 'in/lstm_data')
rundir <- file.path(getwd(), 'out/lstm_runs')

#insert i/o directories into model/config text files
cfgs <- list.files('.', pattern = '.*\\.yml$', recursive = TRUE, full.names = TRUE)
cfgs <- c(cfgs, list.files('.', pattern = 'pretrained_model_loc',
                           recursive = TRUE, full.names = TRUE))
for(cfg_ in cfgs){

    read_file(cfg_) %>%
        str_replace_all('PLACEHOLDER3', datadir) %>%
        str_replace_all('PLACEHOLDER2', rundir) %>%
        str_replace_all('PLACEHOLDER', confdir) %>%
        write_file(cfg_)
}

#crude way to pass variables into python script from R
r2pyenv_template <- new.env()
r2pyenv_template$wdir <- getwd()
r2pyenv_template$confdir <- confdir
r2pyenv_template$rundir <- rundir

#run generalists (batch 1)
r2pyenv_template$strategy <- 'generalist'
r2pyenv_template$runset <- 'runs_1468-1507'
r2pyenv <- as.list(r2pyenv_template)
py_run_file('src/lstm_dungeon/run_lstms.py')

#run generalists (batch 2)
r2pyenv_template$strategy <- 'generalist'
r2pyenv_template$runset <- 'runs_1508-1520'
r2pyenv <- as.list(r2pyenv_template)
py_run_file('src/lstm_dungeon/run_lstms.py')

#run specialists (batch 1)
r2pyenv_template$strategy <- 'specialist'
r2pyenv_template$runset <- 'runs_2293-2307'
r2pyenv <- as.list(r2pyenv_template)
py_run_file('src/lstm_dungeon/run_lstms.py')

#run specialists (batch 2)
r2pyenv_template$strategy <- 'specialist'
r2pyenv_template$runset <- 'runs_2308-2422'
r2pyenv <- as.list(r2pyenv_template)
py_run_file('src/lstm_dungeon/run_lstms.py')

#run process-guided specialists
r2pyenv_template$strategy <- 'specialist'
r2pyenv_template$runset <- 'runs_2248-2292'
r2pyenv <- as.list(r2pyenv_template)
py_run_file('src/lstm_dungeon/run_lstms.py')

#run ensembles for skilled generalist model (TECR)
r2pyenv_template$strategy <- 'generalist'
r2pyenv_template$runset <- 'runs_2423-2452'
r2pyenv <- as.list(r2pyenv_template)
py_run_file('src/lstm_dungeon/run_lstms.py')

#run ensembles for skilled generalist model (BIGC)
r2pyenv_template$strategy <- 'generalist'
r2pyenv_template$runset <- 'runs_2453-2482'
r2pyenv <- as.list(r2pyenv_template)
py_run_file('src/lstm_dungeon/run_lstms.py')

#run ensembles for skilled generalist model (WALK)
r2pyenv_template$strategy <- 'generalist'
r2pyenv_template$runset <- 'runs_2483-1512'
r2pyenv <- as.list(r2pyenv_template)
py_run_file('src/lstm_dungeon/run_lstms.py')

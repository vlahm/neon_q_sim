# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-03-24

library(readr)
library(reticulate)

#see step 2 in src/lstm_dungeon/README.txt for installing conda environment
use_condaenv('nh2')
xr <- import("xarray")
pd <- import("pandas")
np <- import("numpy")

#establish i/o directory locations
confdir <- file.path(getwd(), 'in/lstm_configs')
datadir <- file.path(getwd(), 'in/lstm_data')
rundir <- file.path(getwd(), 'out/lstm_runs')

#replace i/o directories across model/config text files


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
py_run_file('src/lstm_dungeon/run_lstms.py')

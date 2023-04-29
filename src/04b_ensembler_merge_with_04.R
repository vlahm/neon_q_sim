# library(glue)
# library(ncdf4)
# library(purrr)





#PREPARING LOCAL CONFIG AND RUNS FOR PUBLICATION
#unzip the directories with PLACEHOLDERS,
#ADD run_dir: PLACEHOLDER2 to all configs
#good lord, also conform all finetune_special.txt to finetune.txt and adjust configs accordingly
    #CONFIG: find . -name finetune_special.txt -exec rename -f 's/finetune_special.txt/finetune.txt/' {} \;
    #CONFIG: find . -name '*.yml' | xargs sed -e 's/finetune_special/finetune/g' -i
    #RUNS: find . -name 'config.yml' | xargs sed -e 's/finetune_special/finetune/g' -i
# set save_weights_every: 1 to 10 (this is already fixed locally)
#leave test dirs? probs.
#only modify the new ensemble directories as below, then
#REMOVE ALL non-terminal run files
#add them to the unzipped dirs and re-zip
#MAKE SURE ".", "in/lstm_configs", and "out/lstm_runs" BELOW ARE CHANGED TO THE LOCATION OF THE NEW ENSEMBLE DIRS
ftgrep .yml 'workers: 48'
find . -name '*validation*'
ftgrep .yml run_configs
# find . -name '*.yml' | xargs sed -e 's/num_workers: 48/num_workers: 8/g' -i #shouldn't be necessary, but template if needed

#(DCC config locations)
system(paste0("find in/lstm_configs -name '*.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_data|PLACEHOLDER3|g' -i"))
system(paste0("find in/lstm_configs -name '*.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_runs|PLACEHOLDER2|g' -i"))
system(paste0("find in/lstm_configs -name '*.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_configs|PLACEHOLDER|g' -i"))
# system(paste0("find in/lstm_configs -name 'pretrained_model_loc.txt' | xargs sed -e 's|/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/runs|PLACEHOLDER2|g' -i"))

#(DCC run location)
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_data|PLACEHOLDER3|g' -i"))
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_runs|PLACEHOLDER2|g' -i"))
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/runs|PLACEHOLDER2|g' -i"))
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_configs|PLACEHOLDER|g' -i"))

# system(paste0("find out/ergh -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_data|PLACEHOLDER3|g' -i"))
# system(paste0("find out/ergh -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_runs|PLACEHOLDER2|g' -i"))
# system(paste0("find out/ergh -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/runs|PLACEHOLDER2|g' -i"))
# system(paste0("find out/ergh -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_configs|PLACEHOLDER|g' -i"))

cfgs <- list.files('out/ergh', pattern = 'config.yml', recursive = TRUE, full.names = TRUE)
for(cfg_ in cfgs){
    read_file(cfg_) %>%
        str_replace_all('PLACEHOLDER3', datadir) %>%
        str_replace_all('PLACEHOLDER2', rundir) %>%
        str_replace_all('PLACEHOLDER', confdir) %>%
        write_file(cfg_)
}


#for dcc runs (change cfg parent dir!)
system(paste0("find in/lstm_configs/to_dcc -name '*.yml' | xargs sed -e 's|/home/mike/git/macrosheds/papers/q_sim/in/lstm_data|/hpc/home/mjv22/q_sim/lstm_data|g' -i"))
system(paste0("find in/lstm_configs/to_dcc -name '*.yml' | xargs sed -e 's|/home/mike/git/macrosheds/papers/q_sim/out/lstm_runs|/hpc/home/mjv22/q_sim/lstm_runs|g' -i"))
system(paste0("find in/lstm_configs/to_dcc -name '*.yml' | xargs sed -e 's|/home/mike/git/macrosheds/papers/q_sim/in/lstm_configs|/hpc/home/mjv22/q_sim/lstm_configs|g' -i"))







# system = function(x, ...) system(x, ..., wait = T)

#jk. none of the rest of this file needs to be merged
#but start by running the top of 04_

### COPY FILES, ADJUST TEST PERIODS, RE-PICKLE, RE-EVAL (not for production) ####

RUNTYPE = 'run'#"run" for generalists, 'finetune' for specialists
neon_site = 'TECR' #when you get to MART, need to change name as well
qjj_ = ensembles[[neon_site]][[1]]

#copy configs and rundirs
qjj = paste(range(qjj_), collapse = '-')
# file.copy(paste0('~/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs/runs_', qjj),
#           '~/git/macrosheds/papers/q_sim/in/lstm_configs', recursive = TRUE)
zzj = list.files('~/git/macrosheds/qa_experimentation/imputation/src/nh_methods/runs', pattern = RUNTYPE, full.names = TRUE)
zi_ <- grep(paste(paste0('^', RUNTYPE, qjj_, collapse = '|')), basename(zzj))
rundirs0 <- zzj[zi_]
for(r in rundirs0){
    file.copy(r, 'out/lstm_runs/', recursive = TRUE)
}

#adjust test periods
runset_parent_dir = paste0('in/lstm_configs/runs_', qjj)
rundirs1 = list.files(runset_parent_dir, pattern = 'run', full.names = FALSE)
# zzj = list.files('out/lstm_runs', pattern = RUNTYPE, full.names = TRUE)
# zi_ <- grep(paste(paste0('^', RUNTYPE, qjj_, collapse = '|')), basename(zzj))
# rundirs2 <- zzj[zi_]
for(rd in rundirs1){
    testrng_path = file.path(runset_parent_dir, rd, 'test_ranges.csv')
    read_csv(testrng_path) %>%
        mutate(start_dates = '2015-06-01',
               end_dates = '2024-01-01') %>%
        write_csv(testrng_path)

    # qrq = grep(paste0('^', rd), basename(rundirs2))
    # thisrun = rundirs2[qrq]
    # testrng_path = file.path('out/lstm_runs', thisrun, 'test_ranges.csv')
    # read_csv(testrng_path)
}

#re-pickle
r2pyenv_template <- new.env()
r2pyenv_template$confdir <- confdir
r2pyenv_template$runset <- paste0('runs_', qjj)
r2pyenv <- as.list(r2pyenv_template)
py_run_file('src/lstm_dungeon/dungeon_lvl_2/pickle_train_periods.py')

#replace paths with PLACEHOLDERs
#(new config locations)
system(paste0("find in/lstm_configs -name '*.yml' | xargs sed -e 's|/home/mike/git/macrosheds/papers/q_sim/in/lstm_data|PLACEHOLDER3|g' -i"))
system(paste0("find in/lstm_configs -name '*.yml' | xargs sed -e 's|/home/mike/git/macrosheds/papers/q_sim/out/lstm_runs|PLACEHOLDER2|g' -i"))
system(paste0("find in/lstm_configs -name '*.yml' | xargs sed -e 's|/home/mike/git/macrosheds/papers/q_sim/in/lstm_configs|PLACEHOLDER|g' -i"))
system(paste0("find in/lstm_configs -name 'pretrained_model_loc.txt' | xargs sed -e 's|/home/mike/git/macrosheds/papers/q_sim/out/lstm_runs|PLACEHOLDER2|g' -i"))

#(old config locations)
system(paste0("find in/lstm_configs -name '*.yml' | xargs sed -e 's|/home/mike/git/macrosheds/qa_experimentation/data/CAMELS_macrosheds_combined|PLACEHOLDER3|g' -i"))
system(paste0("find in/lstm_configs -name '*.yml' | xargs sed -e 's|/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/runs|PLACEHOLDER2|g' -i"))
system(paste0("find in/lstm_configs -name '*.yml' | xargs sed -e 's|/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs|PLACEHOLDER|g' -i"))
system(paste0("find in/lstm_configs -name 'pretrained_model_loc.txt' | xargs sed -e 's|/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/runs|PLACEHOLDER2|g' -i"))

#(old run location)
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/home/mike/git/macrosheds/qa_experimentation/data/CAMELS_macrosheds_combined|PLACEHOLDER3|g' -i"))
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/runs|PLACEHOLDER2|g' -i"))
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs|PLACEHOLDER|g' -i"))

#(new run location)
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/home/mike/git/macrosheds/papers/q_sim/in/lstm_data|PLACEHOLDER3|g' -i"))
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/home/mike/git/macrosheds/papers/q_sim/out/lstm_runs|PLACEHOLDER2|g' -i"))
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/home/mike/git/macrosheds/papers/q_sim/in/lstm_configs|PLACEHOLDER|g' -i"))

#(DCC config locations)
system(paste0("find in/lstm_configs -name '*.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_data|PLACEHOLDER3|g' -i"))
system(paste0("find in/lstm_configs -name '*.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_runs|PLACEHOLDER2|g' -i"))
system(paste0("find in/lstm_configs -name '*.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_configs|PLACEHOLDER|g' -i"))
# system(paste0("find in/lstm_configs -name 'pretrained_model_loc.txt' | xargs sed -e 's|/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/runs|PLACEHOLDER2|g' -i"))

#(DCC run location)
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_data|PLACEHOLDER3|g' -i"))
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_runs|PLACEHOLDER2|g' -i"))
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/runs|PLACEHOLDER2|g' -i"))
system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_configs|PLACEHOLDER|g' -i"))

#checks
if(length(suppressWarnings(system('grep -rl --include="*.yml" "qa_experimentation" out/lstm_runs', intern=T)))) stop('oi')
if(length(suppressWarnings(system('grep -rl --include="*.yml" "qa_experimentation" in/lstm_configs', intern=T)))) stop('oi')
if(length(suppressWarnings(system('grep -rl --include="pretrained_model_loc.txt" "qa_experimentation" in/lstm_configs', intern=T)))) stop('oi')
if(length(suppressWarnings(system('grep -rl --include="*.yml" "q_sim" out/lstm_runs', intern=T)))) stop('oi')
if(length(suppressWarnings(system('grep -rl --include="*.yml" "q_sim" in/lstm_configs', intern=T)))) stop('oi')
if(length(suppressWarnings(system('grep -rl --include="pretrained_model_loc.txt" "q_sim" in/lstm_configs', intern=T)))) stop('oi')

#set paths (same as chunk starting line 25 in 04_run_lstms.R)
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

#re-evaluate
r2pyenv_template$wdir <- getwd()
r2pyenv_template$runrange <- as.integer(qjj_)
r2pyenv <- as.list(r2pyenv_template)
py_run_file('src/lstm_dungeon/dungeon_lvl_2/re-evaluate_models.py')


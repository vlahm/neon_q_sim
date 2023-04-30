# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-05

library(tidyverse)
library(reticulate)
library(glue)
library(lubridate)

#see step 2 in src/lstm_dungeon/README.txt for installing conda environment
use_condaenv('nh2')
xr <- import("xarray")
pd <- import("pandas")
np <- import("numpy")

#pre-bundled in/out data available at: [**]
if(! exists('ts_plot')) source('src/00_helpers.R')
if(! exists('ms_areas')) source('src/01_data_retrieval.R')
if(! dir.exists('out/lm_out')) source('src/02_regression.R', local = new.env())
if(! dir.exists('in/lstm_data')) source('src/03_organize_camels_macrosheds_nhm.R', local = new.env())
if(! dir.exists('out/lstm_runs')) stop("you need to run src/04_run_lstms.R. It will take many days unless run on a cluster. Or use our bundled results.")

#might want to run src/07_barplot.R to decide how you want to rank models
ranks <- read_csv('cfg/model_ranks.csv')

# 1. helpers ####

load_q_neon <- function(site){

    q <- read_csv(glue('in/NEON/neon_continuous_Q/{site}.csv')) %>%
        mutate(source = 'NEON') %>%
        select(datetime, discharge_Ls = discharge, discharge_lower95_Ls = discharge_lower,
               discharge_upper95_Ls = discharge_upper, source)

    return(q)
}

load_q_lm <- function(site){

    lm_res <- read_csv('out/lm_out/results.csv') %>%
        filter(site_code == site) %>%
        pull(kge)

    lm_sq_res <- read_csv('out/lm_out_specQ/results_specificq.csv') %>%
        filter(site_code == site) %>%
        pull(kge)

    if(is.na(lm_res)) stop('no lm result for this site')

    q <- read_csv(glue('out/{lm}/predictions/{site}.csv',
                       lm = ifelse(lm_res > lm_sq_res, 'lm_out', 'lm_out_specQ'))) %>%
        mutate(source = ifelse(!!lm_res > !!lm_sq_res, 'Linreg', 'Linreg_scaled'))

    if('date' %in% colnames(q)){
        q$datetime <- as_datetime(paste(q$date, '12:00:00'))
    }

    q <- select(q, datetime, discharge_Ls = Q_predicted,
                discharge_lower95_Ls = Q_pred_int_2.5,
                discharge_upper95_Ls = Q_pred_int_97.5, source)

    return(q)
}

load_q_generalist <- function(site){

    fs <- list.files('out/lstm_runs', full.names = TRUE, recursive = TRUE)

    testsite_by_file <- grep('test_metrics\\.csv$', fs, value = TRUE) %>%
        sort() %>%
        read_lines(skip = 1) %>%
        {sub('(_MANUALQ|,).*', '', .)}

    ensemble_runs_ <- rle2(testsite_by_file) %>%
        filter(values == !!site,
               lengths == 30)

    if(nrow(ensemble_runs_)) stop('more than one ensemble detected for this site')

    ensemble_runs <- fs[ensemble_runs_$starts:ensemble_runs_$stops]



    q <- select(q, datetime, discharge_Ls = Q_predicted,
                discharge_lower95_Ls = Q_pred_int_2.5,
                discharge_upper95_Ls = Q_pred_int_97.5, source)

    return(q)
}

# 2. ####

for(i in seq_len(nrow(ranks))){

    s <- ranks$site[i]
    rankvec <- unlist(ranks[i, 2:4])

    composite <- tibble()
    for(r in rankvec){

        if(r == 'N'){
            composite <- bind_rows(composite, load_q_neon(site = s))
        } else if(r == 'L'){
            composite <- bind_rows(composite, load_q_lm(site = s))
        } else if(r == 'G'){
            HERE: was working on load_q_generalist, but actually need to adapt collect_best_q_predictions.R
        } else if(r == 'S'){

        } else if(r == 'P'){

        } else {
            stop('model type ', r, ' not recognized')
        }
    }

    composite %>%
        complete() %>%
        na_seadec(up to 15 m) %>%
        write_csv()
}

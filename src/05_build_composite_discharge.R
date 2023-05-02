# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-05

library(tidyverse)
library(reticulate)
library(glue)
library(lubridate)
library(zoo)
library(xts)
library(imputeTS)

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
results <- read_csv('out/lstm_out/results.csv')

dir.create('figs/smooth_plots', showWarnings = FALSE)

# 1. helpers ####

load_q_neon <- function(site, smooth_plot = FALSE){

    window_size <- 15 #keep it odd
    transient <- window_size %/% 2
    q <- read_csv(glue('in/NEON/neon_continuous_Q/{site}.csv'))
    # q_head <- head(q, transient) %>% select(datetime, discharge)
    # q_tail <- tail(q, transient) %>% select(datetime, discharge)
    q_etc <- rename(q, discharge_orig = discharge) %>% mutate(indicator = 1)
    q <- select(q, datetime, discharge)

    #smooth out bayesian obs error from neon product
    if(site != 'TOMB'){

        q <- complete(q, datetime = seq(min(datetime), max(datetime), by = '1 min'))
        q_dt <- q$datetime
        # q_orig <- q$discharge

        # q$discharge_loess <- q$discharge_lower_loess <- q$discharge_upper_loess <- NA_real_
        # lvals <- ! is.na(q$discharge)
        # lind <- 1:nrow(q)
        #
        # q$discharge_loess[lvals] <- predict(loess(q$discharge ~ I(lind), span = 0.0005, degree = 1))

        q <- xts(select(q, -datetime), q$datetime, tzone = 'UTC')
        q <- rollmean(q, window_size, fill = NA, align = 'center', na.rm = TRUE)
        q <- rollmean(q, window_size, fill = NA, align = 'center', na.rm = TRUE)

        # q <- restore_transient(q, q_head, q_tail, transient)

        q <- suppressWarnings(q) %>%
            as_tibble() %>%
            mutate(datetime = !!q_dt) %>%
            left_join(q_etc, by = 'datetime') %>%
            filter(! is.na(indicator))

        if(smooth_plot){
            require(dygraphs)
            require(htmlwidgets)
            dygraphs::dygraph(xts(x = select(q, discharge_orig, discharge),
            # dygraphs::dygraph(xts(x = select(q, discharge, discharge_loess),
                                  order.by = q$datetime)) %>%
                dyRangeSelector() %>%
                saveWidget(glue('figs/smooth_plots/{site}.html'))
        }

        # q$discharge <- q$discharge_loess
        # q$discharge_lower <- q$discharge_lower_loess
        # q$discharge_upper <- q$discharge_upper_loess
    }

    q <- q %>%
        slice((transient + 1):(nrow(q) - transient)) %>%
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

load_q_lstm <- function(site){

    q <- read_csv(glue('out/lstm_out/predictions/{site}.csv'))

    q$source <- filter(results, site_code == !!site)$strategy
    q$datetime <- as_datetime(paste(q$date, '12:00:00'))

    q <- select(q, datetime, discharge_Ls = Q_predicted,
                discharge_lower95_Ls = Q_pred_int_2.5,
                discharge_upper95_Ls = Q_pred_int_97.5, source)

    return(q)
}

# 2. ####

Show all missing for PRIN, OKSR,
use GUIL N only from good years
downsample neon to 5m
is this set up to deal with como?

for(i in seq_len(nrow(ranks))){

    s <- ranks$site[i]
    rankvec <- unlist(ranks[i, 2:4])

    composite <- tibble()
    for(r in rankvec){

        if(r == 'N'){
            composite <- bind_rows(composite, load_q_neon(site = s, smooth_plot = TRUE))
        } else if(r == 'L'){
            composite <- bind_rows(composite, load_q_lm(site = s))
        } else if(r %in% c('G', 'S', 'PG', 'PS')){
            composite <- bind_rows(composite, load_q_lstm(site = s))
        } else {
            stop('model type ', r, ' not recognized')
        }

        # mode_interval_dt(composite$datetime)
    }
    readLines(n=1)

    # composite %>%
    #     complete() %>%
    #     na_seadec(up to 15 m) %>%
    #     write_csv()
}

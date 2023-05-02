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
    q_bounds <- select(q, -site_code, -discharge)
    q_orig <- select(q, datetime, discharge_orig = discharge)
    # q_etc <- rename(q, discharge_orig = discharge) %>% mutate(indicator = 1)
    q <- select(q, datetime, discharge)

    #smooth out bayesian obs error from neon product
    if(site != 'TOMB'){

        #fill out missing data (make explicit) to one minute interval
        q <- complete(q, datetime = seq(min(datetime), max(datetime), by = '1 min'))

        #pad head and tail for moving average window
        buffer1 <- tibble(datetime = seq(q$datetime[1] - (60 * transient),
                                         q$datetime[1] - 60,
                                         by = '1 min'))
        buffer2 <- tibble(datetime = seq(q$datetime[nrow(q)] + 60,
                                         q$datetime[nrow(q)] + (60 * transient),
                                         by = '1 min'))
        q <- bind_rows(buffer1, q, buffer2)

        q_dt <- q$datetime
        q_na <- is.na(q$discharge)

        #traingular moving average
        q <- xts(q$discharge, q$datetime, tzone = 'UTC')
        q <- rollmean(q, window_size, fill = NA, align = 'center', na.rm = TRUE)
        q <- rollmean(q, window_size, fill = NA, align = 'center', na.rm = TRUE)
        q <- suppressWarnings(as.vector(q))

        #restore missing values
        q[q_na] <- NA_real_

        #interpolate decomposed series with daily periodicity.
        #need to interpolate uncertainty too, so it isn't lost during downsampling
        q <- left_join(tibble(discharge = q, datetime = q_dt), q_bounds, by = 'datetime') %>%
            slice((transient + 1):(length(q) - transient)) %>%
            mutate(across(-datetime, ~na_seadec(ts(., deltat = 1/1440), maxgap = 14)))
        # q <- na_seadec(ts(q, deltat = 1/1440), maxgap = 14)

        #downsample to even 5 minute interval. priority breakdown:
        #   1. preserve observations originally on minutes divisible by 5
        #   2. accept newly interpolated values where the original interval was
        #      2-15 minutes and/or offset from a natural 5-minute sequence
        #   3. shift original values by up to 2 minutes if interpolation was
        #      not performed, e.g. where the sampling interval is 16-60 minutes

        # q <- tibble(discharge = q, datetime = q_dt)
        # q <- left_join(q, q_bounds, by = 'datetime')
        q$m5 <- FALSE
        q$m5[minute(q$datetime) %% 5 == 0] <- TRUE
        q$datetime <- round_date(q$datetime, unit = '5 min')
        q <- group_by(q, datetime) %>%
            slice(if(any(m5)) which(m5)[1] else 1) %>%
            ungroup() %>%
            select(-m5)

        # qt = tibble(datetime = q_dt[1:10], discharge = 1:10, m5 = q$m5[1:10]) %>%
        #     slice(c(1, 6:10)) %>%
        #     mutate(datetime = round_date(datetime, unit = '5 min'))
        # group_by(qt, datetime) %>%
        #     slice(if(any(m5)) which(m5)[1] else 1)

        #downsample to 5 minutes
        # q <- filter(q, minute(datetime) %% 5 == 0)

        #
        # q <- as_tibble(q) %>%
        #     rename(discharge = x) %>%
        #     mutate(datetime = !!q_dt) %>%
        #     left_join(q_etc, by = 'datetime')
            # filter(! is.na(indicator))
            # mutate(discharge = ifelse(is.na(discharge_orig), NA_real_, discharge))

        # zz[is.na(zz$discharge) & !is.na(zz$discharge_orig),]
        # zz[is.na(zz$discharge_orig) & !is.na(zz$discharge),]
        # q[is.na(q$discharge) & !is.na(q$discharge_orig),]
        # q[is.na(q$discharge_orig) & !is.na(q$discharge),]

        if(smooth_plot){

            require(dygraphs)
            require(htmlwidgets)

            zz = full_join(q, q_orig, by = 'datetime')
            dygraphs::dygraph(xts(x = select(zz, discharge_orig, discharge_lower, discharge_upper, discharge),
                                  order.by = zz$datetime)) %>%
                dyRangeSelector() %>%
                saveWidget(glue('figs/smooth_plots/{site}.html'))
        }
    }

    q <- q %>%
        # slice((transient + 1):(nrow(q) - transient)) %>%
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
# downsample neon to 5m
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

        #HERE: need a stronger feel for the intervals present. where do we see 30, 60?
        # also need to know if usgs/hj/niwot all evenly on 5min
        # first, replace NAs after running rollmean
        then, interp N (maxgap 14) and downscale to 5
        do not interpolate uncertainty
        also interp L to 5 (maxgap 14); will require complete to 1 if usgs/hj/niwot not evenly on 5min
        then, just insert r-2 wherever there are NAs in r-1
        ...if there are 30s and/or 60s in an otherwise 1-15 series, consider interping to 5
        for TOMB, leave it at 60
        for COMO?

        # mode_interval_dt(composite$datetime)
    }
    readLines(n=1)

    # composite %>%
    #     complete() %>%
    #     na_seadec(up to 15 m) %>%
    #     write_csv()
}

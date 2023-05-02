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

prepare_q_neon <- function(site, smooth_plot = FALSE){

    window_size <- 15 #keep it odd
    transient <- window_size %/% 2
    q <- read_csv(glue('in/NEON/neon_continuous_Q/{site}.csv'))
    q_bounds <- select(q, -site_code, -discharge)
    q_orig <- select(q, datetime, discharge_orig = discharge)
    q <- select(q, datetime, discharge)

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

        #traingular moving average to smooth over bayesian obs error
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

        #downsample to even 5 minute interval. priority breakdown:
        #   1. preserve observations originally on minutes divisible by 5
        #   2. accept newly interpolated values where the original interval was
        #      2-15 minutes and/or offset from a natural 5-minute sequence
        #   3. shift original values by up to 2 minutes if interpolation was
        #      not performed, e.g. where the sampling interval is 16-60 minutes
        q$m5 <- FALSE
        q$m5[minute(q$datetime) %% 5 == 0] <- TRUE
        q$datetime <- round_date(q$datetime, unit = '5 min')
        q <- group_by(q, datetime) %>%
            slice(if(any(m5)) which(m5)[1] else if(any(! is.na(discharge))) which(! is.na(discharge))[1] else 1) %>%
            ungroup() %>%
            select(-m5)

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
        mutate(source = 'NEON') %>%
        select(datetime, discharge_Ls = discharge, discharge_lower95_Ls = discharge_lower,
               discharge_upper95_Ls = discharge_upper, source)

    return(q)
}

prepare_q_lm <- function(site){

    lm_res <- read_csv('out/lm_out/results.csv') %>%
        filter(site_code == site) %>%
        pull(kge)

    lm_sq_res <- read_csv('out/lm_out_specQ/results_specificq.csv') %>%
        filter(site_code == site) %>%
        pull(kge)

    if(is.na(lm_res)) stop('no lm result for this site')

    strategy <- ifelse(lm_res > lm_sq_res, 'lm_out', 'lm_out_specQ')

    q <- read_csv(glue('out/{strategy}/predictions/{site}.csv'))

    if('date' %in% colnames(q)){
        q$datetime <- as_datetime(paste(q$date, '12:00:00'))
    }

    q <- select(q, datetime, discharge_Ls = Q_predicted,
                discharge_lower95_Ls = Q_pred_int_2.5,
                discharge_upper95_Ls = Q_pred_int_97.5)

    # plot(q$datetime, q$discharge_Ls, type = 'l', ylim = c(0, 1000))
    # points(q$datetime[!q$m5], q$discharge_Ls[!q$m5], col = 'red')

    q <- complete(q, datetime = seq(min(datetime), max(datetime), by = '5 min'))
    q$m5 <- FALSE
    q$m5[minute(q$datetime) %% 5 == 0] <- TRUE
    q$datetime <- round_date(q$datetime, unit = '5 min')
    q <- group_by(q, datetime) %>%
        slice(if(any(m5)) which(m5)[1] else which(! is.na(discharge_Ls))[1]) %>%
        ungroup() %>%
        select(-m5)

    q$discharge_Ls <- na_seadec(ts(q$discharge_Ls, deltat = 1/288), maxgap = 2)

    q$source <- ifelse(strategy == 'lm_out', 'Linreg', 'Linreg_scaled')

    return(q)
}

prepare_q_lstm <- function(site){

    q <- read_csv(glue('out/lstm_out/predictions/{site}.csv'))

    q$source <- filter(results, site_code == !!site)$strategy
    q$datetime <- as_datetime(paste(q$date, '12:00:00'))

    q <- select(q, datetime, discharge_Ls = Q_predicted,
                discharge_lower95_Ls = Q_pred_int_2.5,
                discharge_upper95_Ls = Q_pred_int_97.5, source)

    return(q)
}

# 2. ####

# for(i in 2:nrow(ranks)){
#     s <- ranks$site[i]
#     zz = try(load_q_neon(site = s, smooth_plot = TRUE))
#     if(! inherits(zz, 'try-error')) feather::write_feather(zz, glue('~/Desktop/q_sim_junk/clean_N/{s}.feather'))
# }
# for(i in 2:nrow(ranks)){
#     s <- ranks$site[i]
#     qq = q_eval %>%
#         filter(site == !!s) %>%
#         distinct(final_qual)
#     print(s); print(qq)
# }
for(i in seq_len(nrow(ranks))){

    s <- ranks$site[i]
    rankvec <- unlist(ranks[i, 2:4])

    composite <- tibble()
    for(r in rankvec){

        if(r == 'N'){
            # composite <- bind_rows(composite, prepare_q_neon(site = s, smooth_plot = TRUE))
            composite <- feather::read_feather(glue('~/Desktop/q_sim_junk/clean_N/{s}.feather'))

            good_months <- q_eval %>%
                filter(site == !!s, final_qual %in% c('Tier1', 'Tier2')) %>%
                select(year, month) %>%
                mutate(good = 1)

            composite <- composite %>%
                mutate(year = year(datetime),
                       month = month(datetime)) %>%
                left_join(good_months, by = c('year', 'month'))

            neon_pass <- composite %>%
                filter(! is.na(good)) %>%
                select(-year, -month, -good)

            good_wys <- read_csv(glue('out/neon_wateryear_assess/{s}.csv')) %>%
                filter(nse >= 0.5, kge >= 0.5)

            neon_backup <- composite %>%
                filter(is.na(good)) %>%
                mutate(wateryr = year,
                       wateryr = ifelse(month >= 10, wateryr + 1, wateryr)) %>%
                filter(wateryr %in% good_wys$wateryr) %>%
                select(-year, -month, -good, -wateryr)

        } else if(r == 'L'){
            composite <- bind_rows(composite, prepare_q_lm(site = s))
        } else if(r %in% c('G', 'S', 'PG', 'PS')){
            composite <- bind_rows(composite, prepare_q_lstm(site = s))
        } else {
            stop('model type ', r, ' not recognized')
        }

        if(rankvec[1] == 'N')

        #HERE: need a stronger feel for the intervals present. where do we see 30, 60?
        # also need to know if usgs/hj/niwot all evenly on 5min
        # first, replace NAs after running rollmean
        # then, interp N (maxgap 14) and downscale to 5
        # do not interpolate uncertainty
        # also interp L to 5 (maxgap 14); will require complete to 1 if usgs/hj/niwot not evenly on 5min
        then, just insert r-2 wherever there are NAs in r-1
        ...if there are 30s and/or 60s in an otherwise 1-15 series, consider interping to 5
        for TOMB, leave it at 60
        for COMO?
        Show all missing for PRIN, OKSR lm
        use GUIL N only from good years
        NEED TO RERUN NEON BLOCK (and stop reading feathers in r == N
        PUT LOW QUAL N BACK IF THERE ARE STILL GAPS?
            ONLY IF IT PASSES neon_wateryear_assess

        # mode_interval_dt(composite$datetime)
    }
    readLines(n=1)

    # composite %>%
    #     complete() %>%
    #     na_seadec(up to 15 m) %>%
    #     write_csv()
}

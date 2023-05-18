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
library(dygraphs)
library(htmlwidgets)

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
dir.create('figs/composite_plots', showWarnings = FALSE)
dir.create('out/composite_series', showWarnings = FALSE)

# 1. helpers ####

prepare_q_neon <- function(site, smooth_plot = FALSE){

    window_size <- 15 #keep it odd
    q <- read_csv(glue('in/NEON/neon_continuous_Q/{site}.csv'))

    if(site != 'TOMB'){

        transient <- window_size %/% 2
        q_bounds <- select(q, -site_code, -discharge)
        q_orig <- select(q, datetime, discharge_orig = discharge)
        q <- select(q, datetime, discharge)

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
            mutate(across(-datetime, ~na_interpolation(., , maxgap = 14)))
            # mutate(across(-datetime, ~na_seadec(ts(., deltat = 1/1440), maxgap = 14)))

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

        if(any(duplicated(q$datetime))) warning('dupes in ', site)

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

    if(any(duplicated(q$datetime))) warning('dupes in ', site)

    # q$discharge_Ls <- na_seadec(ts(q$discharge_Ls, deltat = 1/288), maxgap = 2)
    q$discharge_Ls <- na_interpolation(q$discharge_Ls, maxgap = 2)

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

# 2. build composite series ####

for(i in seq_len(nrow(ranks))){

    s <- ranks$site[i]
    rankvec <- unlist(ranks[i, 2:4])

    composite <- tibble(datetime = seq(as.POSIXct('2014-01-01', tz = 'UTC'),
                                       as.POSIXct('2024-01-01', tz = 'UTC'),
                                       by = '5 min'),
                        src_n = NA_integer_)

    pgdl <- unname(grep('G|S|PG|PS', rankvec, value = TRUE))
    rankvec <- grep('G|S|PG|PS', rankvec, value = TRUE, invert = TRUE)
    rankvec <- rankvec[! is.na(rankvec)]

    for(j in seq_along(rankvec)){

        r <- rankvec[j]
        if(r == 'N'){

            comp_ <- prepare_q_neon(site = s, smooth_plot = TRUE)

            good_months <- q_eval %>%
                filter(site == !!s, final_qual %in% c('Tier1', 'Tier2')) %>%
                select(year, month) %>%
                distinct() %>%
                mutate(good = 1)

            comp_ <- mutate(comp_, year = year(datetime), month = month(datetime))

            if(s == 'TOMB'){
                comp_$good <- 1
            } else {
                comp_ <- left_join(comp_, good_months, by = c('year', 'month'))
            }

            good_wys <- read_csv(glue('out/neon_wateryear_assess/{s}.csv')) %>%
                filter(nse >= 0.5, kge >= 0.5)

            neon_backup <- comp_ %>%
                filter(is.na(good)) %>%
                mutate(wateryr = year,
                       wateryr = ifelse(month >= 10, wateryr + 1, wateryr)) %>%
                filter(wateryr %in% good_wys$wateryr) %>%
                select(-year, -month, -good, -wateryr) %>%
                mutate(src_n = 3)

            comp_ <- comp_ %>%
                filter(! is.na(good)) %>%
                select(-year, -month, -good) %>%
                mutate(src_n = j)

        } else if(r == 'L'){

            comp_ <- prepare_q_lm(site = s) %>%
                filter(! is.na(discharge_Ls)) %>%
                mutate(src_n = j)

        } else {
            stop('model type ', r, ' not recognized')
        }

        if(j == 1){
            composite[composite$datetime %in% comp_$datetime, 'src_n'] <- 1
            r1 <- comp_
        } else {
            popul_inds <- get_populatable_indices(composite$src_n, mingap = 144)
            comp_ <- filter(comp_, datetime %in% composite$datetime[popul_inds])
            composite[popul_inds & composite$datetime %in% comp_$datetime, 'src_n'] <- 2
            r2 <- comp_
        }
    }

    if('N' %in% rankvec){
        popul_inds <- get_populatable_indices(composite$src_n, mingap = 144)
        neon_backup <- filter(neon_backup, datetime %in% composite$datetime[popul_inds])
        composite[popul_inds & composite$datetime %in% neon_backup$datetime, 'src_n'] <- 3
    }

    if(length(pgdl)){

        comp_ <- prepare_q_lstm(site = s) %>%
            filter(! is.na(discharge_Ls)) %>%
            mutate(src_n = 4)

        popul_inds <- get_populatable_indices(composite$src_n, mingap = 144)
        comp_ <- filter(comp_, datetime %in% composite$datetime[popul_inds])
        composite[popul_inds & composite$datetime %in% comp_$datetime, 'src_n'] <- 4
    }

    d <- r1
    if(exists('r2')) d <- bind_rows(d, r2)
    if(exists('neon_backup')) d <- bind_rows(d, neon_backup)
    if(length(pgdl)) d <- bind_rows(d, comp_)

    #remove leading/trailing unpopulated rows
    got <- which(! is.na(composite$src_n))
    composite <- composite[do.call(`:`, as.list(range(got))), ]

    composite <- left_join(composite, d, by = c('datetime', 'src_n')) %>%
        select(-src_n) %>%
        arrange(datetime)

    write_csv(composite, glue('out/composite_series/{s}.csv'))

    suppressWarnings(rm(list = c('r1', 'r2', 'neon_backup', 'pgdl')))

    #files too big and render too slowly if bounds included in plot
    q_plot <- read_csv(glue('in/NEON/neon_continuous_Q/{s}.csv')) %>%
        select(datetime, discharge_original = discharge) %>%
               # lower95_original = discharge_lower,
               # upper95_original = discharge_upper) %>%
        full_join(composite, by = 'datetime') %>%
        full_join(select(read_csv(glue('in/NEON/neon_field_Q/{s}.csv')),
                         datetime, discharge_field_measured = discharge),
                  by = 'datetime') %>%
        select(-discharge_lower95_Ls, -discharge_upper95_Ls) %>%
        arrange(datetime)
        # rename(lower95_composite = discharge_lower95_Ls,
        #        upper95_composite = discharge_upper95_Ls)

    sources <- na.omit(unique(q_plot$source))
    for(src in sources){
        newcol <- paste0('discharge_', src)
        if(src == 'NEON') newcol <- paste0(newcol, '_mod')
        popul_inds <- !is.na(q_plot$source) & q_plot$source == src
        q_plot[[newcol]] <- NA_real_
        q_plot[popul_inds, newcol] <- q_plot$discharge_Ls[popul_inds]
    }

    q_plot <- select(q_plot, -source, -discharge_Ls)

    dygraphs::dygraph(xts(x = select(q_plot, -datetime),
                          order.by = q_plot$datetime)) %>%
        dySeries('discharge_field_measured', pointSize = 4, pointShape = 'ex', strokeWidth = 2) %>%
        dyRangeSelector() %>%
        saveWidget(glue('figs/composite_plots/{s}.html'))
}

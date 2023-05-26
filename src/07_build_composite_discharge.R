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

#pre-bundled in/out data available at: 10.6084/m9.figshare.c.6488065
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

# 1. build composite series ####

for(i in seq_len(nrow(ranks))){

    s <- ranks$site[i]
    rankvec <- unlist(ranks[i, 2:4])

    composite <- tibble(datetime = seq(as.POSIXct('2014-01-01', tz = 'UTC'),
                                       as.POSIXct('2021-09-30', tz = 'UTC'),
                                       by = '5 min'),
                        src_n = NA_integer_)

    pgdl <- unname(grep('G|S|PG|PS', rankvec, value = TRUE))
    rankvec <- grep('G|S|PG|PS', rankvec, value = TRUE, invert = TRUE)
    rankvec <- rankvec[! is.na(rankvec)]

    for(j in seq_along(rankvec)){

        r <- rankvec[j]
        if(r == 'N'){

            comp_ <- prepare_q_neon(site = s, smooth_plot = FALSE)
            comp_ <- filter(comp_, ! is.na(discharge_Ls))

            good_months <- q_eval %>%
                filter(site == !!s, final_qual %in% c('Tier1', 'Tier2')) %>%
                select(year, month) %>%
                distinct() %>%
                mutate(good = 1)

            comp_ <- mutate(comp_, year = year(datetime), month = month(datetime))

            if(s %in% c('TOMB', 'ARIK')){
                comp_$good <- 1
            } else {
                comp_ <- left_join(comp_, good_months, by = c('year', 'month'))
            }

            good_wys <- read_csv(glue('out/neon_wateryear_assess/{s}.csv')) %>%
                filter(nse >= 0.4, kge >= 0.5)

            neon_backup <- comp_ %>%
                filter(is.na(good),
                       ! is.na(discharge_Ls)) %>%
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

# 2. supplemental table of best models for each site ####

lstm_codes <- c('S', 'G', 'PG', 'PS')

lm_incl <- ranks %>%
    filter(rank1 == 'L' | rank2 == 'L') %>%
    pull(site)
lstm_incl <- ranks %>%
    filter(rank2 %in% lstm_codes | rank3 %in% lstm_codes) %>%
    pull(site)

res_lm <- read_csv('out/lm_out/results.csv') %>%
    select(site_code, kge, nse, method) %>%
    mutate(method = paste0(method, '_abs'))
res_lms <- read_csv('out/lm_out_specQ/results_specificq.csv') %>%
    select(site_code, kge, nse, method) %>%
    mutate(method = paste0(method, '_spec'))
res_lstm <- read_csv('out/lstm_out/results.csv') %>%
    select(site_code, kge = KGE, nse = NSE, method = strategy) %>%
    filter(site_code %in% lstm_incl)

bind_rows(res_lm, res_lms) %>%
    group_by(site_code) %>%
    filter(kge == max(kge)) %>%
    arrange(desc(kge)) %>%
    filter(site_code %in% lm_incl) %>%
    full_join(res_lstm, by = 'site_code', suffix = c('_linreg', '_lstm')) %>%
    bind_rows(tibble(site_code = c('PRIN', 'OKSR'))) %>%
    write_csv('out/models_used_to_build_composite_series.csv')


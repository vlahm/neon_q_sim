# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-05-16

library(tidyverse)
library(neonUtilities)
library(viridis)
library(lubridate)

options(readr.show_progress = FALSE,
        readr.show_col_types = FALSE,
        timeout = 3000)

#pre-bundled in/out data available at: [**]
if(! exists('ts_plot')) source('src/00_helpers.R')
if(! exists('ms_areas')) source('src/01_data_retrieval.R')
if(! dir.exists('out/lm_out')) source('src/02_regression.R', local = new.env())
if(! dir.exists('in/lstm_data')) source('src/03_organize_camels_macrosheds_nhm.R', local = new.env())
if(! dir.exists('out/lstm_runs')) stop("you need to run src/04_run_lstms.R. It will take many days unless run on a cluster. Or use our bundled results.")


## 1. data prep ####

neon_sites <- read_csv('in/NEON/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    pull(SiteID)

if(! length(list.files('in/NEON/neon_continuous_Q_withflags/'))){
    get_neon_inst_discharge(neon_sites, clean_only = FALSE)
}

ranks <- read_csv('cfg/model_ranks.csv')
results <- read_csv('out/lstm_out/results.csv')

#neon-official continuous Q data
plotd <- list()
for(s in neon_sites){
    plotd[[s]] <- read_csv(paste0('in/NEON/neon_continuous_Q_withflags/', s, '.csv'))
}

#linear regression predictions
plotfill <- list()
for(s in neon_sites){
    plotfill[[s]] <- try(read_csv(paste0('out/lm_out/predictions/', s, '.csv')),
                         silent = TRUE)
}

## 2. plot ####

neon_sites <- sort(neon_sites)

png(width = 8, height = 12, units = 'in', type = 'cairo', res = 300,
    filename = 'figs/gaps.png')
par(mfrow = c(27, 1), mar = c(1, 3, 0, 0), oma = c(2, 2, 2, 10))
for(i in seq_along(neon_sites)){

    s <- neon_sites[i]

    d <- plotd[[s]] %>%
        filter(minute(datetime) %% 15 == 0) %>%
        mutate(discharge = neglog(discharge)) %>%
        tidyr::complete(datetime = seq(min(datetime), max(datetime), by = '15 min')) %>%
        select(datetime, discharge)

    if(! inherits(plotfill[[s]], 'try-error')){

        cmplt <- ifelse(s == 'COMO', '1 day', '15 min')

        dfill <- plotfill[[s]] %>%
            mutate(across(any_of('date'), as_datetime)) %>%
            rename_with(function(x) sub('^date$', 'datetime', x)) %>%
            filter(minute(datetime) %% 15 == 0) %>%
            select(datetime, Q_predicted) %>%
            tidyr::complete(datetime = seq(min(datetime), max(datetime), by = cmplt)) %>%
            mutate(Q_predicted = neglog(Q_predicted)) %>%
            filter(datetime > min(d$datetime),
                   datetime < max(d$datetime))
    } else {
        dfill <- NULL
    }

    ensembled <- any(grepl('G|S|PG|PS', unlist(filter(ranks, site == !!s)[, 2:4])))
    if(ensembled){

        dlstm <- read_csv(glue('out/lstm_out/predictions/{s}.csv')) %>%
            mutate(datetime = as_datetime(date)) %>%
            mutate(Q_predicted = neglog(Q_predicted)) %>%
            filter(datetime > min(d$datetime),
                   datetime < max(d$datetime)) %>%
            select(-date)
    }

    plot(d$datetime, d$discharge, type = 'n',
         ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n')

    if(! is.null(dfill)){

        dboth <- full_join(dfill, d, by = 'datetime') %>%
            mutate(still_missing = is.na(discharge) & is.na(Q_predicted)) %>%
            arrange(datetime)
        dboth$still_missing[dboth$still_missing] <- NA

        polygon_with_gaps2(dboth, 'Q_predicted', 0, 100, '#FFBDF9')
        polygon_with_gaps2(dboth, 'discharge', 0, 100, '#bdd9ff')

    } else {

        polygon_with_gaps2(d, 'discharge', 0, 100, '#FFBDF9', invert = TRUE)
        polygon_with_gaps2(d, 'discharge', 0, 100, '#c09eff')

        if(ensembled){
            polygon_with_gaps2(mutate(d, allna = NA_real_),
                               'allna', 0, 100, 'white', hashed = TRUE)
        }
    }

    if(! is.null(dfill)){

        polygon_with_gaps2(dboth, 'still_missing', 0, 100, '#c09eff')

        if(ensembled){
            polygon_with_gaps2(dboth, 'Q_predicted', 0, 100, 'white', hashed = TRUE)
        }

        if(s == 'COMO'){
            dboth %>%
                group_by(as.Date(datetime)) %>%
                mutate(Q_predicted = if(any(! is.na(Q_predicted))) 1 else NA_real_) %>%
                ungroup() %>%
                polygon_with_gaps2('Q_predicted', 0, 100, 'white', hashed = TRUE, invert = TRUE)
        }
    }

    if(s == 'TOMB'){

        dtomb <- filter(d, ! is.na(discharge))
        lines(dtomb$datetime, dtomb$discharge, col = 'gray70')

    } else {

        lines(d$datetime, d$discharge, col = 'gray70')

        if(ensembled){
            dmask <- if(is.null(dfill)) NULL else rename(dfill, discharge = Q_predicted)
            plot_gap_points(d = d, dfill = dlstm, mingap = 6, dmask = dmask)
        }
    }

    if(! is.null(dfill)){
        plot_gap_points(d = d, dfill = dfill, mingap = 6)
    }

    if(i == 1){

        legend(x = as.POSIXct('2021-09-19'), y = -80, bty = 'n', border = FALSE, cex = 1.5,
               legend = c('NEON\ndischarge', 'Reconstruction\ngapfill'),
                          # 'Reconstruction\ndaily'),
               col = c('gray70', 'black'), lty = 1, lwd = 2, xpd = NA,
               seg.len = 1.2, y.intersp = 2)

        legend(x = as.POSIXct('2021-09-30'), y = -103, bty = 'n', border = FALSE, cex = 1.5,
               legend = c('NEON\nmissing', 'Reconstruction\nmissing',
                          'Both missing', 'Reconstruction\ndaily',
                          'Daily only'),
               fill = c('#bdd9ff', '#FFBDF9', '#c09eff', '#FFBDF9', '#c09eff'),
               density = c(NA, NA, NA, 30, 30), angle = 45,
               xpd = NA, x.intersp = 1.1, y.intersp = 2)
    }

    mtext(s, side = 2, las = 1, line = -1.2)
}

mtext(expression('Log Discharge'), side = 2, outer = TRUE, cex = 1.2,
      line = -0.5)

dev.off()


## 3. gap stats ####

ranks <- read_csv('cfg/model_ranks.csv')

gap_stats <- matrix(
    NA_real_,
    nrow = length(neon_sites),
    ncol = 26,
    dimnames = list(
        neon_sites,
        c('n_tot', 'n_gap', 'n_fill', 'n_nodonor', 'n_noreconst', 'n_doublegap', 'n_triplegap',
          'pct_gap', 'pct_fill', 'pct_nodonor', 'pct_noreconst', 'pct_doublegap', 'pct_triplegap',
          'neon_intvl_orig', 'est_intvl_orig', 'cast_intvl',
          'nday_tot', 'nday_gap', 'nday_fill', 'nday_nodonor', 'nday_noreconst', 'nday_doublegap',
          'nday_triplegap', 'n_fill_lstm', 'nday_fill_lstm', 'pct_fill_lstm')
    )
)

for(i in seq_along(neon_sites)){

    s <- neon_sites[i]
    print(paste(i, s))

    d <- select(plotd[[s]], -site_code)

    gap_stats[i, 'neon_intvl_orig'] <- mode_interval_dt(d$datetime)
    synch_intvl <- ifelse(unname(gap_stats[i, 'neon_intvl_orig']) == 60, 60, 15)
    gap_stats[i, 'cast_intvl'] <- synch_intvl
    #               m/d  h/d
    nsamp_per_day <- 60 * 24 / synch_intvl

    d <- d %>%
        filter(minute(datetime) %% synch_intvl == 0) %>%
        tidyr::complete(datetime = seq(min(datetime), max(datetime),
                                       by = paste(synch_intvl, 'min')))

    gap_stats[i, 'n_tot'] <- nrow(d)
    gap_stats[i, 'n_gap'] <- sum(is.na(d$discharge))
    gap_stats[i, 'nday_tot'] <- round(gap_stats[i, 'n_tot'] / nsamp_per_day, 2)
    gap_stats[i, 'nday_gap'] <- round(gap_stats[i, 'n_gap'] / nsamp_per_day, 2)

    if(! inherits(plotfill[[s]], 'try-error')){

        dfill <- plotfill[[s]] %>%
            mutate(across(any_of('date'), as_datetime)) %>%
            rename_with(function(x) sub('^date$', 'datetime', x)) %>%
            select(datetime, Q_predicted)

        gap_stats[i, 'est_intvl_orig'] <- mode_interval_dt(dfill$datetime)

        dfill <- dfill %>%
            filter(minute(datetime) %% synch_intvl == 0) %>%
            tidyr::complete(datetime = seq(min(datetime), max(datetime),
                                           by = paste(synch_intvl, 'min')))

    } else {

        gap_stats[i, 'n_fill'] <- 0
        gap_stats[i, 'n_nodonor'] <- gap_stats[i, 'n_tot']
        gap_stats[i, 'n_doublegap'] <- gap_stats[i, 'n_gap']
        gap_stats[i, 'nday_fill'] <- 0
        gap_stats[i, 'nday_nodonor'] <- gap_stats[i, 'nday_tot']
        gap_stats[i, 'nday_doublegap'] <- gap_stats[i, 'nday_gap']
        next
    }

    dall <- left_join(d, dfill, by = 'datetime')

    pgdl <- any(grepl('G|S|PG|PS', unlist(ranks[ranks$site == s, 2:4])))
    if(pgdl){

        dall <- prepare_q_lstm(s) %>%
            select(datetime) %>% #, Q_predicted_lstm = discharge_Ls) %>%
            tidyr::complete(datetime = seq(as_datetime(as.Date(datetime[1])),
                                           as_datetime(as.Date(datetime[nrow(.)]) + 1),
                                           by = '15 min')) %>%
            mutate(Q_predicted_lstm = 1) %>%
            right_join(dall, by = 'datetime')
    } else {
        dall$Q_predicted_lstm <- NA_real_
    }

    gap_stats[i, 'n_fill'] <- sum(is.na(dall$discharge) & ! is.na(dall$Q_predicted))
    gap_stats[i, 'n_nodonor'] <- sum(is.na(dall$Q_predicted))
    gap_stats[i, 'n_noreconst'] <- sum(is.na(dall$Q_predicted) & is.na(dall$Q_predicted_lstm) & ! is.na(dall$discharge))
    gap_stats[i, 'n_doublegap'] <- sum(is.na(dall$discharge) & is.na(dall$Q_predicted))
    gap_stats[i, 'nday_fill'] <- round(gap_stats[i, 'n_fill'] / nsamp_per_day, 2)
    gap_stats[i, 'nday_nodonor'] <- round(gap_stats[i, 'n_nodonor'] / nsamp_per_day, 2)
    gap_stats[i, 'nday_noreconst'] <- round(gap_stats[i, 'n_noreconst'] / nsamp_per_day, 2)
    gap_stats[i, 'nday_doublegap'] <- round(gap_stats[i, 'n_doublegap'] / nsamp_per_day, 2)
    gap_stats[i, 'n_fill_lstm'] <- sum(is.na(dall$discharge) & is.na(dall$Q_predicted) & ! is.na(dall$Q_predicted_lstm))
    gap_stats[i, 'nday_fill_lstm'] <- round(gap_stats[i, 'n_fill_lstm'] / nsamp_per_day, 2)
    gap_stats[i, 'n_triplegap'] <- sum(is.na(dall$discharge) & is.na(dall$Q_predicted) & is.na(dall$Q_predicted_lstm))
    gap_stats[i, 'nday_triplegap'] <- round(gap_stats[i, 'n_triplegap'] / nsamp_per_day, 2)
}

gap_stats[, 'pct_gap'] <- round(gap_stats[, 'n_gap'] / gap_stats[, 'n_tot'], 2)
gap_stats[, 'pct_fill'] <- round(gap_stats[, 'n_fill'] / gap_stats[, 'n_gap'], 2)
gap_stats[, 'pct_nodonor'] <- round(gap_stats[, 'n_nodonor'] / gap_stats[, 'n_tot'], 2)
gap_stats[, 'pct_noreconst'] <- round(gap_stats[, 'n_noreconst'] / gap_stats[, 'n_tot'], 2)
gap_stats[, 'pct_doublegap'] <- round(gap_stats[, 'n_doublegap'] / gap_stats[, 'n_tot'], 2)
gap_stats[, 'pct_fill_lstm'] <- round(gap_stats[, 'n_fill_lstm'] / gap_stats[, 'n_gap'], 2)
gap_stats[, 'pct_triplegap'] <- round(gap_stats[, 'n_triplegap'] / gap_stats[, 'n_tot'], 2)


(gap_sums <- colSums(gap_stats[, c('nday_tot', 'nday_gap', 'nday_fill',
                                   'nday_nodonor', 'nday_doublegap', 'nday_fill_lstm',
                                   'nday_triplegap')],
                     na.rm = TRUE))

(gapdays_informed <- gap_sums[3] + gap_sums[6])
(gapmonths_informed <- unname(gapdays_informed / 30))
# days without reconstitutions, but with NEON estimates
sum(gap_stats[, 'nday_noreconst'], na.rm = TRUE)
# % gapdays filled with high-res estimates
unname(gap_sums[3] / gap_sums[2])
# % gapdays filled by LSTM or high-res
unname(gapdays_informed / gap_sums[2])

write.csv(gap_stats, 'out/gap_stats.csv', row.names = TRUE)

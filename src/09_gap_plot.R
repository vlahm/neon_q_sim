# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-03-08

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
if(! dir.exists('out/lm_out')) source('src/02_linear_regression.R', local = new.env())
if(! dir.exists('in/lstm_data')) source('src/03_organize_camels_macrosheds_nhm.R', local = new.env())
if(! dir.exists('out/lstm_runs')) stop("you need to run src/04_run_lstms.R. It will take many days unless run on a cluster. Or use our bundled results.")


## 1. data prep ####

neon_sites <- read_csv('in/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    pull(SiteID)

#neon-official continuous Q data
plotd <- list()
for(s in neon_sites){
    plotd[[s]] <- read_csv(paste0('in/neon_continuous_Q/', s, '.csv'))
}

#linear regression predictions
plotfill <- list()
for(s in neon_sites){

    # if(s == 'TECR'){
    #     plotfill[[s]] <- read_csv(paste0('out///', s, '.csv'))
    # }

    plotfill[[s]] <- try(read_csv(paste0('out/lm_out/predictions/', s, '.csv')),
                         silent = TRUE)
}

#lstm predictions
HERE
gen_res <- retrieve_test_results(generalist_runids)
plotfill$TECR
plotfill$BIGC

## 2. plot ####

neon_sites <- sort(neon_sites)

png(width = 8, height = 12, units = 'in', type = 'cairo', res = 300,
    filename = 'figs/fig4v2.png')
par(mfrow = c(27, 1), mar = c(1, 3, 0, 0), oma = c(2, 2, 2, 10))
for(i in seq_along(neon_sites)){

    s <- neon_sites[i]

    d <- plotd[[s]] %>%
        filter(minute(datetime) %% 5 == 0) %>%
        mutate(discharge = neglog(discharge)) %>%
        tidyr::complete(datetime = seq(min(datetime), max(datetime), by = '5 min'))

    if(! inherits(plotfill[[s]], 'try-error')){
        dfill <- plotfill[[s]] %>%
            mutate(across(any_of('date'), as_datetime)) %>%
            rename_with(function(x) sub('^date$', 'datetime', x)) %>%
            mutate(Q_predicted = neglog(Q_predicted)) %>%
            filter(datetime > min(d$datetime),
                   datetime < max(d$datetime))
    } else {
        dfill <- NULL
    }
        # tidyr::complete(datetime = seq(min(datetime), max(datetime), by = '5 min'))

    plot(d$datetime, d$discharge, type = 'n',
         ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n')

    if(! is.null(dfill)){
        dboth <- left_join(dfill, d, by = 'datetime') %>%
            mutate(still_missing = is.na(discharge) & is.na(Q_predicted))
        dboth$still_missing[dboth$still_missing] <- NA

        polygon_with_gaps2(dfill, 'Q_predicted', -100, 999999999, '#FFBDF9')

        polygon_with_gaps2(d, 'discharge', -100, 999999999, '#bdd9ff')
    } else {
        polygon_with_gaps2(d, 'discharge', -100, 999999999, '#FFBDF9', invert = TRUE)
        polygon_with_gaps2(d, 'discharge', -100, 999999999, '#c09eff')
    }

    if(! is.null(dfill)){
        polygon_with_gaps2(dboth, 'still_missing', -100, 999999999, '#c09eff')
    }

    if(s == 'TOMB'){
        dtomb <- filter(d, ! is.na(discharge))
        lines(dtomb$datetime, dtomb$discharge, col = 'gray60')
    } else {
        lines(d$datetime, d$discharge, col = 'gray60')
    }

    if(! is.null(dfill)){

        rl = rle(! is.na(d$discharge))
        vv = ! rl$values
        chunkfac = rep(cumsum(vv), rl$lengths)
        chunkfac[chunkfac == 0] = 1
        chunks = split(d, chunkfac)
        NAchunks = lapply(chunks, function(x) x[is.na(x$discharge), ])
        NAchunks = Filter(function(x) difftime(x$datetime[nrow(x)], x$datetime[1], units = 'hours') >= 6,
                          NAchunks)

        if(length(NAchunks)){

            for(j in 1:length(NAchunks)){

                missing_dts <- NAchunks[[j]]$datetime[is.na(NAchunks[[j]]$discharge)]
                points(dfill$datetime[dfill$datetime %in% missing_dts],
                       dfill$Q_predicted[dfill$datetime %in% missing_dts],
                       col = 'black', pch = '.')
            }
        }
    }

    if(i == 1){
        legend(x = as.POSIXct('2021-09-30'), y = -90, bty = 'n', border = FALSE, cex = 1.5,
               legend = c('NEON\nmissing', 'Reconstruction\nmissing',
                          'Both\nmissing'),
               fill = c('#bdd9ff', '#FFBDF9', '#c09eff'),
               xpd = NA, x.intersp = 1.1, y.intersp = 2)
        legend(x = as.POSIXct('2021-09-19'), y = -125, bty = 'n', border = FALSE, cex = 1.5,
               legend = c('NEON\ndischarge', 'Reconstruction\ngapfill'),
               col = c('gray60', 'black'), lty = 1, lwd = 2, xpd = NA,
               seg.len = 1.2, y.intersp = 2)
    }

    mtext(s, side = 2, las = 1, line = -1.2)
}

mtext(expression('Log Discharge'), side = 2, outer = TRUE, cex = 1.2,
      line = -0.5)

dev.off()


## 3. gap stats ####

gap_stats <- matrix(
    NA_real_,
    nrow = length(neon_sites),
    ncol = 17,
    dimnames = list(
        neon_sites,
        c('t_tot', 't_gap', 't_fill', 't_nodonor', 't_doublegap',
          'pct_gap', 'pct_fill', 'pct_nodonor', 'pct_doublegap',
          'neon_intvl_orig', 'est_intvl_orig', 'cast_intvl',
          'nday_tot', 'nday_gap', 'nday_fill', 'nday_nodonor', 'nday_doublegap')
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

    gap_stats[i, 't_tot'] <- nrow(d)
    gap_stats[i, 't_gap'] <- sum(is.na(d$discharge))
    gap_stats[i, 'nday_tot'] <- round(gap_stats[i, 't_tot'] / nsamp_per_day, 2)
    gap_stats[i, 'nday_gap'] <- round(gap_stats[i, 't_gap'] / nsamp_per_day, 2)

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
        gap_stats[i, 't_fill'] <- 0
        gap_stats[i, 't_nodonor'] <- gap_stats[i, 't_tot']
        gap_stats[i, 't_doublegap'] <- gap_stats[i, 't_gap']
        gap_stats[i, 'nday_fill'] <- 0
        gap_stats[i, 'nday_nodonor'] <- gap_stats[i, 'nday_tot']
        gap_stats[i, 'nday_doublegap'] <- gap_stats[i, 'nday_gap']
        next
    }

    dboth <- left_join(d, dfill, by = 'datetime')

    gap_stats[i, 't_fill'] <- sum(is.na(dboth$discharge) & ! is.na(dboth$Q_predicted))
    gap_stats[i, 't_nodonor'] <- sum(is.na(dboth$Q_predicted))
    gap_stats[i, 't_doublegap'] <- sum(is.na(dboth$discharge) & is.na(dboth$Q_predicted))
    gap_stats[i, 'nday_fill'] <- round(gap_stats[i, 't_fill'] / nsamp_per_day, 2)
    gap_stats[i, 'nday_nodonor'] <- round(gap_stats[i, 't_nodonor'] / nsamp_per_day, 2)
    gap_stats[i, 'nday_doublegap'] <- round(gap_stats[i, 't_doublegap'] / nsamp_per_day, 2)
}

gap_stats[, 'pct_gap'] <- round(gap_stats[, 't_gap'] / gap_stats[, 't_tot'], 2)
gap_stats[, 'pct_fill'] <- round(gap_stats[, 't_fill'] / gap_stats[, 't_gap'], 2)
gap_stats[, 'pct_nodonor'] <- round(gap_stats[, 't_nodonor'] / gap_stats[, 't_tot'], 2)
gap_stats[, 'pct_doublegap'] <- round(gap_stats[, 't_doublegap'] / gap_stats[, 't_tot'], 2)

(gap_sums <- colSums(gap_stats[, c('nday_tot', 'nday_gap', 'nday_fill', 'nday_nodonor', 'nday_doublegap')]))
unname(gap_sums[3] / gap_sums[2])

write.csv(gap_stats, 'out/gap_stats.csv')

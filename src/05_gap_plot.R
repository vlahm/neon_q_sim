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

# set working directory to same location as in 01_neon_q_sim.R
setwd('~/git/macrosheds/papers/q_sim')

source('src/00_helpers.R')

## 1. load data ####

# NEON discharge data (field measurements and continuous)

if(! length(list.files('in/neon_continuous_Q/'))){
    get_neon_inst_discharge(neon_sites)
}

if(! file.exists('in/neon_site_info.csv')){
    download.file('https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/data/contents/neon_site_info.csv',
                  destfile = 'in/neon_site_info.csv')
}

neon_sites <- read_csv('in/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    pull(SiteID)

## 2. data prep ####

plotd <- list()
for(s in neon_sites){
    plotd[[s]] <- read_csv(paste0('in/neon_continuous_Q/', s, '.csv'))
}

plotfill <- list()
for(s in neon_sites){

    # if(s == 'TECR'){
    #     plotfill[[s]] <- read_csv(paste0('out///', s, '.csv'))
    # }

    plotfill[[s]] <- try(read_csv(paste0('out/lm_out/predictions/', s, '.csv')),
                         silent = TRUE)
}

## 3. plot ####

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

mtext(expression('Log Discharge (Ls'^-1*')'), side = 2, outer = TRUE, cex = 1.2,
      line = -0.5)

dev.off()


## 4. gap stats ####


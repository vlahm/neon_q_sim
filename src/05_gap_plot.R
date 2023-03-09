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

# 2. data prep ####

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

# 3. plot ####

png(width = 8, height = 5, units = 'in', type = 'cairo', res = 300,
    filename = 'figs/fig4.png')
par(mfrow = c(length(neon_sites) / 2, 2), mar = c(2, 3, 1, 1), oma = c(2, 2, 2, 0))
for(s in neon_sites){

    d <- plotd[[s]] %>%
        filter(minute(datetime) %% 5 == 0) %>%
        mutate(discharge = neglog(discharge)) %>%
        tidyr::complete(datetime = seq(min(datetime), max(datetime), by = '5 min'))

    dfill <- plotfill[[s]] %>%
        mutate(across(any_of('date'), as_datetime)) %>%
        rename_with(function(x) sub('^date$', 'datetime', x)) %>%
        mutate(Q_predicted = neglog(Q_predicted)) %>%
        filter(datetime > min(d$datetime),
               datetime < max(d$datetime))
        # tidyr::complete(datetime = seq(min(datetime), max(datetime), by = '5 min'))

    missing_dts <- d$datetime[is.na(d$discharge)] ONLY PLOT POINTS IN BLUE POLYGONS

    plot(d$datetime, d$discharge, type = 'n',
         ylab = '', xlab = '', xaxt = 'n', yaxt = 'n', bty = 'n')

    if(! is.null(dfill)){
        dboth <- left_join(dfill, d, by = 'datetime') %>%
            mutate(still_missing = is.na(discharge) & is.na(Q_predicted))
        dboth$still_missing[dboth$still_missing] <- NA

        polygon_with_gaps2(dfill, 'Q_predicted', -100, 999999999, '#FFBDF9')
    }

    polygon_with_gaps2(d, 'discharge', -100, 999999999, '#bdd9ff')

    if(! is.null(dfill)){
        polygon_with_gaps2(dboth, 'still_missing', -100, 999999999, '#c09eff')
    }

    lines(d$datetime, d$discharge, col = 'gray60')

    if(! is.null(dfill)){
        points(dfill$datetime[dfill$datetime %in% missing_dts],
              dfill$Q_predicted[dfill$datetime %in% missing_dts],
              col = 'black', pch = '.')
    }

    mtext(s, side = 3, line = 0.2, adj = 0)
}
# mtext(expression('Discharge (Ls'^-1*')'), side = 2, outer = TRUE)
dev.off()

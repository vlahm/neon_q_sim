# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-02-28

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

if(! length(list.files('in/neon_field_Q/'))){
    get_neon_field_discharge(neon_sites)
}

# lm predictions


png(width = 8, height = 5, units = 'in', type = 'cairo', res = 300,
    filename = 'figs/fig3.png')
par(mfrow = c(2, 2), mar = c(2, 3, 1, 1), oma = c(2, 2, 0, 0))
ts_plot('REDB', 2019) #0.932
legend('topright', legend = c('Predicted', 'NEON modeled  '),
       col = c('red', 'gray50'), lty = 1, bty = 'n')
legend('topright', legend = c('', '','NEON measured'),
       col = c('transparent', 'transparent', 'black'), pch = 1, bty = 'n')
ts_plot('FLNT', 2018) #0.970
# ts_plot('BLUE', 2020, T)
ts_plot('LEWI', 2019) #0.739
ts_plot('WALK', 2017) #0.886
mtext(expression('Discharge (Ls'^-1*')'), side = 2, outer = TRUE)
dev.off()

# dg = dygraphs::dygraph(xts(x = select(plotd, Q_predicted, Q_neon_field) %>% tail(5e5),
#                            order.by = tail(plotd$datetime, 5e5))) %>%
#     dygraphs::dyRangeSelector()

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

#pre-bundled in/out data available at: [**]
if(! exists('ts_plot')) source('src/00_helpers.R')
if(! exists('ms_areas')) source('src/01_data_retrieval.R')
if(! dir.exists('out/lm_out')) source('src/02_linear_regression.R', local = new.env())
if(! dir.exists('in/lstm_data')) source('src/03_organize_camels_macrosheds_nhm.R', local = new.env())
if(! dir.exists('out/lstm_runs')) source('src/04_run_lstms.R', local = new.env())

# 1. plot ####

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

# 2. site table ####

sitelist <- read_csv('in/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    select(SiteID)

read_csv('in/neon_site_info2.csv') %>%
    right_join(sitelist, by = c('field_site_id' = 'SiteID')) %>%
    arrange(desc(field_watershed_size_km2)) %>%
    mutate(field_site_name = sub(' NEON$', '', field_site_name),
           field_watershed_size_km2 = round(field_watershed_size_km2, 1)) %>%
    select(`Site code` = field_site_id, `Full name` = field_site_name,
           `State (USA)` = field_site_state, `Watershed area (km2)` = field_watershed_size_km2,
           `Mean elevation (m)` = field_mean_elevation_m) %>%
    write_csv('out/site_table.csv')

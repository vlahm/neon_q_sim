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

# 0. setup ####

if(! length(list.files('in/NEON/neon_continuous_Q_withflags/'))){
    get_neon_inst_discharge(neon_sites, clean_only = FALSE)
}

# 1. plot ####

png(width = 8, height = 5, units = 'in', type = 'cairo', res = 300,
    filename = 'figs/ts_plot.png')
par(mfrow = c(2, 2), mar = c(2, 3, 1, 1), oma = c(2, 2, 0, 0))
ts_plot('REDB', 2019)
legend('topright', legend = c('Predicted', 'NEON modeled  '),
       col = c('red', 'gray50'), lty = 1, bty = 'n')
legend('topright', legend = c('', '','NEON measured', 'NEON QC flag'),
       col = c('transparent', 'transparent', 'black', 'black'),
       pch = c(1, 1, 1, 39), bty = 'n')
ts_plot('FLNT', 2018)
ts_plot('LEWI', 2019, border = adjustcolor('red', alpha.f = 0.2))
ts_plot('WALK', 2017, border = adjustcolor('red', alpha.f = 0.2))
mtext(expression('Discharge (Ls'^-1*')'), side = 2, outer = TRUE)
dev.off()

# 2. site table ####

sitelist <- read_csv('in/NEON/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    select(SiteID)

read_csv('in/NEON/neon_site_info2.csv') %>%
    right_join(sitelist, by = c('field_site_id' = 'SiteID')) %>%
    arrange(desc(field_watershed_size_km2)) %>%
    mutate(field_site_name = sub(' NEON$', '', field_site_name),
           field_watershed_size_km2 = round(field_watershed_size_km2, 1)) %>%
    select(`Site code` = field_site_id, `Full name` = field_site_name,
           `State (USA)` = field_site_state, `Watershed area (km2)` = field_watershed_size_km2,
           `Mean elevation (m)` = field_mean_elevation_m) %>%
    write_csv('out/site_table.csv')

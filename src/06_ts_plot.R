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

source('src/00_helpers.R')

## 1. load data ####

# NEON discharge data (field measurements and continuous)

if(! length(list.files('in/neon_continuous_Q/'))){
    get_neon_inst_discharge(neon_sites)
}

if(! length(list.files('in/neon_field_Q/'))){
    get_neon_field_discharge(neon_sites)
}

if(! file.exists('in/neon_site_info.csv')){
    download.file('https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/data/contents/neon_site_info.csv',
                  destfile = 'in/neon_site_info.csv')
}

if(! file.exists('in/neon_site_info2.csv')){
    #filename changes with every update, so might have to modify URL below
    download.file('https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv',
                  destfile = 'in/neon_site_info2.csv')
}

# 2. plot ####

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

# 3. site table ####

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

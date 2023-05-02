# Mike Vlah
# vlahm13@gmail.com
# last data retrieval dates given in comments below:
# last edit: 2023-03-28

library(dataRetrieval)
library(glue)
library(neonUtilities)
library(tidyverse)
library(macrosheds)
library(yaml)

# NOTE: working directory is set by opening q_sim.Rproj. Alternatively, you may set it here.
# setwd('path/to/q_sim')

# ALSO NOTE: you don't need to download and build all these components piecemeal.
# We've bundled it on Figshare ([**]).
# We used a previous version of the CAMELS dataset for
# training LSTMs, so some paths will certainly break if you use the current version.

options(readr.show_progress = FALSE,
        readr.show_col_types = FALSE,
        timeout = 3000)

source('src/00_helpers.R')

## 1. retrieve NEON data ####

dir.create('in/NEON', showWarnings = FALSE)

# NEON site metadata as reproduced by Rhea et al. 2023
# (primary source here: https://www.neonscience.org/field-sites/explore-field-sites)

#last retrieval: #2023-04-14
if(! file.exists('in/NEON/neon_site_info.csv')){
    download.file('https://www.hydroshare.org/resource/1a388391632f4277992889e2de152163/data/contents/neon_site_info.csv',
                  destfile = 'in/NEON/neon_site_info.csv')
}

#last retrieval: #2023-01-31
if(! file.exists('in/NEON/neon_site_info2.csv')){
    #filename changes with every update, so might have to modify URL below
    download.file('https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv',
                  destfile = 'in/NEON/neon_site_info2.csv')
}

#last retrieval: #2023-01-31
neon_areas <- read_csv('in/NEON/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    mutate(ws_area_ha = WSAreaKm2 * 100) %>%
    select(site_code = SiteID, ws_area_ha)

neon_sites <- neon_areas$site_code

# NEON discharge data (field measurements and continuous)

#last retrieval: #2023-05-01
if(! length(list.files('in/NEON/neon_continuous_Q/'))){
    get_neon_inst_discharge(neon_sites)
}

#last retrieval: #2023-01-31
if(! length(list.files('in/NEON/neon_field_Q/'))){
    get_neon_field_discharge(neon_sites)
}

# NEON discharge evaluation results from Rhea et al. 2023

#last retrieval: #2023-04-14
if(! file.exists('in/NEON/neon_q_eval.csv')){
    download.file('https://www.hydroshare.org/resource/1a388391632f4277992889e2de152163/data/contents/neon_q_eval.csv',
                  destfile = 'in/NEON/neon_q_eval.csv')
}

#last retrieval: #2023-01-31
q_eval <- read_csv('in/NEON/neon_q_eval.csv') %>%
    filter(site %in% neon_sites)


## 2. retrieve additional input data for linear regression ####

dir.create('in/macrosheds', showWarnings = FALSE)

# MacroSheds Q data from Niwot domain: used for donor gauges in lieu of USGS data

#last retrieval: #2023-01-31
if(! dir.exists('in/macrosheds/niwot')){
    macrosheds::ms_download_core_data(macrosheds_root = './in/macrosheds',
                                      domains = 'niwot')
}

ms_q <- macrosheds::ms_load_product(macrosheds_root = './in/macrosheds',
                                    prodname = 'discharge',
                                    domains = 'niwot') %>%
    filter(ms_status == 0)

# H. J. Andrews Experimental Forest Q data: used for donor gauges in lieu of USGS data

#last retrieval: #2022-04-14
if(! file.exists('in/hjandrews_q.txt')){
    download.file('https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-and.4341.33&entityid=86490799297ff361b0741b807804c43a',
                  destfile = 'in/hjandrews_q.txt')
}

# relevant MacroSheds site metadata

ms_areas <- ms_load_sites() %>%
    filter(site_type == 'stream_gauge',
           ! is.na(ws_area_ha)) %>%
    select(site_code, ws_area_ha) %>%
    distinct(site_code, .keep_all = TRUE) %>%
    bind_rows(neon_areas)

# donor gauge IDs

donor_gauges <- read_yaml('cfg/donor_gauges.yml')

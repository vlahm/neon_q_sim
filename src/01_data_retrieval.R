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
library(ncdf4)
library(feather)

# NOTE: working directory is set by opening q_sim.Rproj. Alternatively, you may set it here.
# setwd('path/to/q_sim')

# ALSO NOTE: you don't need to download and build all these components piecemeal.
# We've bundled it on Figshare ([]).
# We used a previous version of the CAMELS dataset for
# training LSTMs, so some paths will certainly break if you use the current version.

options(readr.show_progress = FALSE,
        readr.show_col_types = FALSE,
        timeout = 3000)

source('src/00_helpers.R')

## 1. retrieve NEON data ####

# NEON site metadata as reproduced by Rhea et al. 2023
# (primary source here: https://www.neonscience.org/field-sites/explore-field-sites)

#last retrieval: #2023-01-31
if(! file.exists('in/neon_site_info.csv')){
    download.file('https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/data/contents/neon_site_info.csv',
                  destfile = 'in/neon_site_info.csv')
}

#last retrieval: #2023-01-31
if(! file.exists('in/neon_site_info2.csv')){
    #filename changes with every update, so might have to modify URL below
    download.file('https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv',
                  destfile = 'in/neon_site_info2.csv')
}

#last retrieval: #2023-01-31
neon_areas <- read_csv('in/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    mutate(ws_area_ha = WSAreaKm2 * 100) %>%
    select(site_code = SiteID, ws_area_ha)

neon_sites <- neon_areas$site_code

# NEON discharge data (field measurements and continuous)

#last retrieval: #2023-03-09
if(! length(list.files('in/neon_continuous_Q/'))){
    get_neon_inst_discharge(neon_sites)
}

#last retrieval: #2023-01-31
if(! length(list.files('in/neon_field_Q/'))){
    get_neon_field_discharge(neon_sites)
}

# NEON discharge evaluation results from Rhea et al. 2023

#last retrieval: #2023-01-31
if(! file.exists('in/neon_q_eval.csv')){
    download.file('https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/data/contents/neon_q_eval.csv',
                  destfile = 'in/neon_q_eval.csv')
}

#last retrieval: #2023-01-31
q_eval <- read_csv('in/neon_q_eval.csv') %>%
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

ms_areas <- macrosheds::ms_load_sites() %>%
    filter(site_type == 'stream_gauge',
           ! is.na(ws_area_ha)) %>%
    select(site_code, ws_area_ha)

# donor gauge IDs

donor_gauges <- read_yaml('in/donor_gauges.yml')


## 3. retrieve additional input data for LSTM simulation ####

#the CAMELS dataset is available here: https://gdex.ucar.edu/dataset/camels.html,
#but note that the current version is organized differently from the one we used,
#which was downloaded on 2021-12-11. We recommend just downloading our bundled input data from
#[]

#National Hydrologic Model (v1) data were provided to us for each CAMELS
#and MacroSheds site by the USGS. A minimal copy of this dataset (omitting
#files superfluous to this analysis) is included with our bundled input data
#on []

#retrieve MacroSheds core dataset and CAMELS-compliant MacroSheds attributes

ms_sites <- ms_load_sites()

if(! dir.exists('in/macrosheds/arctic')){

    ms_download_core_data(macrosheds_root = './in/macrosheds',
                          domains = 'all', version = '1.0')

    ms_download_ws_attr(macrosheds_root = './in/macrosheds',
                        dataset = 'all', version = '1.0')
}

ms_q <- ms_load_product(macrosheds_root = './in/macrosheds',
                        prodname = 'discharge', warn = FALSE) %>%
    filter(ms_status == 0)

ms_attr_static <- ms_load_product(macrosheds_root = './in/macrosheds',
                                  prodname = 'ws_attr_CAMELS_summaries',
                                  warn = FALSE)

ms_attr2_dyn <- ms_load_product(macrosheds_root = './in/macrosheds',
                                prodname = 'ws_attr_CAMELS_Daymet_forcings',
                                warn = FALSE)

## 4. reshape CAMELS, MacroSheds, NHM datasets for NeuralHydrology

ms_attr_static <- ms_attr_static %>%
    bind_rows(read_csv('in/macrosheds/neon_attrs.csv')) %>%
    rename(q5 = Q5, q95 = Q95, baseflow_index = baseflow_index_landson)

nhm_sites <- read_csv('in/NHMv1/sites_with_segids.csv')

ms_attr_static <- ms_attr_static %>%
    filter(domain == 'neon') %>%
    mutate(site_code = paste0(site_code, '_MANUALQ')) %>%
    bind_rows(.)

ms_attr_static <- ms_attr_static %>%
    filter(domain == 'neon') %>%
    mutate(site_code = paste0(site_code, '_MANUALQ'))

ms_attr_static %>%
    filter(site_code %in% nhm_sites$site_code) %>%
    mutate(site_code = paste0('NHM_', site_code))

    select(-domain, -network) %>%
write_csv(ms_attr_static

#

camels_clim <- read_delim('in/CAMELS/camels_attributes_v2.0/camels_clim.txt',
                          col_types = 'cnnnnnnncnnc', delim = ';')
camels_soil <- read_delim('in/CAMELS/camels_attributes_v2.0/camels_soil.txt',
                          col_types = 'cnnnnnnnnnnn', delim = ';')
camels_geol <- read_delim('in/CAMELS/camels_attributes_v2.0/camels_geol.txt',
                          col_types = 'ccncnnnn', delim = ';') %>%
    rename(geol_porosity = geol_porostiy)
camels_topo <- read_delim('in/CAMELS/camels_attributes_v2.0/camels_topo.txt',
                          col_types = 'cnnnnnn', delim = ';')
camels_vege <- read_delim('in/CAMELS/camels_attributes_v2.0/camels_vege.txt',
                          col_types = 'cnnnnnncnn', delim = ';')
camels_hydro <- read_delim('in/CAMELS/camels_attributes_v2.0/camels_hydro.txt',
                          col_types = 'cnnnnnnnnnnnn', delim = ';')

camels_attrs <- camels_geol %>%
    full_join(camels_soil, by = 'gauge_id') %>%
    full_join(camels_topo, by = 'gauge_id') %>%
    full_join(camels_vege, by = 'gauge_id') %>%
    full_join(camels_clim, by = 'gauge_id') %>%
    full_join(camels_hydro, by = 'gauge_id') %>%
    rename(site_code = gauge_id, area = area_gages2)

# we were unable to exactly replicate some CAMELS methods when generating attributes for
# macrosheds sites. for consistency, we replace those original CAMELS attrs
# with our own recomputed ones (pet_mean, aridity, frac_snow, soil fractions)
camels_mod_attrs <- full_join(
    read_feather('in/CAMELS/recomputed_attributes/clim.feather'),
    read_feather('in/CAMELS/recomputed_attributes/soil.feather'),
    by = 'site_code'
)

replace_cols <- grep('site_code', colnames(camels_mod_attrs), value = TRUE, invert = TRUE)

camels_attrs <- camels_attrs %>%
    select(-all_of(replace_cols)) %>%
    full_join(camels_mod_attrs, by = 'site_code') %>%
    select(-all_of(setdiff(colnames(camels_attrs), colnames(ms_attr_static))))


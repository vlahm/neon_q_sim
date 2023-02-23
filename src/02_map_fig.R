# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-02-22

library(dataRetrieval)
library(tidyverse)
library(macrosheds)
library(tmap)
library(yaml)
library(sf)

# set working directory to same location as in 01_neon_q_sim.R
setwd('~/git/macrosheds/papers/q_sim')

## 1. load data ####

if(! file.exists('in/neon_site_info.csv')){
    download.file('https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/data/contents/neon_site_info.csv',
                  destfile = 'in/neon_site_info.csv')
}

neon_sites <- read_csv('in/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    st_as_sf(coords = c('Longitude', 'Latitude'))

ms_site_codes_in_use <- c('GSWS01', 'GSWS06', 'GSWS07', 'GSWS08', 'GSLOOK',
                          'ALBION', 'SADDLE')

ms_sites <- macrosheds::ms_load_sites() %>%
    filter(site_type == 'stream_gauge',
           site_code %in% ms_site_codes_in_use) %>%
    st_as_sf(coords = c('longitude', 'latitude'))

donor_gauges <- yaml::read_yaml('in/donor_gauges.yml')

usgs_sites <- lapply(donor_gauges, function(x) Filter(function(y) ! grepl('[A-Z]', y), x))
usgs_sites <- usgs_sites[sapply(usgs_sites, length) > 0]
usgs_sites <- lapply(usgs_sites, function(x) sapply(x, function(y) tibble(geometry = gaugeid_to_location(y))))
usgs_sites <- lapply(usgs_sites, function(x) Reduce(bind_rows, x))
usgs_sites_flattened <- Reduce(bind_rows, usgs_sites)

## 2. main map ####

basemap_main <- tmaptools::read_osm(usgs_sites_flattened)
main_map <- tm_shape(neon_sites) + tm_symbols(shape = 1, col = 'red', size = 0.8, border.lwd = 2) +
    tm_shape(ms_sites) + tm_symbols(shape = 1, col = 'blue', size = 0.8, border.lwd = 2) +
    tm_shape(usgs_sites_flattened) + tm_symbols(shape = 1, col = 'green', size = 0.8, border.lwd = 2) +
    tm_shape(basemap_main) + tm_basemap()
    # tm_shape(basemap_main) + tm_raster()
main_map

# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-02-22

library(dataRetrieval)
library(tidyverse)
library(macrosheds)
library(tmap)
library(yaml)
library(sf)
library(terra)
library(osmdata)

#pre-bundled in/out data available at: [**]
if(! exists('ts_plot')) source('src/00_helpers.R')
if(! exists('ms_areas')) source('src/01_data_retrieval.R')

## 1. load data ####

neon_sites <- read_csv('in/NEON/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    st_as_sf(coords = c('Longitude', 'Latitude'))

ms_site_codes_in_use <- c('GSWS01', 'GSWS06', 'GSWS07', 'GSWS08', 'GSLOOK',
                          'ALBION', 'SADDLE')

ms_sites <- macrosheds::ms_load_sites() %>%
    filter(site_type == 'stream_gauge',
           site_code %in% ms_site_codes_in_use) %>%
    st_as_sf(coords = c('longitude', 'latitude'), crs = 4326)

usgs_sites <- lapply(donor_gauges, function(x) Filter(function(y) ! grepl('[A-Z]', y), x))
usgs_sites <- usgs_sites[sapply(usgs_sites, length) > 0]
usgs_sites <- lapply(usgs_sites, function(x) sapply(x, function(y) tibble(geometry = gaugeid_to_location(y))))
usgs_sites <- lapply(usgs_sites, function(x) Reduce(bind_rows, x))
usgs_sites_flattened <- Reduce(bind_rows, usgs_sites)

## 2. Figure 1: map ####

tmap_mode("plot")

neoncol = 'goldenrod3'

bbox <- st_bbox(st_buffer(usgs_sites_flattened, dist = 2e5))
basemap_main <- tmaptools::read_osm(bbox, type = 'osm-public-transport')
main_map <-
    tm_shape(basemap_main) + tm_rgb() +
    tm_shape(neon_sites) + tm_symbols(shape = 19, col = neoncol, size = 0.5, border.lwd = 2, alpha=0.5) +
    tm_compass(type="arrow", position=c("left", "bottom"), show.labels = 0, size = 1.2) +
    tm_scale_bar(position=c(0.1, 0.01), breaks = c(0, 1000, 2000), text.size = 0.6) +
    tm_add_legend(type='symbol', labels = c('NEON', 'USGS', 'MacroSheds'),
                  col = c(neoncol, 'darkorange4', 'purple4'), size = 0.5, shape=19) +
    tm_legend(show=TRUE, position=c('right', 'top'), bg.color='gray97', frame = TRUE, height = -0.115)

tmap_save(main_map, filename='figs/map_components/main_map.png', bg="white",
          dpi = 600, height = 5, width = 5, units = 'in')
# tmap_save(main_map, filename='figs/map_components/main_map.pdf', bg="white",
#           dpi = 300, device = cairo_pdf, height = 5, width = 5, units = 'in')

neon_site <- filter(neon_sites, SiteID == 'REDB')
inset_sites <- bind_rows(usgs_sites$REDB, neon_site)
bbox <- st_bbox(st_buffer(inset_sites, dist = 5e3))
basemap_redb <- tmaptools::read_osm(bbox, type = 'stamen-terrain')
streams_layer <- get_osm_streams(st_as_sf(basemap_redb))
inset_redb <-
    tm_shape(basemap_redb) + tm_rgb() +
    tm_shape(streams_layer) + tm_lines(col = 'steelblue3', lwd = 1) +
    tm_shape(usgs_sites$REDB) + tm_symbols(shape = 19, col = 'darkorange4', size = 0.7, border.lwd = 2, alpha=0.6) +
    tm_shape(neon_site) + tm_symbols(shape = 19, col = neoncol, size = 0.7, border.lwd = 2, alpha=0.9) +
    # tm_compass(type="arrow", position=c("left", "bottom"), show.labels = 0, size = 1) +
    tm_scale_bar(position=c(0.05, 0.01), width = 0.25, text.size = 0.7)

tmap_save(inset_redb, filename='figs/map_components/redb.png', bg="transparent",
          dpi = 300, height = 2.5, width = 2.5, units = 'in')

neon_site <- filter(neon_sites, SiteID == 'GUIL')
inset_sites <- bind_rows(usgs_sites$GUIL, neon_site)
bbox <- st_bbox(st_buffer(inset_sites, dist = 1e4))
basemap_guil <- tmaptools::read_osm(bbox, type = 'stamen-terrain')
streams_layer <- get_osm_streams(st_as_sf(basemap_guil))
inset_guil <-
    tm_shape(basemap_guil) + tm_rgb() +
    tm_shape(streams_layer) + tm_lines(col = 'steelblue3', lwd = 1) +
    tm_shape(usgs_sites$GUIL) + tm_symbols(shape = 19, col = 'darkorange4', size = 0.7, border.lwd = 2, alpha=0.6) +
    tm_shape(neon_site) + tm_symbols(shape = 19, col = neoncol, size = 0.7, border.lwd = 2, alpha=0.9) +
    # tm_compass(type="arrow", position=c("right", "bottom"), show.labels = 0, size = 1) +
    tm_scale_bar(position=c(0.5, 0.14), width = 0.15, text.size = 2) +
    tm_layout(outer.margins = c(0,0,0,0), inner.margins = c(0,0,0,0))

tmap_save(inset_guil, filename='figs/map_components/guil.png', bg="transparent",
          dpi = 300, height = 2.5, width = 2.5, units = 'in')
# tmap_save(inset_guil, filename='figs/map_components/guil.pdf', bg="transparent",
#           dpi = 300, device = cairo_pdf, height = 2.5, width = 2.5, units = 'in')

neon_site <- filter(neon_sites, SiteID == 'MCRA')
inset_sites <- filter(ms_sites, domain == 'hjandrews')
bbox <- st_bbox(st_buffer(inset_sites, dist = 2e3))
basemap_mcra <- tmaptools::read_osm(bbox, type = 'stamen-terrain')
streams_layer <- get_osm_streams(st_as_sf(basemap_mcra))
inset_mcra <-
    tm_shape(basemap_mcra) + tm_rgb() +
    tm_shape(streams_layer) + tm_lines(col = 'steelblue3', lwd = 1) +
    tm_shape(inset_sites) + tm_symbols(shape = 19, col = 'purple4', size = 0.7, border.lwd = 2, alpha=0.5) +
    tm_shape(neon_site) + tm_symbols(shape = 19, col = neoncol, size = 0.7, border.lwd = 2, alpha=0.8) +
    # tm_compass(type="arrow", position=c("left", "top"), show.labels = 0, size = 1.5) +
    tm_scale_bar(position=c(0.05, 0.8), width = 0.25, text.size = 0.7)

tmap_save(inset_mcra, filename='figs/map_components/mcra.png', bg="transparent",
          dpi = 300, height = 2.5, width = 2.5, units = 'in')

# neon_site <- filter(neon_sites, SiteID == 'CARI')
# inset_sites <- bind_rows(usgs_sites$CARI, neon_site)
# bbox <- st_bbox(st_buffer(inset_sites, dist = 8e3))
# basemap_cari <- tmaptools::read_osm(bbox, type = 'stamen-terrain')
# streams_layer <- get_osm_streams(st_as_sf(basemap_cari))
# inset_cari <-
#     tm_shape(basemap_cari) + tm_rgb() +
#     tm_shape(streams_layer) + tm_lines(col = 'steelblue3', lwd = 1) +
#     tm_shape(usgs_sites$CARI) + tm_symbols(shape = 19, col = 'darkorange4', size = 0.9, border.lwd = 2, alpha=0.6) +
#     tm_shape(neon_site) + tm_symbols(shape = 19, col = neoncol, size = 0.9, border.lwd = 2, alpha=0.9) +
#     # tm_compass(type="arrow", position=c("right", "bottom"), show.labels = 0, size = 1) +
#     tm_scale_bar(position=c(0.55, 0.8), text.size = 2)
#
# tmap_save(inset_cari, filename='figs/map_components/cari.pdf', bg="transparent",
#           dpi = 300, device = cairo_pdf, height = 2.5, width = 2.5, units = 'in')

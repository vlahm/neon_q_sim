# Mike Vlah
# vlahm13@gmail.com
# last data retrieval dates given in comments below:
# last edit: 2023-04-12

library(raster)
library(tidyverse)
library(macrosheds)
library(ncdf4)
library(lubridate)
library(imputeTS)
library(foreach)
library(doParallel)
library(glue)
library(sf)

options(readr.show_progress = FALSE,
        readr.show_col_types = FALSE,
        timeout = 3000)

#pre-bundled in/out data available at: [**]
if(! exists('ts_plot')) source('src/00_helpers.R')
source('src/lstm_dungeon/camels_helpers.R')
if(! exists('ms_areas')) source('src/01_data_retrieval.R')

dir.create('in/lstm_data/attributes', recursive = TRUE, showWarnings = FALSE)
dir.create('in/lstm_data/time_series/', recursive = TRUE, showWarnings = FALSE)

## 1. retrieve additional input data for LSTM simulation ####

#NEON watershed boundary shapefiles; 2023-04-12
download.file('https://www.neonscience.org/sites/default/files/NEONAquaticWatershed_0.zip',
              'in/NEON/NEONAquaticWatershed_0.zip')

#CAMELS forcings and USGS streamflow (only extracting necessary files); 2023-04-12
if(! dir.exists('in/CAMELS/basin_dataset_public_v1p2')){

    dir.create('in/CAMELS', showWarnings = FALSE)
    download.file('https://gdex.ucar.edu/dataset/camels/file/basin_timeseries_v1p2_metForcing_obsFlow.zip',
                  destfile = 'in/CAMELS/basin_timeseries_v1p2_metForcing_obsFlow.zip')
    zf <- unzip('in/CAMELS/basin_timeseries_v1p2_metForcing_obsFlow.zip', list = TRUE)
    unzip('in/CAMELS/basin_timeseries_v1p2_metForcing_obsFlow.zip',
          files = c(grep('v1p2/basin_mean_forcing/daymet', zf$Name, value = TRUE),
                    grep('v1p2/usgs_streamflow', zf$Name, value = TRUE)),
          exdir = 'in/CAMELS')
    unlink('in/CAMELS/basin_timeseries_v1p2_metForcing_obsFlow.zip')
}

#CAMELS basin attributes; 2023-04-12
if(! file.exists('in/CAMELS/camels_attributes_v2.0/camels_hydro.txt')){

    dir.create('in/CAMELS/camels_attributes_v2.0', showWarnings = FALSE)
    download.file('https://gdex.ucar.edu/dataset/camels/file/camels_clim.txt',
                  destfile = 'in/CAMELS/camels_attributes_v2.0/camels_clim.txt')
    download.file('https://gdex.ucar.edu/dataset/camels/file/camels_soil.txt',
                  destfile = 'in/CAMELS/camels_attributes_v2.0/camels_soil.txt')
    download.file('https://gdex.ucar.edu/dataset/camels/file/camels_vege.txt',
                  destfile = 'in/CAMELS/camels_attributes_v2.0/camels_vege.txt')
    download.file('https://gdex.ucar.edu/dataset/camels/file/camels_topo.txt',
                  destfile = 'in/CAMELS/camels_attributes_v2.0/camels_topo.txt')
    download.file('https://gdex.ucar.edu/dataset/camels/file/camels_geol.txt',
                  destfile = 'in/CAMELS/camels_attributes_v2.0/camels_geol.txt')
    download.file('https://gdex.ucar.edu/dataset/camels/file/camels_hydro.txt',
                  destfile = 'in/CAMELS/camels_attributes_v2.0/camels_hydro.txt')
}

#CAMELS basin shapefiles; 2023-04-12
if(! file.exists('in/CAMELS/HCDN_nhru_final_671.shp')){
    download.file('https://gdex.ucar.edu/dataset/camels/file/basin_set_full_res.zip',
                  destfile = 'in/CAMELS/basin_set_full_res.zip')
    unzip('in/CAMELS/basin_set_full_res.zip', exdir = 'in/CAMELS/')
    unlink('in/CAMELS/basin_set_full_res.zip')
}

#Daymet model output (for SWE); 2023-04-12
if(! dir.exists('in/CAMELS/basin_dataset_public_v1p2/model_output_daymet')){
    download.file('https://gdex.ucar.edu/dataset/camels/file/basin_timeseries_v1p2_modelOutput_daymet.zip',
                  destfile = 'in/CAMELS/basin_timeseries_v1p2_modelOutput_daymet.zip')
    zf <- unzip('in/CAMELS/basin_timeseries_v1p2_modelOutput_daymet.zip', list = TRUE)
    unzip('in/CAMELS/basin_timeseries_v1p2_modelOutput_daymet.zip',
          files = grep('output/flow_timeseries/daymet/[0-9]{2}/[^\\.]+_model_output\\.txt',
                       zf$Name, value = TRUE),
          exdir = 'in/CAMELS/basin_dataset_public_v1p2/')
    unlink('in/CAMELS/basin_timeseries_v1p2_modelOutput_daymet.zip')
}

#Re-calibration of Priestley-Taylor coefficient apt=1.26 for ETo method
#(Priestley and Taylor, 1972) using ASCE-EWRI method (Allen et al., 2005)
#for short ref.crop [from Aschonitis et al. 2017]; retrieved 2023-04-12
if(! dir.exists('in/CAMELS/aptt1_30s')){
    download.file('https://hs.pangaea.de/model/ESRN-Database/30arc-sec/aptt1_30s.zip',
                  'in/CAMELS/priestly_taylor_recalibration.zip')
    unzip('in/CAMELS/priestly_taylor_recalibration.zip', exdir = 'in/CAMELS/')
    unlink('in/CAMELS/priestly_taylor_recalibration.zip')
}

#National Hydrologic Model (v1) data were provided to us for each CAMELS
#and MacroSheds site by the USGS. A minimal copy of this dataset (omitting
#files superfluous to this analysis) is included with our bundled input data
#at []

#retrieve MacroSheds core dataset and CAMELS-compliant MacroSheds attributes; 2023-04-12

if(! dir.exists('in/macrosheds/arctic')){

    ms_download_core_data(macrosheds_root = './in/macrosheds',
                          domains = 'all', version = '1.0')

    ms_download_ws_attr(macrosheds_root = './in/macrosheds',
                        dataset = 'CAMELS summaries', version = '1.0')
    ms_download_ws_attr(macrosheds_root = './in/macrosheds',
                        dataset = 'CAMELS Daymet forcings', version = '1.0')
}

ms_q <- ms_load_product(macrosheds_root = './in/macrosheds',
                        prodname = 'discharge', warn = FALSE) %>%
    filter(ms_status == 0)

ms_attr_static <- ms_load_product(macrosheds_root = './in/macrosheds',
                                  prodname = 'ws_attr_CAMELS_summaries',
                                  warn = FALSE)

ms_attr_dyn <- ms_load_product(macrosheds_root = './in/macrosheds',
                               prodname = 'ws_attr_CAMELS_Daymet_forcings',
                               warn = FALSE) %>%
    rename_with(function(x) str_extract(x, '^([^\\(]+)'))

#National Hydrologic Model v1 segment IDs
nhm_sites_ms <- read_csv('in/NHMv1/sites_with_segids.csv')
nhm_sites_camels <- read_csv('in/NHMv1/camels_gauge_info_with_segids.csv')


## X. recompute CAMELS forcings; build NEON forcings and attributes ####

source('src/lstm_dungeon/recompute_camels_climate.R', local = new.env())

## 2. prepare MacroSheds static attributes; write CSV ####

ms_attr_static <- read_csv('in/NEON/neon_attrs.csv') %>%
    mutate(domain = 'neon') %>%
    bind_rows(ms_attr_static) %>%
    filter(! is.na(domain), site_code != 'MC_ FLUME') %>%
    rename(q5 = Q5, q95 = Q95, baseflow_index = baseflow_index_landson)

neon_attr_static <- filter(ms_attr_static, domain == 'neon')

ms_attr_static <- bind_rows(ms_attr_static,
                            mutate(neon_attr_static,
                                   site_code = paste0(site_code, '_MANUALQ')))

ms_attr_static <- ms_attr_static %>%
    filter(site_code %in% nhm_sites_ms$site_code) %>%
    mutate(site_code = paste0('NHM_', site_code)) %>%
    bind_rows(ms_attr_static) %>%
    select(-domain, -network)

#if already deployed lstm_data_v1.zip or lstm_data_v2.zip
suppressWarnings(file.remove('in/lstm_data/attributes/ms_attributes.csv'))

write_csv(ms_attr_static, 'in/lstm_data/attributes/ms_attributes.csv')


## 3. prepare CAMELS static attributes; write CSV ####

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
                          col_types = 'cnnnnnncnn', delim = ';') %>%
    mutate(dom_land_cover = str_remove(dom_land_cover, '^ +'))
camels_hydro <- read_delim('in/CAMELS/camels_attributes_v2.0/camels_hydro.txt',
                           col_types = 'cnnnnnnnnnnnn', delim = ';')

camels_attr_static <- camels_geol %>%
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
    read_csv('in/CAMELS/recomputed_attributes/clim.csv'),
    read_csv('in/CAMELS/recomputed_attributes/soil.csv'),
    by = 'site_code'
)

replace_cols <- grep('site_code', colnames(camels_mod_attrs), value = TRUE, invert = TRUE)

camels_attr_static <- camels_attr_static %>%
    select(-all_of(replace_cols)) %>%
    full_join(camels_mod_attrs, by = 'site_code') %>%
    select(-all_of(setdiff(colnames(camels_attr_static), colnames(ms_attr_static)))) %>%
    filter(site_code %in% nhm_sites_camels$GAGE_ID) %>%
    mutate(site_code = paste0('NHM_', site_code)) %>%
    bind_rows(camels_attr_static)

write_csv(camels_attr_static, 'in/lstm_data/attributes/camels_attributes.csv')


## 4. prepare MacroSheds discharge ####

neon_q <- map_dfr(list.files('in/NEON/neon_continuous_Q', full.names = TRUE),
                  function(x){
                      read_csv(x) %>%
                          group_by(date = as.Date(datetime)) %>%
                          summarize(val = mean(discharge, na.rm = TRUE),
                                    site_code = first(site_code))
                  })
neon_q$val[is.na(neon_q$val)] <- NA_real_

ms_q <- ms_q %>%
    select(-any_of('year'), -val_err) %>%
    group_by(site_code) %>%
    filter(sum(! is.na(val)) > 201) %>% #must have n of at least LSTM sequence-length
    ungroup() %>%
    mutate(date = as.Date(datetime),
           val = ifelse(val < 0, 0, val)) %>%
    bind_rows(neon_q) %>%
    select(date, site_code, discharge = val) %>%
    arrange(site_code, date) %>%
    group_by(site_code) %>%
    tidyr::complete(date = seq(min(date), max(date), by = '1 day')) %>%
    ungroup()

#clean neon continuous Q according to Rhea et al. 2023
q_eval <- read_csv('in/neon_q_eval.csv') %>%
    filter(site %in% neon_sites) %>%
    group_by(site, year, month) %>%
    summarize(keep = all(final_qual %in% c('Tier1', 'Tier2')), #dropping <= tier 3
              .groups = 'drop') %>%
    rename(site_code = site)

neon_q_qc <- ms_q %>%
    filter(site_code %in% neon_sites) %>%
    mutate(year = year(date),
           month = month(date)) %>%
    left_join(q_eval, by = c('month', 'year', 'site_code')) %>%
    mutate(discharge = ifelse(! is.na(keep) & ! keep, NA_real_, discharge)) %>%
    select(-year, -month, -keep)

ms_q <- ms_q %>%
    filter(! site_code %in% neon_sites) %>%
    bind_rows(neon_q_qc) %>%
    group_split(site_code)

#include 365-day buffer of NAs before each discharge series, so the first Q
#value van be fully utilized/predicted
ms_q <- lapply(ms_q, function(x){
    first_q <- x$date[1]
    sc <- x$site_code[1]
    pre_buffer <- tibble(site_code = sc,
                         date = seq(first_q - 365, first_q - 1, 'day'),
                         discharge = NA_real_)
    x <- bind_rows(pre_buffer, x)
    return(x)
})

#additional buffer back to 2014 for NEON sites (to ensure prediction of field Q range)
ms_q <- lapply(ms_q, function(x){
    sc <- x$site_code[1]
    if(! sc %in% neon_sites) return(x)
    first_q <- x$date[1]
    pre_buffer <- tibble(site_code = sc,
                         date = seq(as.Date('2014-01-01'), first_q - 1, 'day'),
                         discharge = NA_real_)
    x <- bind_rows(pre_buffer, x)
    return(x)
})

ms_q <- map_dfr(ms_q, bind_rows)

#convert to specific Q (~= runoff)
ms_q <- ms_q %>%
    left_join(ms_areas, by = 'site_code') %>%
    filter(! is.na(ws_area_ha)) %>%
    #      mm/d        L/s         m^3/L   ha/m^2 mm/m   s/d     ha
    mutate(discharge = discharge * 0.001 * 1e-4 * 1000 * 86400 / ws_area_ha) %>%
    select(-ws_area_ha)

## 5. prepare MacroSheds forcings ####

ms_attr_dyn <- read_csv('in/macrosheds/neon_forcings.csv') %>%
    mutate(domain = 'neon') %>%
    bind_rows(ms_attr_dyn) %>%
    filter(! is.na(domain), site_code != 'MC_ FLUME')

neon_attr_dyn <- filter(ms_attr_dyn, domain == 'neon')

ms_attr_dyn <- bind_rows(ms_attr_dyn,
                         mutate(neon_attr_dyn,
                                site_code = paste0(site_code, '_MANUALQ')))

nhm_sites_ms <- read_csv('in/NHMv1/sites_with_segids.csv')
ms_attr_dyn <- ms_attr_dyn %>%
    filter(site_code %in% nhm_sites_ms$site_code) %>%
    mutate(site_code = paste0('NHM_', site_code)) %>%
    bind_rows(ms_attr_dyn) %>%
    select(-domain, -network)


## 6. merge MacroSheds discharge and forcings; write NetCDFs ####

ms_attr_dyn <- semi_join(ms_attr_dyn, ms_q, by = 'site_code')
ms_q <- semi_join(ms_q, ms_attr_dyn, by = 'site_code')

ms_attr_dyn$date <- as.double(as.Date(ms_attr_dyn$date))
ms_q$date <- as.double(ms_q$date) #numeric dates required for NetCDF format

ms_attr_dyn <- left_join(ms_q, ms_attr_dyn, by = c('date', 'site_code')) %>%
    group_split(site_code) %>%
    as.list()

#drop sites with no forcing data
ms_attr_dyn <- Filter(function(x){
    cmplt <- which(complete.cases(select(x, -discharge)))
    if(length(cmplt)) TRUE else FALSE
}, ms_attr_dyn)

#drop leading/trailing runs of missing forcings data
ms_attr_dyn <- lapply(ms_attr_dyn, function(x){
    cmplt <- which(complete.cases(select(x, -discharge)))
    x[min(cmplt):max(cmplt), ]
})

#impute forcings via linear interpolation with seasonal decomposition
for(i in seq_len(length(ms_attr_dyn))){

    ms_chunk <- ms_attr_dyn[[i]]
    site_code <- ms_chunk$site_code[1]

    na_proportions <- apply(select(ms_chunk, -date, -site_code),
                            MARGIN = 2,
                            FUN = function(x) sum(is.na(x)) / length(x))

    f_nas <- na_proportions[names(na_proportions) != 'discharge']

    #impute forcings NAs via linear interpolation with seasonal decomposition
    if(any(f_nas > 0)){

        ts_firstdate <- as.Date(ms_chunk$date, origin = '1970-01-01')[1]
        ts_startdoy <- as.numeric(format(ts_firstdate, format = '%j'))
        ts_startyr <- year(ts_firstdate)

        ts_lastdate <- as.Date(ms_chunk$date, origin = '1970-01-01')[nrow(ms_chunk)]
        ts_enddoy <- as.numeric(format(ts_lastdate, format = '%j'))
        ts_endyr <- year(ts_lastdate)

        ms_chunk$pet[is.infinite(ms_chunk$pet)] <- NA

        tryCatch({
            ms_attr_dyn[[i]] <- mutate(
                ms_chunk,
                across(-all_of(c('date', 'site_code', 'discharge')),
                       ~na_seadec(ts(.,
                                     start = c(!!ts_startyr, !!ts_startdoy),
                                     frequency = 365.25),
                                  maxgap = Inf)))
        }, warning = function(w){

            na_sums <- apply(select(ms_chunk, -date, -site_code, -discharge),
                             MARGIN = 2,
                             FUN = function(x) sum(is.na(x)))

            if(any(na_sums > 1)){
                stop(paste('must use a different interp algorithm for', i))
            }
        })
    }
}

#make sure all site datasets are sorted by date
ms_attr_dyn <- lapply(ms_attr_dyn, function(x) arrange(x, date))

ms_write_netcdf(df_list = ms_attr_dyn,
                path = 'in/lstm_data/time_series')


## 7. prepare NHM estimates for MacroSheds sites; write NetCDFs ####

nhm_ms_ts_splt <- list()
for(i in 1:nrow(nhm_sites_ms)){

    jobid <- pull(nhm_sites_ms[i, 'JobId'])
    segid <- as.character(pull(nhm_sites_ms[i, 'NHM_SEGID']))
    siteid <- pull(nhm_sites_ms[i, 'site_code'])
    ws_area_ha <- pull(ms_areas[ms_areas$site_code == siteid, 'ws_area_ha'])

    nhm_q <- read_csv(glue('in/NHMv1/nhm_output_ms/basin_{b}/model_output/seg_outflow.csv',
                           b = str_pad(jobid,
                                       width = 4,
                                       side = 'left',
                                       pad = '0')),
                      show_col_types = FALSE) %>%
        select(date = time, discharge_nhm = !!segid) %>%
        #      mm/d            ft^3/s          m^3/ft^3    ha/m^2 mm/m   s/d     ha
        mutate(discharge_nhm = discharge_nhm * 0.0283168 * 1e-4 * 1000 * 86400 / !!ws_area_ha)

    nhm_ms_ts_part <- ms_attr_dyn %>%
        filter(site_code == siteid) %>%
        left_join(nhm_q, by = 'date') %>%
        # filter(! is.na(discharge_nhm)) %>%
        arrange(date) %>%
        select(site_code, date, discharge = discharge_nhm, dayl, prcp, srad,
               swe, tmax, tmin, vp) %>%
        mutate(date = as.numeric(as.Date(date)),
               site_code = paste0('NHM_', site_code))

    not_shoulder_na <- ! cumprod(is.na(nhm_ms_ts_part$discharge)) &
        rev(! cumprod(is.na(rev(nhm_ms_ts_part$discharge))))
    nhm_ms_ts_part <- nhm_ms_ts_part[not_shoulder_na, ]

    nhm_ms_ts_splt[[i]] <- try({
        nhm_ms_ts_part %>%
            tidyr::complete(date = seq(min(date), max(date), 1)) %>%
            mutate(across(-all_of(c('date', 'site_code', 'discharge')),
                          ~na_interpolation(., maxgap = 3)))
    }, silent = TRUE)

    if(inherits(nhm_ms_ts_splt[[i]], 'try-error')) nhm_ms_ts_splt[[i]] = NA
}

nhm_ms_ts_splt[is.na(nhm_ms_ts_splt)] <- NULL

ms_write_netcdf(df_list = nhm_ms_ts_splt,
                path = 'in/lstm_data/time_series')

## 8. prepare CAMELS discharge and forcings; write NetCDFs ####

#locate usgs Q and daymet forcings
discharge_files <- list.files('in/CAMELS/basin_dataset_public_v1p2/usgs_streamflow',
                              full.names = TRUE, recursive = TRUE)

daymet_files_agg <- list.files('in/CAMELS/basin_dataset_public_v1p2/basin_mean_forcing/daymet',
                               full.names = TRUE, recursive = TRUE)
daymet_files_agg <- grep('011230\\*|208111310_forcing|344894205_forcing', #rm wonky files
                         daymet_files_agg, value = TRUE, invert = TRUE)

#the aggregate CAMELS dataset is missing all SWE information, so pulling it from the daymet model output
# daymet_dir <- '/home/mike/git/macrosheds/qa_experimentation/data/CAMELS/basin_dataset_public_v1p2/model_output_daymet/model_output/flow_timeseries/daymet'
daymet_dir <- 'in/CAMELS/basin_dataset_public_v1p2/model_output_daymet/model_output/flow_timeseries/daymet'
daymet_files <- list.files(daymet_dir, pattern = 'model_output',
                           full.names = TRUE, recursive = TRUE)

#dynamic PET will be a helpful forcing. not part of standard CAMELS dataset
pet_files <- list.files('in/CAMELS/camels_pet_isolate', full.names = TRUE)

#set up parallel processing
clst_type <- ifelse(.Platform$OS.type == 'windows', 'PSOCK', 'FORK')
ncores <- parallel::detectCores() %/% 1.5

nfls <- length(discharge_files)
nchunks <- nfls %/% ncores
remaind <- nfls %% ncores
megaloop_list <- list()
for(k in 0:(nchunks - 1)){
    megaloop_list[[k + 1]] <- 1:ncores + (k * ncores)
}
megaloop_list[[k + 2]] <- (nfls - remaind + 1):nfls

cmls_ts_splt <- list()
cnt <- 0
for(iterchunk in megaloop_list){

    cnt <- cnt + 1
    print(paste('chunk', cnt, 'of', length(megaloop_list)))

    #attempting to get around unpredictable mid-run failures
    #by registering and unregistering in an outer loop, parallelizing by chunk.
    #at any rate it helps with restarting the overall job
    clst <- parallel::makeCluster(spec = ncores, type = clst_type)
    doParallel::registerDoParallel(clst)

    cmls_ts_splt0 <- foreach(
        i = iterchunk,
        .combine = append) %dopar% {

            q_file <- discharge_files[i]
            basin_id <- str_match(q_file,
                                  'usgs_streamflow/[0-9]+/([0-9]+)_streamflow_qc.txt')[, 2]
            forcing_file <- grep(basin_id, daymet_files_agg, value = TRUE)
            pet_file <- grep(basin_id, pet_files, value = TRUE)

            #watershed area is also in basin_dataset_public_v1p2/basin_metadata/gauge_information.txt,
            #   but that's a different value (when converted km2 -> m2), and not the one NH uses
            ws_area_m2 <- read_lines(forcing_file, skip = 2, n_max = 1) %>%
                as.numeric()

            q_data_cmls <- read.fwf(
                discharge_files[i],
                widths = c(9, 4, 3, 3, 9, 2),
                header = FALSE,
                strip.white = TRUE,
                colClasses = 'character',
                col.names = c('basin_id', 'year', 'month', 'day', 'Q', 'flag')
            ) %>%
                as_tibble() %>%
                mutate(date = as.double(ymd(paste(year, month, day))),
                       Q = as.numeric(Q)) %>%
                select(date, discharge = Q) %>%
                mutate(discharge = ifelse(discharge == -999,
                                          NA_real_, #replace -999 with NA
                                          discharge),
                       #      mm/d        cfs         m^3/ft^3      mm/m   s/d     m^2
                       discharge = discharge * 0.028316846 * 1000 * 86400 / ws_area_m2)

            ffs = grep(basin_id, daymet_files, value = TRUE)
            #seems like some basins were processed multiple times under different HUCs. use the first HUC listed
            huc_code = unique(str_match(ffs, paste0(daymet_dir, '/([0-9]{2})/'))[, 2])[1]
            ffs = grep(paste0('/', huc_code, '/'), ffs, value = TRUE)

            forcing_iterations = list()
            for(k in seq_along(ffs)){

                fd = read.fwf(
                    ffs[k],
                    skip = 1,
                    widths = c(5, 3, 3, 3, 12, 12, 12, 12, 12, 12, 12, 12),
                    strip.white = TRUE,
                    colClasses = 'numeric'
                )

                colnames(fd) = c('year', 'month', 'day', 'hour', 'swe', 'prcp', 'raim',
                                 'tair', 'pet', 'et', 'mod_run', 'obs_run')

                forcing_iterations[[k]] = fd %>%
                    mutate(date = as.double(ymd(paste(year, month, day)))) %>%
                    select(-month, -day, -year, -hour)
            }

            if(length(forcing_iterations)){
                swe_mean <- purrr::map(forcing_iterations, ~.$swe) %>%
                    reduce(function(x, y) x + y) / 10
            }

            f_data_cmls <- read.fwf(
                forcing_file,
                skip = 4,
                widths = c(5, 3, 3, 100),
                header = FALSE,
                strip.white = TRUE,
                colClasses = 'numeric',
                col.names = c('year', 'month', 'day', 'hour', 'dayl', 'prcp', 'srad',
                              'swe', 'tmax', 'tmin', 'vp')
            ) %>%
                as_tibble() %>%
                mutate(date = as.double(ymd(paste(year, month, day))),
                       swe = ifelse(swe < 0, 0, swe)) %>%
                select(-month, -day, -year, -hour)
            print(table(f_data_cmls$swe))

            if(length(forcing_iterations)){
                f_data_cmls = left_join(f_data_cmls,
                                        tibble(swe = swe_mean, date = forcing_iterations[[1]]$date),
                                        copy = TRUE, by = 'date') %>%
                    select(-swe.x) %>% rename(swe = swe.y)
            } else {
                f_data_cmls$swe <- NA
            }

            if(length(pet_file)){

                cmls_pet <- read_csv(pet_file) %>%
                    mutate(date = as.double(date),
                           pet = ifelse(is.infinite(pet), NA, pet)) %>%
                    select(-site_code)

                ts_firstdate <- as.Date(cmls_pet$date, origin = '1970-01-01')[1]
                ts_startdoy <- as.numeric(format(ts_firstdate, format = '%j'))
                ts_startyr <- year(ts_firstdate)

                ts_lastdate <- as.Date(cmls_pet$date, origin = '1970-01-01')[nrow(cmls_pet)]
                ts_enddoy <- as.numeric(format(ts_lastdate, format = '%j'))
                ts_endyr <- year(ts_lastdate)

                tryCatch({
                    cmls_pet$pet = na_seadec(ts(cmls_pet$pet,
                                                start = c(!!ts_startyr, !!ts_startdoy),
                                                frequency = 365.25),
                                             maxgap = Inf)
                }, warning = function(w){

                    if(any(is.na(cmls_pet$pet))){
                        stop(paste('must use a different interp algorithm for', i))
                    }
                })

                cmls_data <- left_join(q_data_cmls, f_data_cmls, by = 'date') %>%
                    left_join(cmls_pet) %>%
                    mutate(site_code = !!basin_id,
                           q_source = 'true')

            } else {

                cmls_data <- left_join(q_data_cmls, f_data_cmls, by = 'date') %>%
                    mutate(site_code = !!basin_id,
                           q_source = 'true')
            }

            return(list(cmls_data))
        }

    cmls_ts_splt <- append(cmls_ts_splt, cmls_ts_splt0)

    parallel::stopCluster(clst)
}

# saveRDS(cmls_ts_splt, '~/Desktop/cmls_ts_splt.rds')
# cmls_ts_splt <- readRDS('~/Desktop/cmls_ts_splt.rds')

ms_write_netcdf(df_list = cmls_ts_splt,
                path = 'in/lstm_data/time_series')

## 9. prepare NHM estimates for CAMELS sites; write NetCDFs ####

#each NHMv1 gauge output file includes estimated Q for all upstream reaches, but
#nowhere does the output indicate which segid corresponds to the gauge itself.
#so, to find out, we get the location of the gage from the CAMELS dataset,
#then find the nearest feature within the NHM. verified on 01013500, and
#assumed to work fine for the other 530 gages.

camels_segids <- st_read('in/NHMv1/GF_nat_reg.gdb',
                         layer = 'nsegmentNationalIdentifier') %>%
    st_transform(crs = 4326)

camels_gauge_info_with_segids <- 'in/NHMv1/camels_gauge_info_with_segids.csv'
if(file.exists(camels_gauge_info_with_segids)){
    camels_gauge_info <- read_csv(camels_gauge_info_with_segids)
} else {

    camels_gauge_info <- read_tsv(
        'in/CAMELS/basin_dataset_public_v1p2/basin_metadata/gauge_information.txt',
        col_names = c('HUC_02', 'GAGE_ID', 'GAGE_NAME', 'LAT', 'LONG', 'DRAINAGE_AREA'),
        skip = 1
    )

    camels_gauge_info$segid <- NA
    for(i in seq_along(cmls_ts_splt)){

        cmls_gageid <- cmls_ts_splt[[i]]$site_code[1]
        cmls_gage_ind <- which(camels_gauge_info$GAGE_ID == cmls_gageid)
        cmls_lat <- camels_gauge_info$LAT[cmls_gage_ind]
        cmls_long <- camels_gauge_info$LONG[cmls_gage_ind]
        cmls_site <- tibble(x = cmls_lat,
                            y = cmls_long) %>%
            sf::st_as_sf(coords = c('y', 'x'),
                         crs = 4326)

        nearest_ind <- st_nearest_feature(cmls_site, camels_segids)
        cmls_segid <- camels_segids[nearest_ind, ]$seg_id_nat

        camels_gauge_info$segid[cmls_gage_ind] <- cmls_segid
    }

    write_csv(camels_gauge_info, camels_gauge_info_with_segids)
}

nhm_cmls_ts_splt <- cmls_ts_splt
for(i in seq_along(cmls_ts_splt)){

    if(i %% 50 == 0) print(paste0(i, '/', length(cmls_ts_splt)))

    cmls_ts_part <- cmls_ts_splt[[i]]
    cmls_gageid <- cmls_ts_part$site_code[1]
    cmls_gage_ind <- which(camels_gauge_info$GAGE_ID == cmls_gageid)
    cmls_segid <- as.character(camels_gauge_info$segid[cmls_gage_ind])

    forcing_file <- grep(cmls_gageid, daymet_files_agg, value = TRUE)
    cmls_ws_area_m2 <- read_lines(forcing_file, skip = 2, n_max = 1) %>%
        as.numeric()

    nhm_q <- try({
        read_csv(glue('in/NHMv1/nhm_output_camels/{gid}/model_output/seg_outflow.csv',
                      gid = cmls_gageid),
                 show_col_types = FALSE) %>%
            select(date = time,
                   discharge_nhm = !!cmls_segid) %>%
            #      mm/d            cfs             m^3/ft^3      mm/m   s/d     m^2
            mutate(discharge_nhm = discharge_nhm * 0.028316846 * 1000 * 86400 / cmls_ws_area_m2)
    }, silent = TRUE)

    if(inherits(nhm_q, 'try-error')){
        # print(paste('error in', cmls_gageid))
        next
    }

    cmls_ts_part <- cmls_ts_part %>%
        mutate(date_as_date = as.Date(date, origin = '1970-01-01')) %>%
        left_join(nhm_q, by = c('date_as_date' = 'date')) %>%
        select(site_code, date, discharge = discharge_nhm, dayl, prcp, srad,
               swe, tmax, tmin, vp, any_of('pet')) %>%
        mutate(q_source = 'NHM',
               site_code = paste0('NHM_', site_code)) %>%
        arrange(date)

    not_shoulder_na <- ! cumprod(is.na(cmls_ts_part$discharge)) &
        rev(! cumprod(is.na(rev(cmls_ts_part$discharge))))
    cmls_ts_part <- cmls_ts_part[not_shoulder_na, ]

    nhm_cmls_ts_splt[[i]] <- cmls_ts_part
}

ms_write_netcdf(df_list = nhm_cmls_ts_splt,
                path = 'in/lstm_data/time_series')

## 10. create additional NetCDFs with field Q for each NEON site ####

ncdfs <- list.files('in/lstm_data/time_series', full.names = TRUE)

for(i in seq_along(neon_sites)){

    s <- neon_sites[i]
    f <- grep(paste0('(?<!NHM_)', s, '\\.nc'), ncdfs, value = TRUE, perl = TRUE)
    if(! length(f)) next
    f_newname <- gsub('\\.nc$', '_MANUALQ.nc', f)
    file.copy(f, f_newname, overwrite = TRUE)

    nc_con <- ncdf4::nc_open(f_newname, write = TRUE)
    nc_q <- ncdf4::ncvar_get(nc_con, 'discharge')
    nc_date <- ncdf4::ncvar_get(nc_con, 'date')

    field_q <- read_csv(paste0('in/neon_field_Q/', s, '.csv')) %>%
        mutate(date = as.Date(datetime)) %>%
        left_join(ms_areas, by = 'site_code') %>%
        filter(! is.na(ws_area_ha)) %>%
        #      mm/d        L/s         m^3/L   ha/m^2 mm/m   s/d     ha
        mutate(discharge = discharge * 0.001 * 1e-4 * 1000 * 86400 / ws_area_ha) %>%
        select(-ws_area_ha, -site_code, -datetime) %>%
        distinct(date, .keep_all = TRUE)

    nc_q_new <- left_join(tibble(date = nc_date),
                          mutate(field_q, date = as.numeric(date)), by = 'date')
    nc_q_new$discharge[is.na(nc_q_new$discharge)] <- NaN

    if(nrow(nc_q_new) != length(nc_q)) stop('!!')

    ncvar_put(nc_con, 'discharge', nc_q_new$discharge)
    nc_close(nc_con)
}


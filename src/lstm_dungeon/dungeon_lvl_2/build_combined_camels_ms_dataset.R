library(segmented)
# library(splines)
# library(lme4)
library(feather)
library(nhdplusTools)
library(data.table)
library(imputeTS)
library(lubridate)
library(glue)
library(dygraphs)
library(xts)
library(errors)
library(ncdf4)
library(foreach)
library(doParallel)
library(purrr)
library(sf)
library(ggplot2)
library(gridExtra)
library(htmlwidgets)
library(DAAG)
library(forecast)
library(tidyverse)

reticulate::use_condaenv('nh2')
xr <- reticulate::import("xarray")
pd <- reticulate::import("pandas")
np <- reticulate::import("numpy")

#TODO: update list of arctic sites to drop, once dataset updates
#switch on addition of USGS Q record for TECR if we want to predict that site better in the past (~1980)
#once NHC is in ms for reals, remove custom NHC code below

options(readr.show_progress = FALSE)
options(readr.show_col_types = FALSE)

#for macrosheds sites that are commensurate with the CAMELS dataset,
#   this script formats their timeseries data (discharge and forcings)
#   and static attributes according to the documentation for a "generic" dataset
#   in the neuralhydrology package. It saves the formatted data to
#   /home/mike/git/macrosheds/qa_experimentation/data/CAMELS_macrosheds_combined.
#   It does the same for CAMELS data, creating a combined dataset of CAMELS and
#   macrosheds basins that can be used to train neural nets.

setwd('/home/mike/git/macrosheds/qa_experimentation/data')

### 1 setup ####

ms_sitelist_loc <- '/home/mike/git/macrosheds/qa_experimentation/data/informative_stuff/site_data.csv'

conf <- jsonlite::fromJSON('../../data_acquisition/config.json',
                           simplifyDataFrame = FALSE)

site_data <- googlesheets4::read_sheet(
    conf$site_data_gsheet,
    na = c('', 'NA'),
    col_types = 'ccccccccnnnnnccccc'
)

# write_csv(site_data, ms_sitelist_loc)

neon_sites <- c('FLNT', 'WALK', 'WLOU', 'LEWI', 'HOPB', 'MART', 'CUPE',
                'BIGC', 'REDB', 'COMO', 'BLDE', 'BLWA', 'PRIN', 'POSE',
                'MCDI', 'MCRA', 'CARI', 'ARIK', 'GUIL', 'KING',
                'MAYF', 'LECO', 'BLUE', 'SYCA', 'TECR', 'TOMB', 'OKSR')

neon_neighbors = read_csv('informative_stuff/nearest_neon_neighbors.csv')
neon_neighbors = distinct(neon_neighbors, nearest_neighbor, .keep_all = TRUE)

ms_areas <- site_data %>%
    filter(site_type == 'stream_gauge',
           ! is.na(ws_area_ha)) %>%
    select(site_code, ws_area_ha)

#keepers: Trevor_Creek_Main, Kuparuk_River_1, Oksrukuyik_Creek_1, Roche_Moutonnee_Creek_Main, Toolik_Inlet_Main
arctic_sites_to_drop <- c('Kuparuk_River_-0.1', 'Kuparuk_River_-0.177', 'Kuparuk_River_-0.3','Kuparuk_River_-0.47',
                          'Kuparuk_River_-0.7', 'Kuparuk_River_0.3', 'Kuparuk_River_0.56', 'Kuparuk_River_0.74',
                          'Kuparuk_River_1.39', 'Kuparuk_River_1.5', 'Kuparuk_River_1.8', 'Kuparuk_River_0', 'Kuparuk_River_2.5',
                          'Kuparuk_River_2', 'Kuparuk_River_3', 'Kuparuk_River_4.1', 'Kuparuk_River_4', 'Oksrukuyik_Creek_-0.1',
                          'Oksrukuyik_Creek_-0.3', 'Oksrukuyik_Creek_-0.7', 'Oksrukuyik_Creek_-1.2', 'Oksrukuyik_Creek_0.23',
                          'Oksrukuyik_Creek_0.48', 'Oksrukuyik_Creek_1.06', 'Oksrukuyik_Creek_1.37', 'Oksrukuyik_Creek_1.7',
                          'Oksrukuyik_Creek_2.5', 'Oksrukuyik_Creek_2.7',
                          'Kuparuk_River_-0', 'Oksrukuyik_Creek_-0', 'Oksrukuyik_Creek_-1',
                          'Oksrukuyik_Creek_0', 'Oksrukuyik_Creek_2')

q_source_map <- list(
    'true' = 0, #measured Q data
    'NHM' = 1 #process model estimates
)

ms_version <- 1
ms_domains_included <- list.files('ms_in_camels_format')
ms_domains_included <- ms_domains_included[! grepl('camels', ms_domains_included)]

dir.create('CAMELS_macrosheds_combined/time_series',
           showWarnings = FALSE,
           recursive = TRUE)
dir.create('CAMELS_macrosheds_combined/attributes',
           showWarnings = FALSE)

source('../../data_acquisition/src/output_dataset_convenience_functions/load_product.R')

ms_write_netcdf <- function(df_list, path){

    vars_units0 <- c(discharge = 'cms', #if you change cms, change section 6b below
                     dayl = 's', prcp = 'mm/day',
                    srad = 'W/m2', swe = 'mm', tmax = 'C', tmin = 'C', vp = 'Pa',
                    pet = 'mm', #units?
                    q_source = '')

    for(i in seq_along(df_list)){

        d_chunk <- df_list[[i]]

        if(! 'pet' %in% colnames(d_chunk)){
            vars_units <- vars_units0[names(vars_units0) != 'pet']
        } else {
            vars_units <- vars_units0
        }

        ncdf_loc <- glue('{pth}/{s}.nc',
                         pth = path,
                         s = d_chunk$site_code[1])

        dim_time <- ncdim_def('date',
                              units = 'days since 1970-01-01 00:00',
                              calendar = 'standard',
                              vals = d_chunk$date)

        ncdf_varlist <- list()
        for(j in seq_along(vars_units)){

            ncdf_varlist[[j]] <- ncvar_def(name = names(vars_units)[j],
                                           units = unname(vars_units)[j],
                                           dim = list(dim_time),
                                           missval = NULL)
        }

        con <- nc_create(filename = ncdf_loc,
                         vars = ncdf_varlist)

        for(j in seq_along(vars_units)){

            varname <- ncdf_varlist[[j]]$name
            ncvar_put(con, varname, d_chunk[[varname]])
        }

        nc_close(con)
    }
}

nhm_segids <- read_csv('site_id_lists/sites_COMIDS_SEGIDS.csv') %>%
    select(domain, site_code, NHM_SEGID) %>%
    filter(! NHM_SEGID %in% c('non-NHDPlus', 'too small'))

nhm_segids <- read_csv('NHMv1/nhm_output_ms/seg_id_nats.csv',
                       col_types = 'cc') %>%
    right_join(nhm_segids, by = 'NHM_SEGID')

Mode <- function(x, na.rm = TRUE){

    if(na.rm){
        x <- na.omit(x)
    }

    ux <- unique(x)
    mode_out <- ux[which.max(tabulate(match(x, ux)))]
    return(mode_out)

}

neglog <- function(x){
    sgn <- sign(x)
    neglog_ <- log(abs(x) + 1) * sgn
    return(neglog_)
}
inv_neglog <- function(x){
    sgn <- sign(x)
    inv_neglog_ <- (exp(x) - 1) * sgn
    return(inv_neglog_)
}
boxcox_write = function(vec, neon_site, write = TRUE){

    lambda = forecast::BoxCox.lambda(vec)
    vec_trans = forecast::BoxCox(vec, lambda = lambda)

    if(write){
        write_lines(lambda, glue('../imputation/data/boxcox_lambdas/{neon_site}.txt'))
    }

    return(vec_trans)
}
boxcox_inv_read = function(vec_trans, neon_site){

    lambda = write_lines(glue('../imputation/data/boxcox_lambdas/{neon_site}.txt'))
    vec_backtrans = forecast::InvBoxCox(vec_trans, lambda = lambda)

    return(vec_backtrans)
}

# functions for estimating reach-distance between points
get_fl <- function(hl, net){
    filter(net,
           REACHCODE == hl$REACHCODE &
               FromMeas < hl$REACH_meas &
               ToMeas > hl$REACH_meas)
}
get_lengths <- function(hl, net = NULL, fl = NULL){

    if(is.null(fl)) {
        fl <- get_fl(hl, net)
    }

    meas <- rescale_measures(measure = hl$REACH_meas,
                             from = fl$FromMeas,
                             to = fl$ToMeas) / 100

    list(dn = fl$LENGTHKM * meas,
         up = fl$LENGTHKM * 1 - meas)
}
get_distance <- function(hl1, hl2, net, check_connected = TRUE){
    fl1 <- get_fl(hl1, net)

    fl2 <- get_fl(hl2, net)

    x <- list(upfl = fl2, uphl = hl2, dnfl = fl1, dnhl = hl1)

    if(fl1$Hydroseq > fl2$Hydroseq) x <- list(upfl = fl1, uphl = hl1,
                                              dnfl = fl2, dnhl = hl2)

    if(check_connected) {
        if(!x$dnfl$COMID %in% get_DM(net, x$upfl$COMID))
            stop("hydro locations not connected.")
    }

    x$upfl$Pathlength -
        x$dnfl$Pathlength -
        get_lengths(x$uphl, fl = x$upfl)$up -
        get_lengths(x$dnhl, fl = x$dnfl)$dn
}

#lm funcs
get_neon_field_discharge = function(neon_sites){

    for(i in seq_along(neon_sites)){

        s = neon_sites[i]
        print(s)

        #field discharge measurements
        qd = neonUtilities::loadByProduct('DP1.20048.001', site = s, check.size = FALSE)

        q1 = q2 = tibble()
        try({
            q1 = select(qd$dsc_fieldDataADCP, discharge = totalDischarge, date = startDate, totalDischargeUnits)
        }, silent = TRUE)
        try({
            q2 = select(qd$dsc_fieldData, discharge = totalDischarge, date = startDate, totalDischargeUnits)
        }, silent = TRUE)
        if(nrow(q1) && nrow(q2)){
            q = bind_rows(q1, q2)
        } else if(nrow(q1)){
            q = q1
        } else {
            q = q2
        }

        if(any(! q$totalDischargeUnits %in% c('cubicMetersPerSecond', 'litersPerSecond'))) stop('new unit detected. account for this.')

        q = q %>%
            mutate(discharge = ifelse(totalDischargeUnits == 'cubicMetersPerSecond', discharge * 1000, discharge),
                   site_code = s) %>%
            rename(datetime = date) %>%
            select(-totalDischargeUnits) %>%
            as_tibble()

        write_csv(q, glue('../imputation/data/neon_field_Q/{s}.csv'))
    }
}
get_neon_inst_discharge = function(neon_sites){

    for(i in seq_along(neon_sites)){

        s = neon_sites[i]
        print(s)

        #continuous discharge measurements
        qd = neonUtilities::loadByProduct('DP4.00130.001', site = s, check.size = FALSE)

        q1 = q2 = tibble()
        try({
            q1 = select(qd$csd_continuousDischarge, discharge = maxpostDischarge, datetime = endDate)
        }, silent = TRUE)
        try({
            q2 = select(qd$csd_continuousDischargeUSGS, discharge = usgsDischarge, datetime = endDate)
        }, silent = TRUE)
        if(nrow(q1) && nrow(q2)){
            q = bind_rows(q1, q2)
        } else if(nrow(q1)){
            q = q1
        } else {
            q = q2
        }

        q = as_tibble(q) %>%
            mutate(site_code = s)

        write_csv(q, glue('../imputation/data/neon_continuous_Q/{s}.csv'))
    }
}
assemble_q_lm_df = function(neon_site, nearby_usgs_gages = NULL, ms_Q_data = NULL,
                            datetime_snapdist_hrs = 12, overwrite = TRUE, scale_q_by_area = T){

    #ms_Q_data: data.frame with posixct datetime column, and Q columns. Each Q column must
    #have its site_code as header. transformed columns will be created

    neon_q_manual = read_csv(glue('../imputation/data/neon_field_Q/{neon_site}.csv')) %>%
        mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
        rename(discharge_manual = discharge) %>%
        distinct(datetime, .keep_all = TRUE)
    earliest_date = as.character(date(min(neon_q_manual$datetime)))

    if(scale_q_by_area){
        wsa = filter(ms_areas, site_code == neon_site) %>% pull(ws_area_ha)
        neon_q_manual$discharge_manual = neon_q_manual$discharge_manual / wsa * 1000
    }

    if(! file.exists(glue('../imputation/data/usgs_ms_q/{neon_site}.csv')) || overwrite){

        site_nearbyA = NULL
        if(! is.null(nearby_usgs_gages)){

            usgsq = dataRetrieval::readNWISdata(
                sites = nearby_usgs_gages,
                service = 'iv', #instantaneous values
                parameterCd = '00060', #discharge (cfs)
                startDate = earliest_date,
                endDate = '2099-01-01'
            )

            if(! nrow(usgsq)) stop('no instantaneous Q?')

            # if(any(! usgsq$X_00060_00000_cd %in% c('A', 'P', 'A e', 'A R'))) stop()
            if(any(grepl('>', usgsq$X_00060_00000_cd))) warning(glue('">" detected in {neon_site}'))
            if(any(grepl('<', usgsq$X_00060_00000_cd))) warning(glue('"<" detected in {neon_site}'))
            if(any(usgsq$tz_cd != 'UTC')) stop()

            site_nearbyA = usgsq %>%
                as_tibble() %>%
                select(site_no, datetime = dateTime, discharge = X_00060_00000) %>%
                pivot_wider(names_from = site_no, values_from = discharge)

            if(scale_q_by_area){

                for(g in nearby_usgs_gages){

                    # ddd = nhdplusTools::get_nldi_basin(list(featureSource = "nwissite", featureID = "USGS-50128907"),
                    #                                    simplify = FALSE, split = TRUE)
                    # ga = st_area(ddd)
                    # units(ga) = 'ha'

                    if(g == '06190540'){
                        site_nearbyA[[g]] = site_nearbyA[[g]] / 51089.73 * 1000
                        next
                    }

                    nwissite = dataRetrieval::readNWISsite(g)

                    if(is.null(nwissite)) stop('cannot find nwis watershed area')

                    if(is.na(nwissite$contrib_drain_area_va)){
                        wsa = nwissite$drain_area_va * 258.999 #mi^2 -> ha
                    } else {
                        wsa = nwissite$contrib_drain_area_va * 258.999 #mi^2 -> ha
                    }

                    if(is.na(wsa)) stop(paste('need to look up area for', g))

                    site_nearbyA[[g]] = site_nearbyA[[g]] / wsa * 1000
                }
            }

            site_nearbyA = site_nearbyA %>%
                mutate(across(matches('^[0-9]+$'), ~ . * 28.3168)) %>%
                mutate(across(matches('^[0-9]+$'),
                              neglog,
                              # ~boxcox_write(., !!neon_site, write = FALSE),
                              .names = '{.col}_log')) %>%
                arrange(datetime)
        }

        if(! is.null(ms_Q_data)){

            ms_sites = grep('datetime', colnames(ms_Q_data), invert = TRUE, value = TRUE)

            if(scale_q_by_area){

                for(g in ms_sites){

                    wsa = filter(ms_areas, site_code == g) %>% pull(ws_area_ha)
                    ms_Q_data[[g]] = ms_Q_data[[g]] / wsa * 1000
                }
            }

            ms_Q_data = ms_Q_data %>%
                mutate(across(-datetime, neglog, .names = '{.col}_log'))
        }

        if(! is.null(site_nearbyA) && ! is.null(ms_Q_data)){
            site_nearby = full_join(site_nearbyA, ms_Q_data, by = c('datetime'))
        } else if(! is.null(site_nearbyA)){
            site_nearby = site_nearbyA
        } else  {
            site_nearby = ms_Q_data
        }

        write_csv(site_nearby, glue('../imputation/data/usgs_ms_q/{neon_site}.csv'))

    } else {
        site_nearby = read_csv(glue('../imputation/data/usgs_ms_q/{neon_site}.csv'))
    }

    #rolling join usgs data to neon field measurements by datetime

    x = rename(neon_q_manual, datetime_x = datetime) %>% as.data.table()
    y = rename(site_nearby, datetime_y = datetime) %>% as.data.table()

    rollmax = 60 * 60 * datetime_snapdist_hrs #join up to this number of hrs off
    x[, `:=` (datetime_min = datetime_x - rollmax,
              datetime_max = datetime_x + rollmax)]
    y[, `:=` (datetime_y_orig = datetime_y)] #this datetime col will be dropped

    #join x rows to y if y's datetime falls within the x range
    joined = y[x, on = .(datetime_y <= datetime_max,
                         datetime_y >= datetime_min)]
    joined = na.omit(joined, cols = 'datetime_y_orig') #drop rows without matches

    # joined = data.table(x=as.Date(c('2019-12-31', '2020-01-01', '2020-01-02',
    #                                 '2019-12-31', '2020-01-01', '2020-01-02',
    #                                 '2019-12-31', '2020-01-01', '2020-01-02')),
    #                     y=as.Date(c('2020-01-01', '2020-02-01', '2020-03-01',
    #                                 '2020-01-01', '2020-01-02', '2020-01-03',
    #                                 '2020-03-01', '2020-02-01', '2020-01-02')),
    #            a = c(NA, 1, 2, NA,NA,NA, NA,NA,3), b=c(1,1,1,2,2,2,3,3,3))
    # joined[, `:=` (datetime_match_diff = abs(x - y))]
    # joined[order(datetime_match_diff), lapply(.SD, function(xx) first(na.omit(xx))), by = x]

    #for any datetimes in x or y that were matched more than once, keep only
    #the nearest non-NA match by time
    joined[, `:=` (datetime_match_diff = abs(datetime_x - datetime_y_orig))]
    joined = joined[order(datetime_match_diff),
                    lapply(.SD, function(z) first(na.omit(z))),
                    by = datetime_x]

    joined[, c('datetime_y', 'datetime_y.1', 'datetime_y_orig', 'datetime_match_diff') := NULL]
    setnames(joined, 'datetime_x', 'datetime')

    joined = joined %>%
        as_tibble() %>%
        arrange(datetime) %>%
        rename(discharge = discharge_manual) %>%
        mutate(season = factor(lubridate::quarter(datetime)),
               discharge_log = neglog(discharge)) %>%
               # discharge_log = boxcox_write(discharge, !!neon_site)) %>%
        select(site_code, datetime, discharge, discharge_log, everything())

    return(joined)
}
assemble_q_lm_df_daily = function(neon_site, ms_Q_data = NULL, scale_q_by_area = T,
                                  overwrite = TRUE){

    #ms_Q_data: data.frame with date column, and Q columns. Each Q column must
    #have its site_code as header. transformed columns will be created

    neon_q_manual = read_csv(glue('../imputation/data/neon_field_Q/{neon_site}.csv')) %>%
        mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
        rename(discharge_manual = discharge) %>%
        distinct(datetime, .keep_all = TRUE) %>%
        mutate(date = as.Date(datetime)) %>%
        select(-datetime)

    if(scale_q_by_area){
        wsa = filter(ms_areas, site_code == neon_site) %>% pull(ws_area_ha)
        neon_q_manual$discharge_manual = neon_q_manual$discharge_manual / wsa * 1000
    }

    if(! file.exists(glue('../imputation/data/usgs_ms_q/{neon_site}.csv')) || overwrite){

        ms_sites = grep('date', colnames(ms_Q_data), invert = TRUE, value = TRUE)

        if(scale_q_by_area){

            for(g in ms_sites){

                wsa = filter(ms_areas, site_code == g) %>% pull(ws_area_ha)
                ms_Q_data[[g]] = ms_Q_data[[g]] / wsa * 1000
            }
        }

        ms_Q_data = ms_Q_data %>%
            mutate(across(-date, neglog, .names = '{.col}_log'))

        write_csv(ms_Q_data, glue('../imputation/data/usgs_ms_q/{neon_site}.csv'))

    } else {
        ms_Q_data = read_csv(glue('../imputation/data/usgs_ms_q/{neon_site}.csv'))
    }

    joined = left_join(neon_q_manual, ms_Q_data, by = 'date') %>%
        arrange(date) %>%
        rename(discharge = discharge_manual) %>%
        mutate(season = factor(lubridate::quarter(date)),
               discharge_log = neglog(discharge)) %>%
        select(site_code, date, discharge, discharge_log, everything())

    return(joined)
}
assemble_q_lm_df_daily_ms = function(site, ms_Q_data = NULL, nearby_usgs_gages = NULL,
                                     scale_q_by_area = T, overwrite = TRUE){

    #ms_Q_data: data.frame with date column, and Q columns. Each Q column must
    #have its site_code as header. transformed columns will be created

    if(! file.exists(glue('../imputation/data/usgs_ms_q/{site}.csv')) || overwrite){

        usgsq = dataRetrieval::readNWISdata(
            sites = nearby_usgs_gages,
            service = 'dv', #daily values
            parameterCd = '00060', #discharge (cfs)
            startDate = '1900-01-01',
            endDate = '2099-01-01'
        )

        if(! nrow(usgsq)) stop('no daily Q?')

        if(any(grepl('>', usgsq$X_00060_00000_cd))) warning(glue('">" detected in {neon_site}'))
        if(any(grepl('<', usgsq$X_00060_00000_cd))) warning(glue('"<" detected in {neon_site}'))
        if(any(usgsq$tz_cd != 'UTC')) stop()

        site_nearbyA = usgsq %>%
            as_tibble() %>%
            select(site_no, datetime = dateTime, discharge = X_00060_00003) %>%
            pivot_wider(names_from = site_no, values_from = discharge)

        if(scale_q_by_area){

            for(g in nearby_usgs_gages){

                nwissite = dataRetrieval::readNWISsite(g)

                if(is.null(nwissite)) stop('cannot find nwis watershed area')

                if(is.na(nwissite$contrib_drain_area_va)){
                    wsa = nwissite$drain_area_va * 258.999 #mi^2 -> ha
                } else {
                    wsa = nwissite$contrib_drain_area_va * 258.999 #mi^2 -> ha
                }

                if(is.na(wsa)) stop(paste('need to look up area for', g))

                site_nearbyA[[g]] = site_nearbyA[[g]] / wsa * 1000
            }
        }

        site_nearbyA = site_nearbyA %>%
            mutate(across(matches('^[0-9]+$'), ~ . * 28.3168)) %>%
            mutate(across(matches('^[0-9]+$'),
                          neglog,
                          .names = '{.col}_log')) %>%
            mutate(date = as.Date(datetime)) %>%
            select(-datetime) %>%
            arrange(date)

        if(scale_q_by_area){
            wsa = filter(ms_areas, site_code == site) %>% pull(ws_area_ha)
            ms_Q_data[[site]] = ms_Q_data[[site]] / wsa * 1000
        }

        ms_Q_data = ms_Q_data %>%
            mutate(across(-date, neglog, .names = '{.col}_log'))

        if(! is.null(site_nearbyA) && ! is.null(ms_Q_data)){
            site_nearby = full_join(site_nearbyA, ms_Q_data, by = c('date'))
        } else if(! is.null(site_nearbyA)){
            site_nearby = site_nearbyA
        } else  {
            site_nearby = ms_Q_data
        }

        site_nearby = rename_with(site_nearby, function(x) sub(site, 'discharge', x)) %>%
            mutate(season = factor(lubridate::quarter(date)),
                                   site_code = !!site) %>%
            select(site_code, date, discharge, discharge_log, everything())

        write_csv(site_nearby, glue('../imputation/data/usgs_ms_q/{site}.csv'))

    } else {
        site_nearby = read_csv(glue('../imputation/data/usgs_ms_q/{site}.csv'))
    }

    return(site_nearby)
}
generate_nested_formulae = function(full_spec, d, interactions = TRUE, through_origin = TRUE, min_points_per_param = 15,
                                    max_interaction = Inf){

    #this one allows all interactions if interactions =+ TRUE

    trms = rownames(attributes(terms(full_spec))$factors)

    y = trms[1]
    y_nobacktick = gsub('`', '', y)
    trms = trms[! trms == y]
    ntrms = length(trms)

    indeps = list()
    for(i in seq_len(ntrms)){
        indeps = append(indeps, combn(trms, i, simplify = FALSE))
    }

    if(interactions){
        int_trms = Filter(function(x) length(x) > 1, indeps) %>%
            map(~paste(., collapse = ':')) %>%
            unlist()
    }

    max_params = sum(! is.na(d[[y_nobacktick]])) %/% min_points_per_param #https://statisticsbyjim.com/regression/overfitting-regression-models/
    indeps = Filter(function(x) length(x) <= max_params, indeps)

    mods = list()
    modind = 0
    for(i in seq_along(indeps)){

        modind = modind + 1
        mods[[modind]] = as.formula(glue('{y} ~ {paste(indeps[[i]], collapse = "+")}'))

        if(through_origin){
            modind = modind + 1
            mods[[modind]] = as.formula(glue('{y} ~ 0 + {paste(indeps[[i]], collapse = "+")}'))
        }

        if(interactions){

            valid_ints_for_indep = map(int_trms, ~str_split(., ':')[[1]]) %>%
                map(~all(. %in% indeps[[i]])) %>%
                unlist()
            ints_to_tack_on = int_trms[valid_ints_for_indep]
            ints_to_tack_on = Filter(function(x) length(strsplit(x, ':')[[1]]) < max_params, ints_to_tack_on)
            ints_to_tack_on = Filter(function(x) length(strsplit(x, ':')[[1]]) <= max_interaction, ints_to_tack_on)

            jj = min(length(ints_to_tack_on), (max_params - 2))
            for(j in seq_len(jj)){

                int_trm_grps = combn(ints_to_tack_on, j, simplify = FALSE)
                int_trm_grps = Filter(function(x) length(x) < (max_params - 2), int_trm_grps)

                for(grp in int_trm_grps){

                    modind = modind + 1
                    mods[[modind]] = as.formula(glue('{y} ~ {xs} + {xs_int}',
                                                xs = paste(indeps[[i]], collapse = '+'),
                                                xs_int = paste(grp, collapse = '+')))

                    if(through_origin){
                        modind = modind + 1
                        mods[[modind]] = as.formula(glue('{y} ~ 0 + {xs} + {xs_int}',
                                                         xs = paste(indeps[[i]], collapse = '+'),
                                                         xs_int = paste(grp, collapse = '+')))
                    }
                }
            }
        }
    }

    param_counts = sapply(mods, function(x) sum(grepl('season', attr(terms(x), 'term.labels')) * 2 + 1)) + 1
    is_through_origin = sapply(mods, function(x) grepl('^0 +', as.character(x)[3]))
    param_counts = param_counts - is_through_origin
    mods = mods[param_counts <= max_params]

    return(mods)
}
eval_model_set = function(data, model_list, metric, log = 'xy', unscale_q_by_area = T){

    #data: a data frame with a date column, and all columns referenced in model_list
    #model_list: a list of lm formulae, as generated by generate_nested_formulae
    #metric: [IGNORED - uses NSE] a function for evaluating the models. the first argument to this function
    #   will be the predictions, and the second will be the observations for the dependent variable

    best_mod = NA
    best_score = -Inf
    # do_not_center = c('date', 'site_code', 'season', 'discharge', 'discharge_log', 'discharge_neon_orig', 'discharge_neon_cont')
    for(i in seq_along(model_list)){

        try({

            trms = rownames(attributes(terms(model_list[[i]]))$factors)
            indeps = gsub('`', '', trms[-1])
            dep = gsub('`', '', trms[1])
            dep_transformed = sub('_log', '', dep)

            dd = filter(data, if_all(all_of(c(dep, indeps)), ~ ! is.na(.)))
            # dd = mutate(dd, across(-any_of(do_not_center), #for centering indeps. not necessary, since multicollinearity does not affect predictions
            #                   ~ . - mean(., na.rm = TRUE)))
            ndata_per_fold = sum(! is.na(dd[[dep]])) %/% 30
            ndata_per_fold = ifelse(ndata_per_fold == 0, 1, ndata_per_fold)

            tryCatch({
                dd = CVlm(data = dd, model_list[[i]], m = 3,
                          plotit = FALSE, printit = FALSE)
            }, warning = function(w) stop('rank-deficiency detected. moving on.'))

            if(log == 'xy'){

                neon_site = data$site_code[1]
                if(unscale_q_by_area){
                    wsa = filter(ms_areas, site_code == neon_site) %>% pull(ws_area_ha)
                    nsenum = sum((inv_neglog(dd$cvpred) * wsa / 1000 - inv_neglog(dd[[dep]]) * wsa / 1000)^2)
                } else {
                    nsenum = sum((inv_neglog(dd$cvpred) - inv_neglog(dd[[dep]]))^2)
                }

            } else if(log == 'y'){
                stop('do not use this')
                nsenum = sum((inv_neglog(dd$cvpred) - dd[[dep]])^2)
            } else if(log == FALSE){
                stop('do not use this')
                nsenum = sum((dd$cvpred - dd[[dep]])^2)
            }

            if(unscale_q_by_area){
                obs = pull(dd[dep_transformed]) * wsa / 1000
            } else {
                obs = pull(dd[dep_transformed])
            }

            nsedenom = sum((obs - mean(obs))^2)
            nse = 1 - (nsenum / nsedenom)

            if(nse > best_score){
                best_score = nse
                best_mod = i
            }
        })
    }

    trms = rownames(attributes(terms(model_list[[best_mod]]))$factors)
    dep = gsub('`', '', trms[1])
    indeps = gsub('`', '', trms[-1])
    if(length(indeps) == 1 && indeps == 'season') stop('season is the only indep. rubbish site(s)')
    site_indeps = grep('season', indeps, invert = TRUE, value = TRUE)

    #for un-centering (not right yet. since the means are based on dd, and data has more obs)
    # dd = filter(data, if_all(all_of(c(dep, indeps)), ~ ! is.na(.)))
    # ddmeans = summarize(dd, across(-any_of(do_not_center),
    #                                ~mean(., na.rm = TRUE))) %>%
    #     rename_with(~paste0(., '_mean'))
    # dd = mutate(dd, across(-any_of(do_not_center),
    #                        ~ . - mean(., na.rm = TRUE)))
    # for(colm in setdiff(colnames(data), c(do_not_center, 'fold', 'Predicted', 'cvpred'))){
    #     data[[colm]] = data[[colm]] - ddmeans[[paste0(colm, '_mean')]]
    # }

    dd = filter(data, if_all(all_of(c(dep, indeps)), ~ ! is.na(.)))
    m = lm(model_list[[best_mod]], data = dd)

    ggps = list()
    yvar = ifelse(log %in% c('y', 'xy'), 'discharge_log', 'discharge')
    for(i in seq_along(site_indeps)){
        d2 = filter(data, ! is.na(!!sym(yvar)) & ! is.na(!!sym(site_indeps[i])))
        ggps[[i]] = ggplot(d2, aes(x = !!sym(site_indeps[i]), y = !!sym(yvar))) +
            geom_point() +
            stat_smooth(method = "lm", col = "red")
    }
    gd = do.call("grid.arrange", c(ggps))
    print(gd)

    newdata = select(data, all_of(indeps))

    if('season' %in% indeps){
        modeled_seasons = m$xlevels$season
        modseas_inds = as.character(newdata$season) %in% modeled_seasons
    } else {
        modeled_seasons = as.character(1:4)
        modseas_inds = rep(TRUE, nrow(newdata))
    }

    data$lm = NA_real_

    if(log %in% c('y', 'xy')){

        if(unscale_q_by_area){
            data$lm[modseas_inds] = inv_neglog(predict(m, newdata = newdata[modseas_inds, ])) * wsa / 1000
            data$discharge = data$discharge * wsa / 1000
        } else {
            data$lm[modseas_inds] = inv_neglog(predict(m, newdata = newdata[modseas_inds, ]))
        }

    } else {
        stop('unused')
        data$lm[modseas_inds] = predict(m, newdata = newdata[modseas_inds, ])
    }

    data$lm[data$lm < 0] = 0

    nse_out = hydroGOF::NSE(data$lm, c(data$discharge))

    first_non_na = Position(function(x) ! is.na(x), data$lm)
    last_non_na = nrow(data) - Position(function(x) ! is.na(x), rev(data$lm)) + 1
    plot_data = data[first_non_na:last_non_na, ]

    site_indeps = c(site_indeps, sub('_log', '', site_indeps))
    drop_cols = grep('^[0-9]', colnames(plot_data), value = TRUE)
    drop_cols = drop_cols[! drop_cols %in% site_indeps]
    plot_data = select(plot_data, -any_of(drop_cols), -ends_with('_log')) %>%
        select(site_code, any_of(c('date', 'datetime')), Q_neon_field = discharge, Q_predicted = lm,
               everything())

    if(unscale_q_by_area){

        nearby_usgs_gages = grep('^[0-9]+$', colnames(plot_data), value = TRUE)
        for(g in nearby_usgs_gages){

            if(g == '06190540'){
                plot_data[[g]] = plot_data[[g]] * 51089.73 #ha
                next
            }

            nwissite = dataRetrieval::readNWISsite(g)

            if(is.null(nwissite)) stop('cannot find nwis watershed area')

            if(is.na(nwissite$contrib_drain_area_va)){
                wsa = nwissite$drain_area_va * 258.999 #mi^2 -> ha
            } else {
                wsa = nwissite$contrib_drain_area_va * 258.999 #mi^2 -> ha
            }

            plot_data[[g]] = plot_data[[g]] * wsa / 1000
        }

        nearby_ms_gages = select(plot_data, -site_code, -any_of(c('date', 'datetime')),
                                 -starts_with('Q_'), -season, -any_of(nearby_usgs_gages)) %>%
            colnames()

        for(g in nearby_ms_gages){
            wsa = filter(ms_areas, site_code == g) %>% pull(ws_area_ha)
            plot_data[[g]] = plot_data[[g]] * wsa / 1000
        }
    }

    # site_code = data$site_code[1]
    #
    # if(file.exists(ncdf)){
    #
    #     # xx = reticulate::py_load_object('../imputation/src/nh_methods/runs/run1360_1907_011320/test/model_epoch030/test_results.p')
    #     xx = reticulate::py_load_object('../imputation/src/nh_methods/runs/run1361_1907_032422/test/model_epoch030/test_results.p')
    #
    #     pred = xx[[paste0(site, '_GAPPED')]]$`1D`$xr$discharge_sim$to_pandas()
    #     pred = tibble(date = as.Date(rownames(pred)), Q = pred$`0`)
    # }

    # plot_data = select(plot_data, date, Q_predicted = lm, Q_used_in_regression = discharge,
    #        Q_neon_continuous_filtered = discharge_neon_cont,
    #        Q_neon_continuous_raw = discharge_neon_orig, Q_neon_manual = discharge_manual_forreals)
    #
    # dg = dygraphs::dygraph(xts(x = select(plot_data, -date), order.by = plot_data$date)) %>%
    #     dyRangeSelector()

    out = list(best_model = model_list[[best_mod]],
               best_model_object = m,
               prediction = unname(data$lm),
               lm_data = plot_data,
               score = nse_out,
               score_crossval = best_score,
               # plot = dg,
               fits = gd)

    return(out)
}
plots_and_results = function(neon_site, best, lm_df, results, return_plot = FALSE, unscale_q_by_area = T){

    if(length(best$prediction) != nrow(lm_df)) stop('oi')

    #load corroborating usgs/ms site data

    sites_nearby = read_csv(glue('../imputation/data/usgs_ms_q/{neon_site}.csv'))

    #filter usgs/ms site data that didn't end up in the model
    trms = rownames(attributes(terms(best$best_model))$factors)
    dep = gsub('`', '', trms[1])
    indeps = gsub('`', '', trms[-1])
    indeps = gsub('_log', '', indeps)
    if(length(indeps) == 1 && indeps == 'season') stop('season is the only indep. rubbish site(s)')
    site_indeps = grep('season', indeps, invert = TRUE, value = TRUE)
    site_indeps_log = paste0(site_indeps, '_log')

    sites_nearby = sites_nearby %>%
        select(datetime, all_of(site_indeps), all_of(site_indeps_log))

    if('season' %in% indeps){
        sites_nearby$season = factor(lubridate::quarter(sites_nearby$datetime))
    }

    #assemble neon sensor data, filtered neon sensor data, neon field data
    #into one frame and plot it

    neon_q_auto = read_csv(glue('../imputation/data/neon_continuous_Q/{neon_site}.csv')) %>%
        filter(! is.na(discharge)) %>%
        rename(discharge_auto = discharge)

    q_eval = read_csv('../neon_discharge_eval/data/neon_q_eval_final.csv') %>%
        filter(site == neon_site)

    check1 = is.na(q_eval$regression_status) | q_eval$regression_status %in% c('good')
    check2 = is.na(q_eval$drift_status) | q_eval$drift_status %in% c('likely_no_drift', 'not_assessed')
    check3 = is.na(q_eval$rating_curve_status) | q_eval$rating_curve_status %in% c('Tier1')

    q_eval$keep = check1 & check2 & check3
    q_eval = q_eval %>%
        group_by(site, year, month) %>%
        summarize(keep = any(keep), #not being strict here. won't be using sensor data in regression anyway
                  .groups = 'drop')

    neon_q_auto_qc = neon_q_auto %>%
        mutate(year = year(datetime),
               month = month(datetime)) %>%
        left_join(select(q_eval, year, month, keep),
                  by = c('year', 'month')) %>%
        filter(keep) %>%
        select(-keep, -year, -month) %>%
        rename(discharge_auto_qc = discharge_auto)

    #predict Q for all datetimes with predictor data

    qall = left_join(sites_nearby, neon_q_auto, by = 'datetime') %>%
        left_join(select(neon_q_auto_qc, -site_code), by = 'datetime')

    if(inherits(best$best_model_object, 'segmented')){
        qall = predict(best$best_model_object,
                       newdata = select(qall, x1 = site_indeps_log),
                       interval = 'predict') %>%
            as_tibble() %>%
            mutate(across(everything(), ~ifelse(. < 0, 0, .))) %>%
            mutate(across(everything(), inv_neglog)) %>%
            bind_cols(qall)
    } else {
        qall = predict(best$best_model_object,
                       newdata = select(qall, all_of(site_indeps_log), any_of('season')),
                       interval = 'predict') %>%
            as_tibble() %>%
            mutate(across(everything(), ~ifelse(. < 0, 0, .))) %>%
            mutate(across(everything(), inv_neglog)) %>%
            bind_cols(qall)
    }

    if(unscale_q_by_area){
        wsa = filter(ms_areas, site_code == neon_site) %>% pull(ws_area_ha)
        qall = mutate(qall, fit = fit * wsa / 1000, lwr = lwr * wsa / 1000, upr = upr * wsa / 1000)
        # qall = mutate(qall, across(c(fit, lwr, upr), . * wsa))
    }

    out_data = qall %>%
        select(-ends_with('_log'), -any_of('season')) %>%
        full_join(select(best$lm_data, datetime, Q_neon_field), by = 'datetime') %>%
        select(datetime, Q_predicted = fit, #Q_used_in_regression = discharge,
               Q_pred_int_2.5 = lwr, Q_pred_int_97.5 = upr, Q_neon_field,
               Q_neon_continuous_filtered = discharge_auto_qc,
               Q_neon_continuous_raw = discharge_auto) %>%
        filter(if_any(-datetime, ~cumsum(! is.na(.)) != 0)) %>%  #remove leading NA rows
        arrange(datetime)

    dg = dygraphs::dygraph(xts(x = select(out_data, -datetime, -Q_pred_int_2.5, -Q_pred_int_97.5) %>% tail(5e5),
                               order.by = tail(out_data$datetime, 5e5))) %>%
        dyRangeSelector()

    saveWidget(dg, glue('../imputation/out/q_lm_plots/pred/{neon_site}_log.html'))

    if(inherits(best$fits, 'grob')){

        #make diagnostic and fit plots
        png(glue('../imputation/out/q_lm_plots/diag/{neon_site}_diag_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
        defpar = par(mfrow=c(2,2))
        plot(best$best_model_object)
        par(defpar)
        dev.off()

        png(glue('../imputation/out/q_lm_plots/fit/{neon_site}_fit_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
        plot(best$fits)
        dev.off()
    }

    #plot predictions versus field measurements. need to round field meas to Q interval
    time_diffs = diff(sites_nearby$datetime)
    units(time_diffs) = 'mins'
    time_interval = Mode(time_diffs)

    zz = out_data %>%
        mutate(datetime = round_date(datetime, paste(time_interval, 'min')))

    field_dts = filter(zz, ! is.na(Q_neon_field)) %>% pull(datetime)

    zz = zz %>%
        filter(datetime %in% field_dts) %>%
        group_by(datetime) %>%
        summarize(across(everything(), ~mean(., na.rm = TRUE))) %>%
        ungroup()

    axlim = c(0, max(c(zz$Q_predicted, zz$Q_neon_field), na.rm = TRUE))

    png(glue('../imputation/out/q_lm_plots/val/{neon_site}_obs_v_pred.png'), 6, 6, 'in', type = 'cairo', res = 300)
    plot(zz$Q_neon_field, zz$Q_predicted, xlab = 'NEON Field Discharge (L/s)',
         ylab = 'Predicted Discharge (L/s)', main = glue('Site: {neon_site}; NSE: {nse1}; NSE crossval: {nse2}',
                                                         nse1 = round(best$score, 2),
                                                         nse2 = round(best$score_crossval, 2)),
         xlim = axlim, ylim = axlim, xaxs = 'i', yaxs = 'i')
    abline(a = 0, b = 1, col = 'blue')
    legend('topleft', legend = '1:1', lty = 1, col = 'blue', bty = 'n')
    dev.off()

    out_data = filter(out_data, is.na(Q_neon_field)) %>% select(-Q_neon_field)

    #save predictions as CSV
    out_data = left_join(out_data,
                         select(sites_nearby, datetime, all_of(site_indeps), any_of('season')),
                         by = 'datetime')
    write_csv(out_data, glue('../imputation/out/q_lm_outdata/predictions/{neon_site}.csv'))

    #save fit data as CSV
    write_csv(best$lm_data, glue('../imputation/out/q_lm_outdata/fit/{neon_site}.csv'))

    #save model summary
    sink(glue('../imputation/out/q_lm_outdata/summary/{neon_site}.txt'))
    print(summary(best$best_model_object))
    sink()

    #return results frame, updated with this site's results
    results$bestmod_logq[results$site_code == neon_site] = as.character(best$best_model)[3]
    results$nse_logq[results$site_code == neon_site] = best$score
    results$nse_cv_logq[results$site_code == neon_site] = best$score_crossval
    results$adj_r_squared[results$site_code == neon_site] = summary(best$best_model_object)$adj.r.squared

    if(return_plot){
        return(list(plot = dg, results = results))
    } else {
        return(results)
    }
}
plots_and_results_daily_composite = function(neon_site, best1, best2, lm_df1, lm_df2, results, unscale_q_by_area = T){

    #best2 and lm_df2 should represent the model with more terms included

    if(length(best1$prediction) != nrow(lm_df1)) stop('oi')
    if(length(best2$prediction) != nrow(lm_df2)) stop('oi')

    #load corroborating usgs/ms site data

    sites_nearby = read_csv(glue('../imputation/data/usgs_ms_q/{neon_site}.csv'))

    #filter usgs/ms site data that didn't end up in the model
    trms1 = rownames(attributes(terms(best1$best_model))$factors)
    dep1 = gsub('`', '', trms1[1])
    indeps1 = gsub('`', '', trms1[-1])
    indeps1 = gsub('_log', '', indeps1)
    if(length(indeps1) == 1 && indeps1 == 'season') stop('season is the only indep. rubbish site(s)')
    site_indeps1 = grep('season', indeps1, invert = TRUE, value = TRUE)
    site_indeps_log1 = paste0(site_indeps1, '_log')

    trms2 = rownames(attributes(terms(best2$best_model))$factors)
    dep2 = gsub('`', '', trms2[1])
    indeps2 = gsub('`', '', trms2[-1])
    indeps2 = gsub('_log', '', indeps2)
    if(length(indeps2) == 1 && indeps2 == 'season') stop('season is the only indep. rubbish site(s)')
    site_indeps2 = grep('season', indeps2, invert = TRUE, value = TRUE)
    site_indeps_log2 = paste0(site_indeps2, '_log')

    if(length(site_indeps2) <= length(site_indeps1)) stop('model 2 should have more sites than 1')

    sites_nearby = sites_nearby %>%
        select(date, all_of(site_indeps2), all_of(site_indeps_log2))

    if('season' %in% indeps1 | 'season' %in% indeps2){
        sites_nearby$season = factor(lubridate::quarter(sites_nearby$date))
    }

    #assemble neon sensor data, filtered neon sensor data, neon field data
    #into one frame and plot it

    neon_q_daily = read_csv(glue('../imputation/data/neon_continuous_Q/{neon_site}.csv')) %>%
        group_by(date = date(datetime)) %>%
        summarize(discharge = mean(discharge, na.rm = TRUE)) %>%
        filter(! is.na(discharge)) %>%
        rename(discharge_daily = discharge)

    q_eval = read_csv('../neon_discharge_eval/data/neon_q_eval_final.csv') %>%
        filter(site == neon_site)

    check1 = is.na(q_eval$regression_status) | q_eval$regression_status %in% c('good')
    check2 = is.na(q_eval$drift_status) | q_eval$drift_status %in% c('likely_no_drift', 'not_assessed')
    check3 = is.na(q_eval$rating_curve_status) | q_eval$rating_curve_status %in% c('Tier1')

    q_eval$keep = check1 & check2 & check3
    q_eval = q_eval %>%
        group_by(site, year, month) %>%
        summarize(keep = any(keep), #not being strict here. won't be using sensor data in regression anyway
                  .groups = 'drop')

    neon_q_daily_qc = neon_q_daily %>%
        mutate(year = year(date),
               month = month(date)) %>%
        left_join(select(q_eval, year, month, keep),
                  by = c('year', 'month')) %>%
        filter(keep) %>%
        select(-keep, -year, -month) %>%
        rename(discharge_daily_qc = discharge_daily)

    #predict Q for all datetimes with predictor data

    qall = left_join(sites_nearby, neon_q_daily, by = 'date') %>%
        left_join(neon_q_daily_qc, by = 'date')

    # qall$pred1 = predict(best1$best_model_object,
    #                      newdata = select(qall, all_of(site_indeps_log1), any_of('season'))) %>%
    #     inv_neglog()
    # qall$pred = predict(best2$best_model_object,
    #                      newdata = select(qall, all_of(site_indeps_log2), any_of('season'))) %>%
    #     inv_neglog()
    # missing_preds = is.na(qall$pred)
    # qall$pred[missing_preds] = qall$pred1[missing_preds]
    # qall$pred1 = NULL

    pred1 = predict(best1$best_model_object,
                    newdata = select(qall, all_of(site_indeps_log1), any_of('season')),
                    interval = 'predict') %>%
        as_tibble() %>%
        mutate(across(everything(), ~ifelse(. < 0, 0, .))) %>%
        mutate(across(everything(), inv_neglog))

    pred = predict(best2$best_model_object,
                   newdata = select(qall, all_of(site_indeps_log2), any_of('season')),
                   interval = 'predict') %>%
        as_tibble() %>%
        mutate(across(everything(), ~ifelse(. < 0, 0, .))) %>%
        mutate(across(everything(), inv_neglog))

    missing_preds = is.na(pred$fit)
    pred[missing_preds, ] = pred1[missing_preds, ]
    qall = bind_cols(qall, pred)

    if(unscale_q_by_area){
        wsa = filter(ms_areas, site_code == neon_site) %>% pull(ws_area_ha)
        qall = mutate(qall, fit = fit * wsa / 1000, lwr = lwr * wsa / 1000, upr = upr * wsa / 1000)
    }

    out_data = qall %>%
        select(-ends_with('_log'), -any_of('season')) %>%
        full_join(select(best2$lm_data, date, Q_neon_field), by = 'date') %>%
        select(date, Q_predicted = fit,
               Q_pred_int_2.5 = lwr, Q_pred_int_97.5 = upr, Q_neon_field,
               Q_neon_continuous_filtered = discharge_daily_qc,
               Q_neon_continuous_raw = discharge_daily) %>%
        filter(if_any(-date, ~cumsum(! is.na(.)) != 0)) %>%  #remove leading NA rows
        arrange(date)

    dg = dygraphs::dygraph(xts(x = select(out_data, -date, Q_pred_int_2.5, -Q_pred_int_97.5) %>% tail(5e5),
                               order.by = tail(out_data$date, 5e5))) %>%
        dyRangeSelector()

    saveWidget(dg, glue('../imputation/out/q_lm_plots/pred/{neon_site}_log.html'))

    #make diagnostic and fit plots
    png(glue('../imputation/out/q_lm_plots/diag/{neon_site}_diag_log_modA.png'), 6, 6, 'in', type = 'cairo', res = 300)
    defpar = par(mfrow=c(2,2))
    plot(best1$best_model_object)
    par(defpar)
    dev.off()

    png(glue('../imputation/out/q_lm_plots/diag/{neon_site}_diag_log_modB.png'), 6, 6, 'in', type = 'cairo', res = 300)
    defpar = par(mfrow=c(2,2))
    plot(best2$best_model_object)
    par(defpar)
    dev.off()

    png(glue('../imputation/out/q_lm_plots/fit/{neon_site}_fit_log_modA.png'), 6, 6, 'in', type = 'cairo', res = 300)
    plot(best1$fits)
    dev.off()

    png(glue('../imputation/out/q_lm_plots/fit/{neon_site}_fit_log_modB.png'), 6, 6, 'in', type = 'cairo', res = 300)
    plot(best2$fits)
    dev.off()

    #plot predictions versus field measurements. need to round field meas to Q interval
    nse = hydroGOF::NSE(out_data$Q_predicted, out_data$Q_neon_field)
    axlim = c(0, max(c(out_data$Q_predicted, out_data$Q_neon_field), na.rm = TRUE))

    png(glue('../imputation/out/q_lm_plots/val/{neon_site}_obs_v_pred.png'), 6, 6, 'in', type = 'cairo', res = 300)
    plot(out_data$Q_neon_field, out_data$Q_predicted, xlab = 'NEON Field Discharge (L/s)',
         ylab = 'Predicted Discharge (L/s)', main = glue('Site: {neon_site}; NSE: {nse1}; NSE crossval: {nse2}',
                                                         nse1 = round(nse, 2),
                                                         nse2 = NA),
         xlim = axlim, ylim = axlim, xaxs = 'i', yaxs = 'i')
    abline(a = 0, b = 1, col = 'blue')
    legend('topleft', legend = '1:1', lty = 1, col = 'blue', bty = 'n')
    dev.off()

    #save predictions as CSV
    out_data = filter(out_data, is.na(Q_neon_field)) %>% select(-Q_neon_field)
    out_data = left_join(out_data,
                         select(sites_nearby, date, all_of(site_indeps2), any_of('season')),
                         by = 'date')
    write_csv(out_data, glue('../imputation/out/q_lm_outdata/predictions/{neon_site}.csv'))

    #save fit data as CSV
    write_csv(best2$lm_data, glue('../imputation/out/q_lm_outdata/fit/{neon_site}.csv'))

    #save model summary
    sink(glue('../imputation/out/q_lm_outdata/summary/{neon_site}.txt'))
    print(summary(best1$best_model_object))
    print('---')
    print(summary(best2$best_model_object))
    sink()

    #return results frame, updated with this site's results
    results$bestmod_logq[results$site_code == neon_site] = paste(as.character(best1$best_model)[3],
                                                                 as.character(best1$best_model)[3],
                                                                 sep = ' & ')
    results$nse_logq[results$site_code == neon_site] = nse
    results$nse_cv_logq[results$site_code == neon_site] = NA_real_
    results$adj_r_squared[results$site_code == neon_site] = mean(summary(best1$best_model_object)$adj.r.squared,
                                                                 summary(best2$best_model_object)$adj.r.squared)

    return(results)
}

#obsolete helpers
plots_and_results_daily = function(neon_site, best, lm_df, results){

    stop('not updated')
    if(length(best$prediction) != nrow(lm_df)) stop('oi')

    d = lm_df %>%
        select(-ends_with('_log'), -season, -matches('[0-9]+'))

    #load corroborating usgs/ms site data

    sites_nearby = read_csv(glue('../imputation/data/usgs_ms_q/{neon_site}.csv'))

    #filter usgs/ms site data that didn't end up in the model
    trms = rownames(attributes(terms(best$best_model))$factors)
    dep = gsub('`', '', trms[1])
    indeps = gsub('`', '', trms[-1])
    indeps = gsub('_log', '', indeps)
    if(length(indeps) == 1 && indeps == 'season') stop('season is the only indep. rubbish site(s)')
    site_indeps = grep('season', indeps, invert = TRUE, value = TRUE)
    site_indeps_log = paste0(site_indeps, '_log')

    sites_nearby = sites_nearby %>%
        select(date, all_of(site_indeps), all_of(site_indeps_log))

    if('season' %in% indeps){
        sites_nearby$season = factor(lubridate::quarter(sites_nearby$date))
    }

    #assemble neon sensor data, filtered neon sensor data, neon field data
    #into one frame and plot it

    neon_q_daily = read_csv(glue('../imputation/data/neon_continuous_Q/{neon_site}.csv')) %>%
        group_by(date = date(datetime)) %>%
        summarize(discharge = mean(discharge, na.rm = TRUE)) %>%
        filter(! is.na(discharge)) %>%
        rename(discharge_daily = discharge)

    q_eval = read_csv('../neon_discharge_eval/data/neon_q_eval_final.csv') %>%
        filter(site == neon_site)

    check1 = is.na(q_eval$regression_status) | q_eval$regression_status %in% c('good')
    check2 = is.na(q_eval$drift_status) | q_eval$drift_status %in% c('likely_no_drift', 'not_assessed')
    check3 = is.na(q_eval$rating_curve_status) | q_eval$rating_curve_status %in% c('Tier1')

    q_eval$keep = check1 & check2 & check3
    q_eval = q_eval %>%
        group_by(site, year, month) %>%
        summarize(keep = any(keep), #not being strict here. won't be using sensor data in regression anyway
                  .groups = 'drop')

    neon_q_daily_qc = neon_q_daily %>%
        mutate(year = year(date),
               month = month(date)) %>%
        left_join(select(q_eval, year, month, keep),
                  by = c('year', 'month')) %>%
        filter(keep) %>%
        select(-keep, -year, -month) %>%
        rename(discharge_daily_qc = discharge_daily)

    #predict Q for all datetimes with predictor data


    qall = left_join(sites_nearby, neon_q_daily, by = 'date') %>%
        left_join(neon_q_daily_qc, by = 'date')

    if(inherits(best$best_model_object, 'segmented')){
        qall$pred = predict(best$best_model_object,
                            newdata = select(qall, x1 = site_indeps_log)) %>%
            inv_neglog()
    } else {
        qall$pred = predict(best$best_model_object,
                            newdata = select(qall, all_of(site_indeps_log), any_of('season'))) %>%
            inv_neglog()
    }

    out_data = qall %>%
        select(-ends_with('_log'), -any_of('season')) %>%
        full_join(select(lm_df, date, discharge), by = 'date') %>%
        select(date, Q_predicted = pred, #Q_used_in_regression = discharge,
           Q_neon_field = discharge,
           Q_neon_daily_filtered = discharge_daily_qc,
           Q_neon_daily_raw = discharge_daily) %>%
        filter(if_any(-date, ~cumsum(! is.na(.)) != 0)) %>%  #remove leading NA rows
        arrange(date)

    dg = dygraphs::dygraph(xts(x = select(out_data, -date) %>% tail(5e5),
                               order.by = tail(out_data$date, 5e5))) %>%
        dyRangeSelector()

    saveWidget(dg, glue('../imputation/out/q_lm_plots/pred/{neon_site}_log.html'))

    if(inherits(best$fits, 'grob')){

        #make diagnostic and fit plots
        png(glue('../imputation/out/q_lm_plots/diag/{neon_site}_diag_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
        defpar = par(mfrow=c(2,2))
        plot(best$best_model_object)
        par(defpar)
        dev.off()

        png(glue('../imputation/out/q_lm_plots/fit/{neon_site}_fit_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
        plot(best$fits)
        dev.off()
    }

    # #plot predictions versus field measurements. need to round field meas to Q interval
    # time_diffs = diff(sites_nearby$date)
    # units(time_diffs) = 'mins'
    # time_interval = Mode(time_diffs)
    #
    # zz = out_data %>%
    #     mutate(date = round_date(date, paste(time_interval, 'min')))
    #
    # field_dts = filter(zz, ! is.na(Q_neon_field)) %>% pull(date)
    #
    # zz = zz %>%
    #     filter(date %in% field_dts) %>%
    #     group_by(date) %>%
    #     summarize(across(everything(), ~mean(., na.rm = TRUE))) %>%
    #     ungroup()

    axlim = c(0, max(c(out_data$Q_predicted, out_data$Q_neon_field), na.rm = TRUE))

    png(glue('../imputation/out/q_lm_plots/val/{neon_site}_obs_v_pred.png'), 6, 6, 'in', type = 'cairo', res = 300)
    plot(out_data$Q_neon_field, out_data$Q_predicted, xlab = 'NEON Field Discharge (L/s)',
         ylab = 'Predicted Discharge (L/s)', main = glue('Site: {neon_site}; NSE: {nse1}; NSE crossval: {nse2}',
                                                         nse1 = round(best$score, 2),
                                                         nse2 = round(best$score_crossval, 2)),
         xlim = axlim, ylim = axlim, xaxs = 'i', yaxs = 'i')
    abline(a = 0, b = 1, col = 'blue')
    legend('topleft', legend = '1:1', lty = 1, col = 'blue', bty = 'n')
    dev.off()

    #save combined frame
    out_data = left_join(out_data,
                         select(sites_nearby, date, all_of(site_indeps), any_of('season')),
                         by = 'date')
    write_csv(out_data, glue('../imputation/out/q_lm_outdata/by_site/{neon_site}.csv'))

    #return results frame, updated with this site's results
    results$bestmod_logq[results$site_code == neon_site] = as.character(best$best_model)[3]
    results$nse_logq[results$site_code == neon_site] = best$score
    results$nse_cv_logq[results$site_code == neon_site] = best$score_crossval

    return(results)
}
assemble_output_df = function(neon_site, nearby_usgs_gages = NULL, nearby_ms_gages = NULL,
                              include_pct_highvals = NULL, include_sensor_data = FALSE, include_sensor_daterange = NULL){

    neon_q_manual = read_csv(glue('../imputation/data/neon_field_Q/{neon_site}.csv')) %>%
        mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
        rename(discharge_manual = discharge)
    earliest_date = as.character(date(min(neon_q_manual$datetime)))

    # neon_q_auto = read_csv(glue('../imputation/data/neon_continuous_Q/{neon_site}.csv')) %>%
    #     filter(! is.na(discharge)) %>%
    #     rename(discharge_auto = discharge)

    site_nearbyA = NULL
    if(! is.null(nearby_usgs_gages)){

        usgsq = dataRetrieval::readNWISdata(
            sites = nearby_usgs_gages,
            service = 'iv', #instantaneous values
            parameterCd = '00060', #discharge (cfs)
            startDate = earliest_date,
            endDate = '2099-01-01'
        )

        if(any(! usgsq$X_00060_00000_cd %in% c('A', 'P', 'A e'))) stop()
        if(any(usgsq$tz_cd != 'UTC')) stop()

        site_nearbyA = usgsq %>%
            as_tibble() %>%
            select(site_no, datetime = dateTime, discharge = X_00060_00000) %>%
            pivot_wider(names_from = site_no, values_from = discharge) %>%
            # mutate(site_code = neon_site) %>%
            mutate(across(matches('^[0-9]+$'), ~ . * 28.3168)) %>%
                   # date = as.Date(date)) %>%
            mutate(across(matches('^[0-9]+$'), neglog, .names = '{.col}_log')) %>%
            arrange(datetime)
    }

    site_nearbyB = NULL
    if(! is.null(nearby_ms_gages)){

        stop('use inst q for these sites if available')

        site_nearbyB = Filter(function(x) x$site_code[1] %in% nearby_ms_gages, q_data) %>%
            reduce(bind_rows) %>%
            pivot_wider(names_from = site_code, values_from = discharge) %>%
            mutate(across(-date, neglog, .names = '{.col}_log')) %>%
            arrange(date) %>%
            mutate(site_code = neon_site)
    }

    if(! is.null(site_nearbyA) && ! is.null(site_nearbyB)){
        site_nearby = full_join(site_nearbyA, site_nearbyB, by = c('site_code', 'date'))
    } else if(! is.null(site_nearbyA)){
        site_nearby = site_nearbyA
    } else  {
        site_nearby = site_nearbyB
    }

    # q_eval = read_csv('../neon_discharge_eval/data/neon_q_eval_final.csv') %>%
    #     filter(site == neon_site)
    #
    # check1 = is.na(q_eval$regression_status) | q_eval$regression_status %in% c('good')
    # check2 = is.na(q_eval$drift_status) | q_eval$drift_status %in% c('likely_no_drift', 'not_assessed')
    # check3 = is.na(q_eval$rating_curve_status) | q_eval$rating_curve_status %in% c('Tier1')
    #
    # q_eval$keep = check1 & check2 & check3
    # q_eval = q_eval %>%
    #     group_by(site, year, month) %>%
    #     summarize(keep = any(keep), #not being strict here. won't be using sensor data in regression anyway
    #               .groups = 'drop')
    #
    # neon_q_auto_qc = neon_q_auto %>%
    #     mutate(year = year(datetime),
    #            month = month(datetime)) %>%
    #     left_join(select(q_eval, year, month, keep),
    #               by = c('year', 'month')) %>%
    #     filter(keep) %>%
    #     select(-keep, -year, -month) %>%
    #     rename(discharge_auto_qc = discharge_auto)
    #
    # # if(neon_site == 'CARI') neon_q_auto_qc$discharge[neon_q_auto_qc$date > as.Date('2020-02-23') & neon_q_auto_qc$date < as.Date('2020-03-31')] = NA
    # # if(neon_site == 'CUPE') neon_q_auto_qc = mutate(neon_q_auto_qc, discharge = ifelse(date >= as.Date('2020-01-01'), discharge, NA))
    # # if(neon_site == 'KING') neon_q_auto_qc$discharge[neon_q_auto_qc$date < as.Date('2018-12-01') | neon_q_auto_qc$date > as.Date('2019-05-06')] = NA
    # # if(neon_site == 'LEWI') neon_q_auto_qc$discharge[neon_q_auto_qc$date > as.Date('2020-03-05') & neon_q_auto_qc$date < as.Date('2020-07-25')] = NA
    # # if(neon_site == 'COMO') neon_q_auto_qc$discharge[neon_q_auto_qc$date > as.Date('2021-01-16') & neon_q_auto_qc$date < as.Date('2021-03-06')] = NA
    # # if(neon_site == 'MCDI') neon_q_auto_qc = mutate(neon_q_auto_qc, discharge = ifelse(date <= as.Date('2020-07-13'), discharge, NA))
    # # if(neon_site == 'MCRA') neon_q_auto_qc = mutate(neon_q_auto_qc, discharge = ifelse(date <= as.Date('2019-05-06'), discharge, NA))
    # # if(neon_site == 'BLDE') neon_q_auto_qc = mutate(neon_q_auto_qc, discharge = ifelse(date <= as.Date('2020-10-01'), discharge, NA))


    # neon_q_orig = select(original_neon_q[[which(names(original_neon_q) == neon_site)]],
    #                      date, discharge_neon_orig = discharge)

    # qall = left_join(site_nearby, neon_q_auto) %>%
    #     left_join(select(neon_q_auto_qc, -site_code))

    x = rename(neon_q_manual, datetime_x = datetime) %>% as.data.table()
    y = rename(site_nearby, datetime_y = datetime) %>% as.data.table()

    rollmax = 60 * 60 * 12 #join up to 12 hrs off
    x[, `:=` (datetime_min = datetime_x - rollmax,
              datetime_max = datetime_x + rollmax)]
    y[, `:=` (datetime_y_orig = datetime_y)] #this datetime col will be dropped

    #join x rows to y if y's datetime falls within the x range
    joined = y[x, on = .(datetime_y <= datetime_max,
                         datetime_y >= datetime_min)]
    joined = na.omit(joined, cols = 'datetime_y_orig') #drop rows without matches

    #for any datetimes in x or y that were matched more than once, keep only
    #the nearest match
    joined[, `:=` (datetime_match_diff = abs(datetime_x - datetime_y_orig))]
    joined = joined[, .SD[which.min(datetime_match_diff)], by = datetime_x]
    joined = joined[, .SD[which.min(datetime_match_diff)], by = datetime_y_orig]

    joined[, c('datetime_y', 'datetime_y.1', 'datetime_y_orig', 'datetime_match_diff') := NULL]
    setnames(joined, 'datetime_x', 'datetime')

    # stop('need to load inst neon q for this')

    # if(! is.null(include_pct_highvals)){
    #
    #     prop = include_pct_highvals / 100
    #
    #     if(! all(is.na(qq_auto))){
    #
    #         #get the largest p*n sensor values, where n is the number of manual measurements and p is a proportion.
    #         #filter any timepoints that are already represented by manual measurements
    #         pn = round(sum(! is.na(qq_manual)) * prop, 0)
    #         qtop = order(qq_auto, decreasing = TRUE)[1:pn]
    #         qtop_df = tibble(qtop = qtop, ddtop = dd[qtop], qqauto_top = qq_auto[qtop])
    #         dd_manual = dd[! is.na(qq_manual)]
    #         qtop_df = filter(qtop_df, ddtop %in% setdiff(ddtop, dd_manual)) %>%
    #             transmute(date = as.Date(as.character(as.Date(ddtop, origin = '1970-01-01'))),
    #                       discharge = qqauto_top)
    #
    #         ddqq = tibble(date = as.Date(as.character(as.Date(dd, origin = '1970-01-01'))),
    #                       discharge = qq_manual) %>%
    #             filter(! date %in% qtop_df$date)
    #
    #         #join those pn values with the manual series and then join it with the rest
    #         # ddqq = bind_rows(ddqq, qtop_df) %>% arrange(date) %>% mutate(site_code = neon_site)
    #         ddqq = full_join(ddqq, qtop_df, by = 'date') %>%
    #             mutate(discharge = case_when(! is.na(discharge.x) ~ discharge.x,
    #                                          ! is.na(discharge.y) ~ discharge.y,
    #                                          TRUE ~ NA_real_)) %>%
    #             arrange(date) %>% mutate(site_code = neon_site) %>%
    #             select(-discharge.x, -discharge.y)
    #
    #         if(nrow(ddqq) != length(qq_auto)) stop('oi')
    #
    #         ddqq$discharge_neon_cont = qq_auto
    #     } else {
    #        message('no sensor values')
    #        ddqq = tibble(date = as.Date(as.character(as.Date(dd, origin = '1970-01-01'))),
    #                      discharge_neon_cont = qq_auto,
    #                      discharge = qq_manual) %>%
    #            mutate(site_code = neon_site)
    #    }
    #
    #} else {
    #    ddqq = tibble(date = as.Date(as.character(as.Date(dd, origin = '1970-01-01'))),
    #                  discharge_neon_cont = qq_auto,
    #                  discharge = qq_manual) %>%
    #        mutate(site_code = neon_site)
    #}

    # nearby_join = full_join(ddqq, site_nearby, by = c('date', 'site_code')) %>%
    #     arrange(date) %>%
    #     full_join(neon_q_orig) %>%
    #     full_join(tibble(date = as.Date(as.character(as.Date(dd[! is.na(qq_manual)], origin = '1970-01-01'))),
    #                      discharge_manual_forreals = qq_manual[! is.na(qq_manual)])) %>%
    #     mutate(site_code = neon_site)

    joined %>%
        as_tibble() %>%
        arrange(datetime) %>%
        rename(discharge = discharge_manual) %>%
        mutate(season = factor(lubridate::quarter(datetime)),
               discharge_log = neglog(discharge)) %>%
        select(site_code, datetime, discharge, discharge_log, everything())

    # nearby_join$season = factor(lubridate::quarter(nearby_join$date))
    # nearby_join$discharge_log = neglog(nearby_join$discharge)

    # nearby_join = arrange(nearby_join, date)
    # nearby_join = lapply(nearby_join, c) %>% as_tibble()

    return(joined)
}
generate_nested_formulaeB = function(full_spec, interaction_term = TRUE, through_origin = FALSE){

    #this one only allows interactions with the specified interaction_term
    #(go back to ce4de33 to get the code)

    return(mods)
}
clean_up = function(lm_df, best){

    trms = rownames(attributes(terms(best$best_model))$factors)
    sites_used = gsub('`', '', trms[-1])
    sites_used = gsub('_log', '', sites_used)
    sites_used = sites_used[! sites_used %in% c('0', 'season')]

    first_non_na = Position(function(x) ! is.na(x), lm_df$best_prediction)
    last_non_na = nrow(lm_df) - Position(function(x) ! is.na(x), rev(lm_df$best_prediction)) + 1

    lm_df = lm_df[first_non_na:last_non_na, ] %>%
        select(site_code, date,
               discharge_predicted = best_prediction,
               discharge_neon_manual = discharge_manual_forreals,
               discharge_neon_continuous_filtered = discharge_neon_cont,
               discharge_neon_continuous_raw = discharge_neon_orig,
               contains(sites_used))

    return(lm_df)
}
eval_model_lmer = function(data, formula, metric, saveloc){

    m = lmer(formula, data = lm_df)

    trms = rownames(attributes(terms(formula))$factors)
    dep = trms[1]
    indeps = gsub('`', '', trms[-1])
    indeps = gsub('1 \\| ', '', indeps)
    site_indeps = grep('season|flow', indeps, invert = TRUE, value = TRUE)

    # pred = predict(m, newdata = select(lm_df, `09510200_log`, flow_level))
    pred = predict(m, newdata = select(data, all_of(indeps)))

    png(saveloc, 6, 6, 'in', type = 'cairo', res = 300)
    plot(data[[site_indeps]], data$discharge_log, xlab = site_indeps, ylab = 'discharge_log')
    o = order(data[[site_indeps]])
    lines(data[[site_indeps]][o], pred[o])
    dev.off()

    data$lm = inv_neglog(pred)
    data$lm[data$lm < 0] = 0

    nse_out = hydroGOF::NSE(data$lm, c(data$discharge))

    first_non_na = Position(function(x) ! is.na(x), data$lm)
    last_non_na = nrow(data) - Position(function(x) ! is.na(x), rev(data$lm)) + 1
    plot_data = data[first_non_na:last_non_na, ]

    dg = dygraphs::dygraph(xts(x = select(plot_data, lm, starts_with(!!dep_transformed), -ends_with('log')),
                               order.by = plot_data$date)) %>%
        dyRangeSelector()

    out = list(formula = formula,
               score = nse_out,
               prediction = unname(data$lm),
               plot = dg)

    return(out)
}

### 2a prepare static attributes of MS basins ####

ms_all <- tibble(site_code = character())
for(d in ms_domains_included){

    ms_clim <- read_feather(paste0('ms_in_camels_format/', d, '/clim.feather'))
    ms_geol <- read_feather(paste0('ms_in_camels_format/', d, '/geol.feather'))
    ms_soil <- try(read_feather(paste0('ms_in_camels_format/', d, '/soil.feather')),
                   silent = TRUE)
    ms_topo <- read_feather(paste0('ms_in_camels_format/', d, '/topo.feather'))
    ms_vege <- try(read_feather(paste0('ms_in_camels_format/', d, '/vege.feather')),
                   silent = TRUE)

    ms_all0 <- ms_geol %>%
        full_join(ms_topo, by = 'site_code') %>%
        # full_join(ms_vege, by = 'site_code') %>%
        full_join(ms_clim, by = 'site_code')

    if(exists('ms_vege') && ! inherits(ms_vege, 'try-error')){
        ms_all0 <- full_join(ms_all0, ms_vege, by = 'site_code')
    }
    if(exists('ms_soil') && ! inherits(ms_soil, 'try-error')){
        ms_all0 <- full_join(ms_all0, ms_soil, by = 'site_code')
    }

    ms_all <- bind_rows(ms_all, ms_all0)

    try(rm(ms_vege), silent = TRUE)
    try(rm(ms_soil), silent = TRUE)
}

ms_all <- filter(ms_all, ! site_code %in% arctic_sites_to_drop)

# ms_all <- ms_all %>%
#     rename(carbonate_rocks_frac = carb_rocks_frac,
#            frac_snow = frac_snow_daily,
#            frac_forest = forest_frac) %>%

### 2b prepare static attributes of CAMELS basins ####

camels_clim <- read_delim('CAMELS/camels_attributes_v2.0/camels_clim.txt',
                          col_types = 'cnnnnnnncnnc', delim = ';')
camels_soil <- read_delim('CAMELS/camels_attributes_v2.0/camels_soil.txt',
                          col_types = 'cnnnnnnnnnnn', delim = ';')
camels_geol <- read_delim('CAMELS/camels_attributes_v2.0/camels_geol.txt',
                          col_types = 'ccncnnnn', delim = ';')
camels_topo <- read_delim('CAMELS/camels_attributes_v2.0/camels_topo.txt',
                          col_types = 'cnnnnnn', delim = ';')
camels_vege <- read_delim('CAMELS/camels_attributes_v2.0/camels_vege.txt',
                          col_types = 'cnnnnnncnn', delim = ';')
# camels_hydro <- read_csv2('CAMELS/camels_attributes_v2.0/camels_hydro.txt',
                         # col_types = cols(.default = 'c'))

camels_attrs <- camels_geol %>%
    full_join(camels_soil, by = 'gauge_id') %>%
    full_join(camels_topo, by = 'gauge_id') %>%
    full_join(camels_vege, by = 'gauge_id') %>%
    full_join(camels_clim, by = 'gauge_id') %>%
    rename(site_code = gauge_id,
           area = area_gages2)

## we were unable to replicate CAMELS methods when generating some attributes for
##  macrosheds sites. for consistency, we replace those original CAMELS attrs
##  with ours.
camels_clim_mod <- read_feather('ms_in_camels_format/camels/clim.feather')
    # filter(! duplicated(.))
camels_soil_mod <- read_feather('ms_in_camels_format/camels/soil.feather')
    # select(-pctCellErr) %>%
    # pivot_wider(names_from = var,
    #             values_from = val)
camels_mod_attrs <- full_join(camels_clim_mod, camels_soil_mod, by = 'site_code')

matched_attrs_bool <- colnames(camels_mod_attrs) %in% colnames(camels_attrs)
# colnames(camels_attrs)
# colnames(camels_mod_attrs)
unmatched_attrs <- colnames(camels_mod_attrs)[! matched_attrs_bool]

camels_attrs <- camels_mod_attrs %>%
    select(-any_of(unmatched_attrs)) %>%
    full_join(camels_attrs,
              by = 'site_code',
              suffix = c('_orig', '_mod')) %>%
    select(-ends_with('_orig')) %>%
    rename_with(function(x) sub('_mod$', '', x))
    # mutate(across(where(is.character),
    #               ~paste0('"', site_code, '"')))

### 2c for sites with NHM data, copy attribute rows and prepend "NHM_" to the site_code ####

# nhm_segids <- filter(nhm_segids, domain == 'neon')
nhm_ms_sites <- ms_all %>%
    filter(site_code %in% nhm_segids$site_code) %>%
    mutate(site_code = paste0('NHM_', site_code))

ms_all <- bind_rows(ms_all, nhm_ms_sites)

nhm_cmls_sites <- camels_attrs %>%
    mutate(site_code = paste0('NHM_', site_code))

camels_attrs <- bind_rows(camels_attrs, nhm_cmls_sites)

### 2d for extrapolation sites, copy attribute rows and append "_extrapolate" to the site_code ####

neon_sites <- read_csv('informative_stuff/neon_basin_areas.csv')$site_code

ms_extrap_sites <- ms_all %>%
    filter(site_code %in% neon_sites) %>%
    mutate(site_code = paste0(site_code, '_extrapolate'))

ms_all <- bind_rows(ms_all, ms_extrap_sites)

### 3a prepare macrosheds timeseries data ####

ms_clim_ts <- tibble()
for(d in ms_domains_included){

    ms_clim_ts0 <- read_feather(paste0('ms_in_camels_format/', d,
                                       '/daymet_full_climate.feather')) %>%
        # select(-pet) %>%
        distinct(date, site_code,
                 .keep_all = TRUE)

    if(d == 'arctic') ms_clim_ts0 <- filter(ms_clim_ts0, ! site_code %in% arctic_sites_to_drop)

    # apply(ms_clim_ts, 2, function(x) sum(is.na(x)) / length(x))

    ms_clim_ts <- ms_clim_ts0 %>%
        group_by(site_code) %>%
        filter(n() > 180) %>%
        ungroup() %>%
        bind_rows(ms_clim_ts)
}

## duplicate ts data for _extrapolate versions of neon sites
ms_clim_ts <- ms_clim_ts %>%
    filter(site_code %in% neon_sites) %>%
    mutate(site_code = paste0(site_code, '_extrapolate')) %>%
    bind_rows(ms_clim_ts)



## prepare Q
q_data <- load_product(
        macrosheds_root = paste0('../../data_acquisition/macrosheds_dataset_v',
                                 as.character(ms_version)),
        prodname = 'discharge'
    ) %>%
    filter(domain %in% ms_domains_included,
           substr(var, 1, 1) == 'I', #only sensor/weir Q measurements for now
           ! site_code %in% arctic_sites_to_drop)

if(length(unique(q_data$var)) > 1) stop('IN and IS discharge encountered. address this.')

q_data <- q_data %>%
    filter(ms_status == 0,
           ms_interp == 0,
           ! site_code %in% setdiff(unique(q_data$site_code),
                                    unique(ms_clim_ts$site_code))) %>%
    group_by(site_code) %>%
    filter(sum(! is.na(val)) > 157) %>%
    ungroup() %>%
    mutate(date = as.Date(datetime),
           val = errors::drop_errors(val),
           val = ifelse(val < 0, 0, val)) %>%
    select(date, site_code, discharge = val) %>%
    arrange(site_code, date) %>%
    group_by(site_code) %>%
    tidyr::complete(date = seq(min(date),
                               max(date),
                               by = '1 day')) %>%
    ungroup() %>%
    group_split(site_code)

# zz = purrr::map_dfr(q_data, bind_rows)
# sum(is.na(zz$discharge))

# #fill small Q gaps with seadec (naw, too risky)
# for(i in seq_along(q_data)){
#
#     qs <- q_data[[i]]
#     q_firstdate <- as.Date(qs$date[1])
#     q_startdoy <- as.numeric(format(q_firstdate, format = '%j'))
#     q_startyr <- year(q_firstdate)
#
#     q_data[[i]]$discharge <- tryCatch({
#
#          # forecast::findfrequency(na_interpolation(qs$discharge))
#         na_seadec(ts(qs$discharge,
#                      start = c(q_startyr, q_startdoy),
#                      frequency = 365.25),
#                   maxgap = 3) %>%
#             as.numeric()
#
#         # nas = is.na(qs$discharge)
#         # ts(qs$discharge,
#         #     start = c(q_startyr, q_startdoy),
#         #     frequency = 365.25)->qq
#         # plot(decompose(na_interpolation(qq)))
#         #
#         #  na_seadec(qs$discharge,
#         #            find_frequency = TRUE,
#         #            maxgap = 3) %>%
#         #     as.numeric()->zz
#         #
#         # dygraphs::dygraph(xts::xts(bind_cols(zz, select(qs, discharge)),
#         #                                  order.by = qs$date)) %>%
#         #     dygraphs::dyRangeSelector()
#
#
#         }, warning = function(w) na_interpolation(qs$discharge, maxgap = 3))
# }

## clean neon Q according to Spencer-Amanda-Gubbins evaluation
q_names <- sapply(q_data, function(x) x$site_code[1])

# nswnn <- read_csv('/home/mike/git/macrosheds/qa_experimentation/data/informative_stuff/ms_basins_with_netcdfs.csv') %>%
#     filter(domain == 'neon') %>%
#     pull(site_code)
nswnn = c('ARIK','BIGC','BLDE','BLUE','BLWA','CARI','COMO','CUPE','FLNT','GUIL','HOPB','KING','LECO','LEWI','MART','MAYF','MCDI','MCRA','OKSR','POSE','PRIN','REDB','SYCA','TECR','TOMB','WALK','WLOU')
neon_sites <- read_csv('informative_stuff/neon_basin_areas.csv')$site_code
q_eval <- read_csv('../neon_discharge_eval/data/neon_q_eval_final.csv') %>%
    filter(site %in% nswnn)

check1 <- is.na(q_eval$regression_status) | q_eval$regression_status %in% c('good', 'fair')
check2 <- is.na(q_eval$drift_status) | q_eval$drift_status %in% c('likely_no_drift', 'not_assessed')
check3 <- is.na(q_eval$rating_curve_status) | q_eval$rating_curve_status %in% c('Tier1', 'Tier2')

neon_sites_remaining <- q_eval %>%
    filter(site %in% nswnn) %>%
    filter(check1, check2, check3) %>%
    pull(site)
setdiff(nswnn, unique(neon_sites_remaining))
table(neon_sites_remaining)
# sort(table(neon_sites_remaining))

# exsite = 'BLUE'
# qq = filter(q_eval, site == exsite) %>%
#     select(-site, -curveID, -ends_with('status'),
#            -rating_curve_above_range, -rating_curve_under_range, -drift_mean_uncertaintiy) %>%
#     arrange(year, month)
# filter(q_eval, site == exsite) %>% arrange(year, month) %>% select(ends_with('status')) %>% print(n = 50)
# print(qq, n = 100)
# paste(sum(neon_sites_remaining == exsite), '/', sum(q_eval$site == exsite), 'for', exsite)

q_eval$keep <- check1 & check2 & check3
q_eval <- q_eval %>%
    group_by(site, year, month) %>%
    summarize(keep = any(keep), #not being strict here, since so few instances and plenty of visual subsetting anyway
              .groups = 'drop')

original_neon_q = list()
for(i in seq_along(neon_sites)){
    nsit = neon_sites[i]
    ni <- which(q_names == nsit)

    qqq = q_data[[ni]] %>%
        mutate(year = year(date),
               month = month(date))

    original_neon_q[[i]] = select(qqq, date, site_code, discharge)

    qqq = qqq %>%
        left_join(select(q_eval, site_code = site, year, month, keep),
                  by = c('site_code', 'year', 'month'))

    if(! nsit %in% c('GUIL', 'WLOU', 'BLDE', 'TOMB')){
        qqq = mutate(qqq, discharge = ifelse(keep, discharge, NA))
    }

    # if(nsit == 'BIGC') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2019-10-01'), discharge, NA))
    if(nsit == 'CARI') qqq$discharge[qqq$date > as.Date('2020-02-23') & qqq$date < as.Date('2020-03-31')] = NA
    if(nsit == 'CUPE') qqq = mutate(qqq, discharge = ifelse(date >= as.Date('2020-01-01'), discharge, NA))
    if(nsit == 'KING') qqq$discharge[qqq$date < as.Date('2018-12-01') | qqq$date > as.Date('2019-05-06')] = NA
    if(nsit == 'LEWI') qqq$discharge[qqq$date > as.Date('2020-03-05') & qqq$date < as.Date('2020-07-25')] = NA
    if(nsit == 'COMO') qqq$discharge[qqq$date > as.Date('2021-01-16') & qqq$date < as.Date('2021-03-06')] = NA
    # if(nsit == 'MAYF') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2019-07-19'), discharge, NA))
    if(nsit == 'MCDI') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2020-07-13'), discharge, NA))
    if(nsit == 'MCRA') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2019-05-06'), discharge, NA))
    # if(nsit == 'POSE') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2019-08-05'), discharge, NA))
    # if(nsit == 'GUIL') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2019-10-01'), discharge, NA))
    if(nsit == 'BLDE') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2020-10-01'), discharge, NA))
    # if(nsit == 'HOPB') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2020-12-28'), discharge, NA)) #it's actually P that's messed up here


    # print(nsit)
    # print(dygraphs::dygraph(xts(select(qqq, original_neon_q[[i]]$discharge, discharge), order.by = qqq$date)))
    # readLines(n = 1)

    q_data[[ni]] = select(qqq, site_code, date, discharge)
}
names(original_neon_q) = sapply(original_neon_q, function(x) x$site_code[1])

## manually clean a few sites (leading and trailing NAs will get trimmed later)

zzx <- which(q_names == 'GSCC01')
# zz = q_data[[zzx]]
# dygraphs::dygraph(xts::xts(select(zz, discharge),
#                                  order.by = zz$date)) %>%
#     dygraphs::dyRangeSelector()
q_data[[zzx]] <- filter(q_data[[zzx]], date > as.Date('1990-01-01'))
zzx <- which(q_names == 'GSCC02')
q_data[[zzx]] <- filter(q_data[[zzx]], date > as.Date('1990-01-01'))
zzx <- which(q_names == 'GSCC03')
q_data[[zzx]] <- filter(q_data[[zzx]], date > as.Date('1990-01-01'))
zzx <- which(q_names == 'GSCC04')
q_data[[zzx]] <- filter(q_data[[zzx]], date > as.Date('1990-01-01'))
# zzx <- which(q_names == 'WS79')
# q_data[[zzx]] <- filter(q_data[[zzx]], date > as.Date('2000-01-01'))
zzx <- which(q_names == 'QS')
q_data[[zzx]] <- filter(q_data[[zzx]], date > as.Date('2001-01-01'))
zzx <- which(q_names == 'GFVN')
q_data[[zzx]] <- filter(q_data[[zzx]], date > as.Date('1990-01-01'))

## replace some neon sites' Q series with those from nearby USGS gages ####

#scenarios:
#A: colocated sites; just fill neon data with the other gauge
#B: nearby site on different reach; predict neon site via regression.
    #if possible, use seasonal and era factors
#C: bracketing gauges; interpolate the neon site
#D: bracketing gauges, but one of them needs to be predicted first
#E: nearby gauge, but needs to be predicted first
#F: multiple nearby gauges: multiple regression

if(! length(list.files('../imputation/data/neon_continuous_Q/'))){
    get_neon_inst_discharge(neon_sites)
}

if(! length(list.files('../imputation/data/neon_field_Q/'))){
    get_neon_field_discharge(neon_sites)
}

dir.create('../imputation/out/q_lm_plots', showWarnings = FALSE)
dir.create('../imputation/out/q_lm_plots/diag', showWarnings = FALSE)
dir.create('../imputation/out/q_lm_plots/pred', showWarnings = FALSE)
dir.create('../imputation/out/q_lm_plots/fit', showWarnings = FALSE)
dir.create('../imputation/out/q_lm_plots/val', showWarnings = FALSE)

dir.create('../imputation/out/q_lm_outdata', showWarnings = FALSE)
dir.create('../imputation/out/q_lm_outdata/predictions', showWarnings = FALSE)
dir.create('../imputation/out/q_lm_outdata/fit', showWarnings = FALSE)
dir.create('../imputation/out/q_lm_outdata/summary', showWarnings = FALSE)

# q_data_neon_lm = list()
# formulaA = 'discharge ~ `{paste(gagenums, collapse = "`+`")}` + season'
formulaB = 'discharge_log ~ `{paste(paste0(gagenums, "_log"), collapse = "`+`")}` + season'
results = tibble(site_code = neon_sites, nse_logq = NA, nse_cv_logq = NA, bestmod_logq = NA,
                 adj_r_squared = NA)

# REDB (scenario B)
neon_site = 'REDB'; gagenums = '10172200'
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
lm_df = filter(lm_df, as.Date(datetime) != as.Date('2019-06-12'))                  #oi
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df, interactions = TRUE, through_origin = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

# HOPB (scenario B)
neon_site = 'HOPB'; gagenums = c('01174565')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df, interactions = TRUE, through_origin = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

# KING (scenario B)
# neon_site = 'KING'; usgsgagenums = c('06879650', '06879810', '06879100'); msgagenums = c('N01B', 'N04D', 'N20B', 'N02B'); gagenums = c(usgsgagenums, msgagenums)
# lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = usgsgagenums, nearby_ms_gages = msgagenums)
neon_site = 'KING'; gagenums = c('06879650', '06879810', '06879100', '06878600');
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# lm_df$`06879650`[lm_df$`06879650` > 1200 & ! is.na(lm_df$discharge)] = NA
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

# BLUE (scenario B)
neon_site = 'BLUE'; gagenums = '07332390'
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
out = plots_and_results(neon_site, best, lm_df, results, return_plot = T)
results = out$results

# CUPE (scenario B; multiple regression)
neon_site = 'CUPE'; gagenums = c('50136400', '50138000', '50144000') #'50128907' no area available (canal)
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

# GUIL (scenario B; multiple regression)
neon_site = 'GUIL'; gagenums = c('50028000', '50024950', '50126150', '50026025')#, '50010500')#, '50021700', '50114900')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
                         # include_sensor_daterange = c('2018-04-01', '2020-03-01'))
                         # include_sensor_daterange = c('2018-04-01', '2019-11-01'))
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE,
    min_points_per_param = 15,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

# # TECR (scenario B; using KREW site)
# neon_site = 'TECR'; gagenums = c('11216400') #gagenums = 'T003'
# lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# mods = generate_nested_formulae(
#     full_spec = as.formula(glue(formulaB)),
#     d = lm_df,
#     interactions = TRUE)
# best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
# results = plots_and_results(neon_site, best, lm_df, results)

#SYCA (as scenario B)

# neon_site = 'SYCA'; gagenums = c('09510200', '09499000') 09510150
neon_site = 'SYCA'; gagenums = c('09510200')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
                         # include_sensor_daterange = c('2019-06-01', '2020-07-01'))
# ggplot(lm_df, aes(x=`09510200_log`, y = `discharge_log`)) + geom_point()
lm_df = filter(lm_df, ! as.Date(datetime) == as.Date('2021-07-26'))

mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE,
    min_points_per_param = 15,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

# lm_df$x1 = lm_df$`09510200_log`
# m = lm(discharge_log ~ x1, data = lm_df)
# m_piecewise = segmented(m, seg.Z = ~x1, npsi=1)
# lm_df$best_prediction = inv_neglog(predict(m_piecewise, newdata = select(lm_df, x1)))
# lm_df$best_prediction[lm_df$best_prediction < 0] = 0
# # nse_syca = hydroGOF::NSE(lm_df$best_prediction, lm_df$discharge)
#
# plot_data = lm_df
# first_non_na = Position(function(x) ! is.na(x), plot_data$best_prediction)
# last_non_na = nrow(plot_data) - Position(function(x) ! is.na(x), rev(plot_data$best_prediction)) + 1
# plot_data = plot_data[first_non_na:last_non_na, ]
# plot_data = select(plot_data, -ends_with('_log'), -x1) %>%
#     select(site_code, datetime, Q_neon_field = discharge, Q_predicted = best_prediction,
#            everything())
#
# best = list()
# best$score = hydroGOF::NSE(lm_df$best_prediction, lm_df$discharge)
# best$score_crossval = NA_real_
# best$fits = NA
# best$best_model = discharge_log ~ `09510200_log`
# best$best_model_object = m_piecewise
# best$prediction = unname(lm_df$best_prediction)
# best$lm_data = plot_data
#
# results = plots_and_results(neon_site, best, lm_df, results)
#
# # plot_data = select(lm_df, datetime, Q_predicted = best_prediction, Q_used_in_regression = discharge,
# #        Q_neon_continuous_filtered = discharge_neon_cont,
# #        Q_neon_continuous_raw = discharge_neon_orig, Q_neon_manual = discharge_manual_forreals)
# #
# # dg = dygraphs::dygraph(xts(x = select(plot_data, -date), order.by = plot_data$date)) %>%
# #     dyRangeSelector()
# # saveWidget(dg, glue('../imputation/out/q_lm_plots/pred/{neon_site}_log.html'))
#
# png(glue('../imputation/out/q_lm_plots/fit/{neon_site}_fit_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
# plot(lm_df$`09510200_log`, lm_df$discharge_log, xlab = '09510200_log', ylab = 'discharge_log')
# plot(m_piecewise, add = TRUE, term = 'x1', col = 'red', lwd = 2)
# dev.off()
#
# # results$bestmod_logq[results$site_code == neon_site] = gsub('x1', '09510200_log', paste(as.character(formula(m_piecewise))[c(2, 1, 3)], collapse = ' '))
# # results$nse_logq[results$site_code == neon_site] = nse_syca
#
# # first_non_na = Position(function(x) ! is.na(x), lm_df$best_prediction)
# # last_non_na = nrow(lm_df) - Position(function(x) ! is.na(x), rev(lm_df$best_prediction)) + 1
# # lm_df = lm_df[first_non_na:last_non_na, ] %>%
# #     select(site_code, date,
# #            discharge_predicted = best_prediction,
# #            discharge_neon_manual = discharge_manual_forreals,
# #            discharge_neon_continuous_filtered = discharge_neon_cont,
# #            discharge_neon_continuous_raw = discharge_neon_orig,
# #            contains('09510200'))
# # write_csv(lm_df, glue('../imputation/out/q_lm_outdata/by_site/{neon_site}.csv'))
#
# class(m_piecewise) = 'lm'
# png(glue('../imputation/out/q_lm_plots/diag/{neon_site}_diag_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
# defpar = par(mfrow=c(2,2))
# plot(m_piecewise)
# par(defpar)
# dev.off()

# WALK (scenario B) NO UP-TO-DATE REFERENCE Q
neon_site = 'WALK'; gagenums = c('03535000', '03535400', '03495405') #gagenums = c('east_fork', 'west_fork');
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

#TOMB (scenario B)
# neon_site = 'TOMB'; gagenums = '02469525'
neon_site = 'TOMB'; gagenums = c('02469761', '02469525')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

#ARIK (scenario B)
neon_site = 'ARIK'; gagenums = c('06827000', '06823000')#, '06821500' produces rank deficiencies
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
                         # include_sensor_data = TRUE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

#MCDI (scenario B)
neon_site = 'MCDI'; gagenums = c('06888500', '06879650')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
                         # include_pct_highvals = 5)
                         # include_sensor_daterange = c('2018-11-01', '2019-10-01'))
# ggplot(lm_df, aes(x=`06888500_log`, y = `discharge_log`)) + geom_point()
# ggplot(lm_df, aes(x=`06879650_log`, y = `discharge_log`)) + geom_point()
# lm_df = filter(lm_df, ! datetime %in% ymd_hms(c('2018-08-27 14:19:00', '2017-07-05 13:19:00', '2018-09-26 15:30:00')))
lm_df = filter(lm_df, ! datetime %in% ymd_hms(c('2018-08-27 14:19:00', '2022-03-23 19:00:00')))

mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

# lm_df$x1 = lm_df$`06888500_log`
# m = lm(discharge_log ~ x1, data = lm_df)
# m_piecewise = segmented(m, seg.Z = ~x1, npsi = 1)
# lm_df$best_prediction = inv_neglog(predict(m_piecewise, newdata = select(lm_df, x1)))
# lm_df$best_prediction[lm_df$best_prediction < 0] = 0
# # nse_mcdi = hydroGOF::NSE(lm_df$best_prediction, lm_df$discharge)
#
# plot_data = lm_df
# first_non_na = Position(function(x) ! is.na(x), plot_data$best_prediction)
# last_non_na = nrow(plot_data) - Position(function(x) ! is.na(x), rev(plot_data$best_prediction)) + 1
# plot_data = plot_data[first_non_na:last_non_na, ]
# plot_data = select(plot_data, -ends_with('_log'), -x1) %>%
#     select(site_code, datetime, Q_neon_field = discharge, Q_predicted = best_prediction,
#            everything())
#
# best = list()
# best$score = hydroGOF::NSE(lm_df$best_prediction, lm_df$discharge)
# best$score_crossval = NA_real_
# best$fits = NA
# best$best_model = discharge_log ~ `06888500_log`
# best$best_model_object = m_piecewise
# best$prediction = unname(lm_df$best_prediction)
# best$lm_data = plot_data
#
# results = plots_and_results(neon_site, best, lm_df, results)
# # plot_data = select(lm_df, date, Q_predicted = best_prediction, Q_used_in_regression = discharge,
# #                    Q_neon_continuous_filtered = discharge_neon_cont,
# #                    Q_neon_continuous_raw = discharge_neon_orig, Q_neon_manual = discharge_manual_forreals)
# #
# # dg = dygraphs::dygraph(xts(x = select(plot_data, -date), order.by = plot_data$date)) %>%
# #     dyRangeSelector()
# # saveWidget(dg, glue('../imputation/out/q_lm_plots/pred/{neon_site}_log.html'))
#
# png(glue('../imputation/out/q_lm_plots/fit/{neon_site}_fit_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
# plot(lm_df$`06888500_log`, lm_df$discharge_log, xlab = '06888500_log', ylab = 'discharge_log')
# plot(m_piecewise, add = TRUE, term = 'x1', col = 'red', lwd = 2)
# dev.off()
#
# # results$bestmod_logq[results$site_code == neon_site] = gsub('x1', '06888500_log', paste(as.character(formula(m_piecewise))[c(2, 1, 3)], collapse = ' '))
# # results$nse_logq[results$site_code == neon_site] = nse_mcdi
#
# # first_non_na = Position(function(x) ! is.na(x), lm_df$best_prediction)
# # last_non_na = nrow(lm_df) - Position(function(x) ! is.na(x), rev(lm_df$best_prediction)) + 1
# # lm_df = lm_df[first_non_na:last_non_na, ] %>%
# #     select(site_code, date,
# #            discharge_predicted = best_prediction,
# #            discharge_neon_manual = discharge_manual_forreals,
# #            discharge_neon_continuous_filtered = discharge_neon_cont,
# #            discharge_neon_continuous_raw = discharge_neon_orig,
# #            contains('06888500'))
# # write_csv(lm_df, glue('../imputation/out/q_lm_outdata/by_site/{neon_site}.csv'))
#
# class(m_piecewise) = 'lm'
# png(glue('../imputation/out/q_lm_plots/diag/{neon_site}_diag_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
# defpar = par(mfrow=c(2,2))
# plot(m_piecewise)
# par(defpar)
# dev.off()


#LECO (scenario B)
neon_site = 'LECO'; gagenums = '03497300'
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

#LEWI (scenario B)
neon_site = 'LEWI'; gagenums = c('01636316', '01616100', '01636464')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

#PRIN (scenario B) MAYBE SOME KIND OF INTERP THING IS POSSIBLE?
# neon_site = 'PRIN'; gagenums = c('08044000', '08042950', '08042800')
neon_site = 'PRIN'; gagenums = c('08044000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

#POSE (scenario B)
neon_site = 'POSE'; gagenums = c('01636316', '01616100', '01662800', '01636464')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
                         # include_pct_highvals = 10)
                         # include_sensor_daterange = c('2017-01-25', '2018-08-01'))
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE,
    min_points_per_param = 15,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

#BLDE (scenario B)
neon_site = 'BLDE'; gagenums = c('06190540', '06188000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
                         # include_sensor_daterange = c('2019-01-01', '2020-09-01'))
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

#BLWA (scenario B)
neon_site = 'BLWA'; gagenums = c('02466030', '02465000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

#MCRA (scenario B)
# download.file('https://portal.edirepository.org/nis/dataviewer?packageid=knb-lter-and.4341.33&entityid=86490799297ff361b0741b807804c43a',
#               destfile = '../imputation/data/macrosheds_continuous_Q/hjandrews_raw.txt')
gagenums = c('GSWS01', 'GSWS06', 'GSWS07', 'GSWS08', 'GSLOOK')
hja = read_csv('../imputation/data/macrosheds_continuous_Q/hjandrews_raw.txt') %>%
    filter(SITECODE %in% gagenums) %>%
    select(datetime = DATE_TIME, SITECODE, INST_Q) %>%
    mutate(INST_Q = INST_Q * 28.317) %>%
    pivot_wider(names_from = SITECODE, values_from = INST_Q) %>%
    arrange(datetime)
neon_site = 'MCRA'; #gagenums = 'GSLOOK'; gagenums2 = c('')
lm_df = assemble_q_lm_df(neon_site = neon_site, ms_Q_data = hja, overwrite = T)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

# #BIGC (scenario B) - needs KREW update (using daily mode)
# neon_site = 'BIGC'; gagenums = '11237500' #gagenums = c('P300', 'P301', 'P304') #gagenums = c('11238250')
# lm_df = assemble_q_lm_df_daily(neon_site = neon_site, nearby_usgs_gages = gagenums)
# mods = generate_nested_formulae(
#     full_spec = as.formula(glue(formulaB)),
#     d = lm_df,
#     interactions = TRUE)
# best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
# results = plots_and_results_daily(neon_site, best, lm_df, results) #STOP: this func needs update

#MAYF (scenario B)
neon_site = 'MAYF'; gagenums = c('02465493', '02465292', '02424000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

#OKSR (scenario B)
neon_site = 'OKSR'; gagenums = c('15905100', '15908000', '15564879', '15875000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

#CARI (scenario B)
neon_site = 'CARI'; gagenums = c('15514000', '15511000', '15493400')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

#COMO (scenario B; composite model)
formulaC = "discharge_log ~ `{paste(paste0(gagenums, \"_log\"), collapse = \"`+`\")}`"

neon_site = 'COMO'; gagenums = c('ALBION', 'SADDLE')
ms_d = Filter(function(x) x$site_code[1] %in% gagenums, q_data) %>%
    reduce(bind_rows) %>%
    pivot_wider(names_from = site_code, values_from = discharge) %>%
    arrange(date)
lm_df1 = assemble_q_lm_df_daily(neon_site = neon_site, ms_Q_data = ms_d)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaC)),
    d = lm_df1,
    interactions = TRUE)
best1 = eval_model_set(data = lm_df1, model_list = mods, metric = hydroGOF::NSE)

neon_site = 'COMO'; gagenums = c('ALBION', 'SADDLE', 'MARTINELLI')
ms_d = Filter(function(x) x$site_code[1] %in% gagenums, q_data) %>%
    reduce(bind_rows) %>%
    pivot_wider(names_from = site_code, values_from = discharge) %>%
    arrange(date)
lm_df2 = assemble_q_lm_df_daily(neon_site = neon_site, ms_Q_data = ms_d)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaC)),
    d = lm_df2,
    interactions = TRUE)
best2 = eval_model_set(data = lm_df2, model_list = mods, metric = hydroGOF::NSE)

results = plots_and_results_daily_composite(neon_site, best1, best2, lm_df1, lm_df2, results)

#FLNT (scenario B)
neon_site = 'FLNT'; gagenums = c('02355662', '02353000', '02356000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

# MART (scenario B)
neon_site = 'MART'; gagenums = c('14138870', '14123500', '14120000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)#, include_pct_highvals = 10)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

# # WLOU (scenario B)
# neon_site = 'WLOU'; gagenums = c('09026500', '09025300', '09027100', '09034900')
# lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# mods = generate_nested_formulae(
#     full_spec = as.formula(glue(formulaB)),
#     d = lm_df,
#     interactions = TRUE)
# best = eval_model_set(data = lm_df, model_list = mods, metric = hydroGOF::NSE)
# results = plots_and_results(neon_site, best, lm_df, results)

#save output

# q_data_neon_lm = lapply(q_data_neon_lm, function(x) x$site_code = x$site_code[1])
# qq = reduce(q_data_neon_lm, bind_rows)
write_csv(results, '../imputation/out/q_lm_outdata/results.csv')


## TEMP: include non-MS Q data for NC-NHC and NC-UNHC ####
nc_q <- read_csv('ms_in_camels_format/NC/NHC_UNHC_good_flow_years_for_Q_estimation.csv') #rm this file too

ms_areas = bind_rows(ms_areas, tibble(site_code = c('NHC', 'UNHC'),
                                      ws_area_ha = c(8106.622487, 5801.544528)))

nc_q <- nc_q %>%
    mutate(date = as.Date(DateTime_UTC)) %>%
    group_by(date) %>%
    summarize(across(starts_with('discharge'),
                     ~mean(., na.rm = TRUE)),
              .groups = 'drop') %>%
    mutate(across(starts_with('discharge'),
                  ~(. * 1000)))

q_data[[length(q_data) + 1]] <- nc_q %>%
    select(date, discharge_nhc) %>%
    arrange(date) %>%
    tidyr::complete(date = seq(min(date), max(date), 'day')) %>%
    mutate(site_code = 'NHC') %>%
    select(site_code, date, discharge = discharge_nhc)

q_data[[length(q_data) + 1]] <- nc_q %>%
    select(date, discharge_unhc) %>%
    arrange(date) %>%
    tidyr::complete(date = seq(min(date), max(date), 'day')) %>%
    mutate(site_code = 'UNHC') %>%
    select(site_code, date, discharge = discharge_unhc)

## predict NHC and UNHC via lm ####

NHC_loc = tibble(lat = 35.9795, lon = -79.0018) %>%
    st_as_sf(coords = c('lon', 'lat'), crs = 4326)
UNHC_loc = tibble(lat = 35.9925, lon = -79.046) %>%
    st_as_sf(coords = c('lon', 'lat'), crs = 4326)
# mv = mapview::mapview
# mv(NHC_loc) + mv(UNHC_loc)

NHC = Filter(function(x) x$site_code[1] == 'NHC', q_data)[[1]] %>%
    arrange(date) %>%
    select(date, NHC = discharge)

site = 'NHC'; gagenums = c('02085000', '02097386', '0209722970', '0208675010')
neon_site = 'NHC' #needed
lm_df = assemble_q_lm_df_daily_ms(site = site, ms_Q_data = NHC, nearby_usgs_gages = gagenums, overwrite = T)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE,
    max_interaction = 3, min_points_per_param = 50)
best = eval_model_set(data = lm_df, model_list = mods[1:1000], metric = hydroGOF::NSE)
results = plots_and_results(neon_site, best, lm_df, results)

#more Q stuff ####
## include 365-day buffer of NAs before each discharge series, so the first Q
## value van be fully utilized/predicted
q_data <- lapply(q_data, function(x){
    first_q <- x$date[1]
    sc <- x$site_code[1]
    pre_buffer <- tibble(site_code = sc,
                         date = seq(first_q - 365, first_q - 1, 'day'),
                         discharge = NA_real_)
    x = bind_rows(pre_buffer, x)
    return(x)
})

#*^*
# plot(qq$date, qq$discharge, type = 'l')

## for neon sites, extend Q NAs back to 1980 and up to present day (full range of desired simulatio)
neon_sites_remaining <- unique(neon_sites_remaining)
site_order <- sapply(q_data, function(x) x$site_code[1])
neon_sites_extrapolate <- q_data[site_order %in% neon_sites_remaining]

#also modify the names of these extrapolated neon sites, so they get separate netcdfs
neon_sites_extrapolate <- lapply(neon_sites_extrapolate, function(x){
    mutate(x, site_code = paste0(site_code, '_extrapolate'))
})

neon_sites_extrapolate <- lapply(neon_sites_extrapolate, function(x){

    sc <- x$site_code[1]
    first_date <- x$date[1]
    historic_fill <- tibble(site_code = sc,
                            date = seq(as.Date('1980-01-01'), Sys.Date(), 'day'))
    x <- x %>%
        full_join(historic_fill, by = c('site_code', 'date')) %>%
        arrange(date)

    return(x)
})

q_data <- append(q_data, neon_sites_extrapolate)

#combine list back to one tibble
q_data <- purrr::map_dfr(q_data, bind_rows)

#duplicate neon watershed areas for the _extrapolate sites
extrap_areas <- ms_areas %>%
    filter(site_code %in% neon_sites_remaining) %>%
    mutate(site_code = paste0(site_code, '_extrapolate'))

ms_areas <- bind_rows(ms_areas, extrap_areas)

## convert to runoff
q_data <- q_data %>%
    left_join(ms_areas, by = 'site_code') %>%
    filter(! is.na(ws_area_ha)) %>%
    #      mm/d        L/s         m^3/L   ha/m^2 mm/m   s/d     ha
    mutate(discharge = discharge * 0.001 * 1e-4 * 1000 * 86400 / ws_area_ha) %>%
    # mutate(discharge = discharge * 0.0353146667) %>% #L/s to cfs
    select(-ws_area_ha)

ms_clim_ts <- filter(ms_clim_ts,
                     ! site_code %in% setdiff(unique(ms_clim_ts$site_code),
                                              unique(q_data$site_code)))

#for best shot at predicting FLNT, doing some extra interpolation.
#DISABLE THIS before final run
# flnt_bool <- q_data$site_code == 'FLNT'
#
# q_firstdate <- as.Date(q_data[flnt_bool, ]$date[1])
# q_startdoy <- as.numeric(format(q_firstdate, format = '%j'))
# q_startyr <- year(q_firstdate)
#
# q_data[flnt_bool, 'discharge'] <-
#     na_seadec(ts(pull(q_data[flnt_bool, 'discharge']),
#                  start = c(!!q_startyr, !!q_startdoy),
#                  deltat = 1/365),
#               maxgap = 11)

# zz = filter(q_data, site_code == 'WALK') %>% mutate(val = drop_errors(val))
# sum(is.na(zz$val))
# sum(zz$val == -999, na.rm = TRUE)

## visualize Q
# q_for_viz <- q_data %>%
#     select(date, site_code, discharge) %>%
#     pivot_wider(names_from = site_code,
#                 values_from = discharge)
#
# dygraphs::dygraph(xts(x = select(q_for_viz, -date),
#                       order.by = q_for_viz$date)) %>%
#     dyRangeSelector()

## one more Q thing (NetCDF needs numeric dates)
q_data <- mutate(q_data, date = as.double(date))

## prepare chem
# chempreds <- predictors[! predictors == 'precipitation']
# if(length(chempreds)){
#     chem <- load_product(prodname = 'stream_chemistry',
#                          location = 'macrosheds_dataset',
#                          filter_domains = rep(domain, length(site_code)),
#                          filter_site_codes = site_code,
#                          filter_vars = chempreds)
#
#     if(any(substr(chem$var, 1, 1) == 'I')){
#         stop('installed chem product encountered. address this')
#     }
#
#     chem <- chem %>%
#         filter(extract_var_prefix(var) == 'GN') %>%
#         mutate(var = drop_var_prefix(var)) %>%
#         filter(var %in% chempreds) %>%
#         select(datetime, site_code, var, val) %>%
#         group_by(site_code, var) %>%
#         tidyr::complete(datetime = seq(min(datetime),
#                                        max(datetime),
#                                        by = '1 day')) %>%
#         ungroup() %>%
#         pivot_wider(names_from = var,
#                     values_from = val) %>%
#         mutate(date = as.double(as.Date(datetime))) %>%
#         select(-datetime)
# }

## prepare precip
# p <- load_product(prodname = 'precipitation',
#                   location = 'macrosheds_dataset',
#                   filter_domains = rep(domain, length(site_code)),
#                   filter_site_codes = site_code) %>%
#     filter(substr(var, 1, 1) == 'I') #only gauge P measurements for now
#
# if(length(unique(p$var)) > 1) stop('IN and IS precip encountered. address this.')
#
# p <- p %>%
#     filter(var == 'IN_precipitation') %>%
#     tidyr::complete(datetime = seq(min(datetime),
#                                    max(datetime),
#                                    by = '1 day')) %>%
#     mutate(date = as.double(as.Date(datetime))) %>%
#     select(date, site_code,
#            precipitation = val)

## merge discharge and climate forcings timeseries
ms_clim_ts <- mutate(ms_clim_ts,
                     date = as.double(as.Date(date)))

ms_ts <- left_join(q_data, ms_clim_ts,
                   by = c('date', 'site_code')) %>%
    mutate(q_source = 'true')

# apply(ms_ts, 2, function(x) sum(is.na(x))/length(x))

# if(length(chempreds)){
#     d <- left_join(d, chem,
#                    by = c('date', 'site_code'))
# }

## split by site (each will be written to a separate NetCDF file)
ms_ts_splt <- group_split(ms_ts, site_code) %>%
    as.list()
# ms_ts_splt->rer
# ms_ts_splt<-rer

## drop sites with no forcing data
ms_ts_splt <- Filter(function(x){
    cmplt <- which(complete.cases(select(x, -discharge)))
    if(length(cmplt)) T else F
}, ms_ts_splt)
ms_ts = filter(ms_ts, site_code %in% sapply(ms_ts_splt, function(x) x$site_code[1]))

# filter(qqq, if_all(everything(), ~! is.na(.)))
# dygraphs::dygraph(xts::xts(select(qqq, -date),
#                            order.by = as.Date(qqq$date, origin = '1970-01-01'))) %>%
#     dygraphs::dyRangeSelector()

## drop leading/trailing runs of forcings NAs
ms_ts_splt <- lapply(ms_ts_splt, function(x){
    cmplt <- which(complete.cases(select(x, -discharge)))
    # print(x$site_code[1])
    x[min(cmplt):max(cmplt), ]
})

## check NAs; impute (forcings); leave discharge alone
site_codes <- unique(ms_ts$site_code)
chuck_these <- c()
chuck_these_names <- c()
for(i in seq_along(site_codes)){

    ms_chunk <- ms_ts_splt[[i]]
    site_code <- ms_chunk$site_code[1]

    na_proportions <- apply(select(ms_chunk, -date, -site_code),
                            MARGIN = 2,
                            FUN = function(x) sum(is.na(x)) / length(x))

    if(na_proportions[names(na_proportions) == 'dayl'] == 1){
        chuck_these <- c(chuck_these, i)
        chuck_these_names <- c(chuck_these_names, site_code)
        warning(paste('dropping', site_code, 'from ms_ts_splt for lack of daymet data'))
        next
    }

    # dygraphs::dygraph(xts(x = select(ms_chunk, -date, -site_code),
    #                       order.by = as.Date(ms_chunk$date,
    #                                          origin = '1970-01-01'))) %>%
    #     dyRangeSelector()

    q_nas <- na_proportions[names(na_proportions) == 'discharge']
    f_nas <- na_proportions[names(na_proportions) != 'discharge']

    #impute forcings NAs via linear interpolation with seasonal decomposition
    if(any(f_nas > 0)){

        print(glue('NA proportions for {s}:\n\t{props}',
                   s = site_code,
                   props = paste(names(f_nas),
                                 round(unname(f_nas), 4),
                                 sep = ': ',
                                 collapse = '\n\t'),
                   .trim = FALSE))

        # message('Replacing NAs with -999')
        #
        # ms_ts_splt[[i]] <- mutate(ms_chunk,
        #                           across(-all_of('date'),
        #                                  ~na_replace(., fill = -999)))

        message('Replacing forcing NAs via linear interpolation with seasonal decomposition')

        ts_firstdate <- as.Date(ms_chunk$date, origin = '1970-01-01')[1]
        ts_startdoy <- as.numeric(format(ts_firstdate, format = '%j'))
        ts_startyr <- year(ts_firstdate)

        ts_lastdate <- as.Date(ms_chunk$date, origin = '1970-01-01')[nrow(ms_chunk)]
        ts_enddoy <- as.numeric(format(ts_lastdate, format = '%j'))
        ts_endyr <- year(ts_lastdate)

        ms_chunk$pet[is.infinite(ms_chunk$pet)] <- NA

        tryCatch({
            ms_ts_splt[[i]] <- mutate(
                ms_chunk,
                across(-all_of(c('date', 'site_code', 'discharge', 'q_source')),
                       ~na_seadec(ts(.,
                                     start = c(!!ts_startyr, !!ts_startdoy),
                                     # end = c(!!ts_endyr, !!ts_enddoy),
                                     frequency = 365.25),
                                  maxgap = Inf)))

                # forecast::findfrequency(na_interpolation(ms_chunk$tmax))

                # dygraphs::dygraph(xts::xts(bind_cols(zz$tmax, ms_chunk$tmax),
                #                                  order.by = as.Date(zz$date))) %>%
                #     dygraphs::dyRangeSelector()
        }, warning = function(w){

            na_sums <- apply(select(ms_chunk, -date, -site_code, -discharge),
                             MARGIN = 2,
                             FUN = function(x) sum(is.na(x)))

            if(any(na_sums > 1)){
                stop(paste('must use a different interp algorithm for', i))
            }
        })
    }

    #replace discharge NAs with -999 (scratch that)
    if(q_nas > 0){

        print(glue('discharge NA proportion for {s}: {props}',
                   s = site_code,
                   props = round(unname(q_nas), 4),
                   .trim = FALSE))

        message('NOT replacing discharge NAs with -999')

        # ms_ts_splt[[i]] <- mutate(ms_ts_splt[[i]],
        #                           discharge = imputeTS::na_replace(discharge,
        #                                                            fill = -999))

        # ms_ts_splt[[i]] <- mutate(
        #     ms_chunk,
        #     across(all_of(c('discharge')),
        #            ~na_interpolation(., maxgap = Inf)))
    }
}

ms_ts_splt[chuck_these] = NULL

## make sure all the tables are sorted by date
ms_ts_splt <- lapply(ms_ts_splt, function(x) arrange(x, date))

### 3b OPTIONAL: introduce gaps in macrosheds Q ####

gap_site <- 'WALK'

site_order <- sapply(ms_ts_splt, function(x) x$site_code[1])

gap_sites_inds <- match(gap_site, site_order)

# for(i in seq_along(gap_sites_inds)){
gs <- ms_ts_splt[[gap_sites_inds]]

#save an unmodified copy
# ms_write_netcdf(df_list = ms_ts_splt[gap_sites_inds],
#                 path = 'unmodified_comparator_netcdfs')
write_csv(mutate(gs, date = as.Date(date, origin = '1970-01-01')),
          file.path('unmodified_comparator_datasets',
                    paste0(gap_site, '.csv')))

# library(dygraphs)
# dygraphs::dygraph(xts::xts(gs$discharge,
#                            order.by = as.Date(gs$date, origin = '1970-01-01'))) %>%
#     dygraphs::dyRangeSelector()

gap_ranges <- list(c('2017-04-22', '2017-08-28'),
                   c('2018-06-01', '2019-04-01'),
                   c('2019-07-19', '2020-01-10'))

# gs$discharge2 = gs$discharge
for(i in seq_along(gap_ranges)){
    gr <- as.numeric(as.Date(gap_ranges[[i]]))
    gs$discharge[gs$date > gr[1] & gs$date < gr[2]] <- NA
}

# dygraphs::dygraph(xts::xts(gs[, c('discharge', 'discharge2')],
dygraphs::dygraph(xts::xts(gs[, c('discharge')],
                           order.by = as.Date(gs$date, origin = '1970-01-01'))) %>%
    dygraphs::dyRangeSelector()

# ms_ts_splt[[gap_sites_inds]] <- gs #UNCOMMENT TO INTRODUCE GAP


### 3c prepare macrosheds timeseries data from NHM (optional plot) ####

# pdf(width = 12, height = 10, file = file.path('..', 'plots', 'nhm_vs_ms_Q.pdf'))
# par(mfrow = c(5, 2))

daymet_files <- paste0('ms_in_camels_format/',
                       unique(nhm_segids$domain),
                       '/daymet_full_climate.feather')
ms_daymet <- map_dfr(daymet_files, read_feather) %>%
    filter(site_code %in% unique(nhm_segids$site_code))

nhm_ms_ts_splt <- list()
for(i in 1:nrow(nhm_segids)){

    jobid <- pull(nhm_segids[i, 'JobId'])
    segid <- pull(nhm_segids[i, 'NHM_SEGID'])
    siteid <- pull(nhm_segids[i, 'site_code'])
    ws_area_ha <- pull(ms_areas[ms_areas$site_code == siteid, 'ws_area_ha'])

    nhm_q <- read_csv(glue('NHMv1/nhm_output_ms/basin_{b}/model_output/seg_outflow.csv',
                      b = str_pad(jobid,
                                  width = 4,
                                  side = 'left',
                                  pad = '0')),
                      show_col_types = FALSE) %>%
        select(date = time,
               discharge_nhm = !!segid) %>%
        # mutate(date_dbl = as.double(date)) %>%
        #      mm/d            ft^3/s          m^3/ft^3    ha/m^2 mm/m   s/d     ha
        mutate(discharge_nhm = discharge_nhm * 0.0283168 * 1e-4 * 1000 * 86400 / !!ws_area_ha)
        # mutate(discharge = discharge * 0.0353146667) %>% #L/s to cfs

    # #timeseries plot of Q
    # temp <- ms_ts %>%
    #     filter(site_code == siteid) %>%
    #     full_join(nhm_q, by = 'date') %>%
    #     select(date, discharge, discharge_nhm) %>%
    #     mutate(date = as.Date(date, origin = '1970-01-01'))
    #
    # plot(temp$date, temp$discharge, type = 'l', col = 'blue', ylab = 'runoff (mm/d)',
    #      xlim = range(temp$date, na.rm = TRUE),
    #      # xlim = as.Date(c('2016-10-01', '2017-02-01')),
    #      xlab = '', main = siteid,
    #      ylim = range(c(temp$discharge, temp$discharge_nhm), na.rm = TRUE))
    #      # ylim = c(0, 0.3))
    # lines(temp$date, temp$discharge_nhm, col = 'red')
    # legend('topright', legend = c('MS', 'NHM'), col = c('blue', 'red'), lty = 1, bty = 'n')
    #
    # temp <- temp[complete.cases(temp), ]
    #
    # #plot NHM Q vs MS Q (or nothing, if no coinciding points)
    # if(nrow(temp)){
    #     rng <- range(c(temp$discharge, temp$discharge_nhm), na.rm = TRUE)
    #     plot(temp$discharge, temp$discharge_nhm, ylab = 'NHM runoff (mm/d)',
    #          xlab = 'MS runoff (mm/d)', xlim = rng, ylim = rng)
    #     abline(0, 1, lty = 3)
    # } else {
    #     plot(1, 1, type = 'n', ann = FALSE, axes = FALSE)
    # }

    # neon_areas <- read_csv('../../portal/data/general/site_data.csv') %>%
    #     filter(site_type == 'stream_gauge',
    #            ! is.na(ws_area_ha),
    #            domain == 'neon') %>%
    #     select(site_code, ws_area_ha) %>%
    #     arrange(desc(ws_area_ha))
    # ff_true = read_feather('../../data_acquisition/data/neon/neon/derived/discharge__ms005/FLNT.feather') %>%
    #     left_join(ms_areas, by = 'site_code') %>%
    #     mutate(discharge_true = val * 1e-4 * 86400 / ws_area_ha) %>%
    #     rename(date = datetime)
    # ff_rodeo = read_csv('~/Downloads/ee-chart(1).csv') %>%
    #     mutate(date = lubridate::mdy(`system:time_start`),
    #            site_code = 'FLNT') %>%
    #     left_join(ms_areas, by = 'site_code') %>%
    #     mutate(Discharge = Discharge * 1e-4 * 86400 * 1000/ ws_area_ha) %>%
    #     select(date, discharge_rodeo = Discharge)

    nhm_ms_ts_part <- ms_daymet %>%
    # zz <- ms_daymet %>%
        filter(site_code == siteid) %>%
        left_join(nhm_q, by = 'date') %>%
        filter(! is.na(discharge_nhm)) %>%
        arrange(date) %>%
        # complete(date = seq(date[1], date[nrow(.)], by = 'day'))

    #     full_join(ff_rodeo, by = 'date') %>%
    #     full_join(ff_true, by = c('date', 'site_code')) %>%
    #     filter(! is.na(discharge_rodeo),
    #            ! is.na(discharge_nhm)) %>%
    #     mutate(dif = abs(discharge_nhm - discharge_rodeo),
    #            dif_pct = dif / ((discharge_nhm + discharge_rodeo) / 2) * 100)
    # mean(zz$dif, na.rm = TRUE)
    # mean(zz$dif_pct, na.rm = TRUE)
    #     # filter(date >= as.Date('2016-01-01'),
    #     #        date <= as.Date('2016-12-31')) %>%
    #     # filter(! is.na(discharge_true)) %>%
    #     select(date, discharge_nhm, discharge_rodeo, discharge_true)

        select(site_code, date, discharge = discharge_nhm, dayl, prcp, srad,
               swe, tmax, tmin, vp) %>%
        mutate(q_source = 'NHM',
               date = as.numeric(as.Date(date)),
               site_code = paste0('NHM_', site_code))

    not_shoulder_na <- ! cumprod(is.na(nhm_ms_ts_part$discharge)) &
        rev(! cumprod(is.na(rev(nhm_ms_ts_part$discharge))))
    nhm_ms_ts_part <- nhm_ms_ts_part[not_shoulder_na, ]

    nhm_ms_ts_splt[[i]] <- try({
    nhm_ms_ts_part %>%
        tidyr::complete(date = seq(min(date), max(date), 1)) %>%
        mutate(across(-all_of(c('date', 'site_code', 'discharge', 'q_source')),
                      ~na_interpolation(., maxgap = 3)),
               q_source = ifelse(is.na(q_source), 'NHM', q_source))
    }, silent = TRUE)

    if(inherits(nhm_ms_ts_splt[[i]], 'try-error')) nhm_ms_ts_splt[[i]] = NA
}

nhm_ms_ts_splt[is.na(nhm_ms_ts_splt)] = NULL

# dev.off()

### 4a prepare CAMELS timeseries data ####

## prepare usgs Q and daymet forcings
discharge_files <- list.files('CAMELS/basin_dataset_public_v1p2/usgs_streamflow',
                              full.names = TRUE,
                              recursive = TRUE)

daymet_files_agg <- list.files('CAMELS/basin_dataset_public_v1p2/basin_mean_forcing/daymet',
                               full.names = TRUE,
                               recursive = TRUE)
daymet_files_agg <- grep('011230\\*|208111310_forcing|344894205_forcing', #rm wonky files
                         daymet_files_agg,
                         value = TRUE,
                         invert = TRUE)

daymet_isolate_path = 'CAMELS_daymet_full/model_output/flow_timeseries/daymet'

#the aggregate CAMELS dataset is missing all SWE information, so pulling it from
#the daymet model output in the parallel loop below
daymet_files <- list.files(daymet_isolate_path,
                           pattern = 'model_output',
                           full.names = TRUE,
                           recursive = TRUE)

pet_files <- list.files('CAMELS_macrosheds_combined/camels_pet_isolate',
                        full.names = TRUE)

if(! file.exists('CAMELS_macrosheds_combined/cmls_ts_splt.rds')){

clst_type <- ifelse(.Platform$OS.type == 'windows', 'PSOCK', 'FORK')
ncores <- parallel::detectCores() %/% 1.5
clst <- parallel::makeCluster(spec = ncores,
                              type = clst_type)
doParallel::registerDoParallel(clst)

#parallelization fails here for some reason. trying to avoid thread load balance issues
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
for(kk in megaloop_list[19:22]){
# for(kk in megaloop_list){

    cnt <- cnt + 1
    print(paste(cnt, 'of', length(megaloop_list)))
    # cmls_ts_splt <- foreach(i = seq_along(discharge_files),
    cmls_ts_splt0 <- foreach(i = kk,
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
        huc_code = unique(str_match(ffs, paste0(daymet_isolate_path, '/([0-9]{2})/'))[, 2])[1]
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

        #also want to know if this pet is the same as what we calculate using gridded alpha parameter
        #(it's not identical, but close)

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
            mutate(date = as.double(ymd(paste(year, month, day)))) %>%
            select(-month, -day, -year, -hour)

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
}

# saveRDS(cmls_ts_splt, 'CAMELS_macrosheds_combined/partial_cmls_ts_splt2.rds')
# cmls_ts_splt00 = readRDS('CAMELS_macrosheds_combined/partial_cmls_ts_splt.rds')
# cmls_ts_splt01 = readRDS('CAMELS_macrosheds_combined/partial_cmls_ts_splt2.rds')
# cmls_ts_splt = append(cmls_ts_splt00, cmls_ts_splt)

# saveRDS(cmls_ts_splt, 'CAMELS_macrosheds_combined/cmls_ts_splt.rds')

parallel::stopCluster(clst)
} else {
    cmls_ts_splt = readRDS('CAMELS_macrosheds_combined/cmls_ts_splt.rds')
}

### 4b prepare CAMELS timeseries data from NHM ####

# mv=mapview::mapview

#each gauge output file includes estimated Q for all upstream reaches, and
#nowhere does the NHM indicate which segid corresponds to the gage. so, to find out,
#we get the location of the gage from the CAMELS dataset, then find the nearest
#feature within the NHM. verified on 01013500, and assumed to work fine for the
#other 530 gages.

# st_layers('NHMv1/GF_nat_reg.gdb')
camels_segids <- st_read('NHMv1/GF_nat_reg.gdb',
                         layer = 'nsegmentNationalIdentifier') %>%
    st_transform(crs = 4326)

camels_gauge_info_with_segids <- 'NHMv1/camels_gauge_info_with_segids.csv'
if(file.exists(camels_gauge_info_with_segids)){
    camels_gauge_info <- read_csv(camels_gauge_info_with_segids)
} else {

    camels_gauge_info <- read_tsv(
        '/home/mike/git/macrosheds/qa_experimentation/data/CAMELS/basin_dataset_public_v1p2/basin_metadata/gauge_information.txt',
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

        # site_buf <- sf::st_buffer(x = cmls_site,
        #                           dist = 10000)
        # site_box <- st_bbox(site_buf)
        # nhm_crop <- suppressWarnings(sf::st_crop(camels_segids, site_box))
        # mv(nhm_crop) + mv(cmls_site)
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
        read_csv(glue('NHMv1/nhm_output_camels/{gid}/model_output/seg_outflow.csv',
                               gid = cmls_gageid),
                 show_col_types = FALSE) %>%
            select(date = time,
                   discharge_nhm = !!cmls_segid) %>%
            #      mm/d            cfs             m^3/ft^3      mm/m   s/d     m^2
            mutate(discharge_nhm = discharge_nhm * 0.028316846 * 1000 * 86400 / cmls_ws_area_m2)
    }, silent = TRUE)

    if(inherits(nhm_q, 'try-error')){
        print(paste('error in', cmls_gageid))
        next
    }

    cmls_ts_part <- cmls_ts_part %>%
        mutate(date_as_date = as.Date(date)) %>%
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

### 5 investigate discharge ranges ####

cmls_names <- sapply(cmls_ts_splt, function(x) x$site_code[1])
ms_names <- sapply(ms_ts_splt, function(x) x$site_code[1])

cmls_ranges <- map_dfr(cmls_ts_splt,
    function(x){
        z = x$discharge
        z[z %in% c(-999, 999)] = NA_real_
        z = range(z, na.rm = TRUE)
        names(z) = c('min', 'max')
        return(z)
    })
ms_ranges <- map_dfr(ms_ts_splt,
    function(x){
        z = x$discharge
        z[z %in% c(-999, 999)] = NA_real_
        z = range(z, na.rm = TRUE)
        names(z) = c('min', 'max')
        return(z)
    })

cmls_Q_ranges <- bind_cols(tibble(site_code = cmls_names), cmls_ranges)
ms_Q_ranges <- bind_cols(tibble(site_code = ms_names), ms_ranges)

options(scipen = 100)
fivenum(cmls_Q_ranges$min)
fivenum(ms_Q_ranges$min)
fivenum(cmls_Q_ranges$max)
fivenum(ms_Q_ranges$max)
options(scipen = 0)

## compare watershed areas across CAMELS and MacroSheds

cmls_areas <- read.fwf(
    'CAMELS/basin_dataset_public_v1p2/basin_metadata/basin_physical_characteristics.txt',
    widths = c(3, 9, 11, 10, 11, 6),
    skip = 1,
    header = FALSE,
    strip.white = TRUE,
    colClasses = c('character', 'character', 'numeric', 'numeric', 'numeric', 'numeric'),
    col.names = c('huc', 'basin_id', 'area_km2', 'elev_m', 'slope_mkm', 'pct_forest')
) %>%
    as_tibble() %>%
    mutate(ws_area_ha = area_km2 * 100) %>%
    select(basin_id, ws_area_ha) %>%
    arrange(desc(ws_area_ha))

ms_areas <- arrange(ms_areas, desc(ws_area_ha))

range(cmls_areas$ws_area_ha)
range(ms_areas$ws_area_ha)

cmls_maxes <- sapply(cmls_ts_splt, function(x) max(x$discharge, na.rm = TRUE))
sum(cmls_maxes > max(q_data$discharge, na.rm = TRUE),
    na.rm = TRUE) #q_data is ms only
sum(cmls_maxes <= max(q_data$discharge, na.rm = TRUE),
    na.rm = TRUE) #q_data is ms only

#okay, there are some huge NEON basins, but we only have a few years of data for
#   them, so their max discharges are still much less than some of the CAMELS basins

## save some useful stuff
dir.create('informative_stuff', showWarnings = FALSE)
write_csv(cmls_areas, 'informative_stuff/camels_basin_areas.csv')
write_csv(ms_areas, 'informative_stuff/ms_basin_areas.csv')

### 6 write netCDF file with forcings and discharge for each basin ####

all_ts_splt <- append(cmls_ts_splt, ms_ts_splt)

all_ts_splt <- lapply(all_ts_splt,
                      function(x){
                          q_src_int <- q_source_map[[x$q_source[1]]]
                          x$q_source <- q_src_int
                          x
                      })

# vapply(all_ts_splt, function(x) length(unique(rle(diff(x$date))$lengths)), 1)

# yy = list.files('CAMELS_macrosheds_combined/time_series') %>%
#     str_replace('\\.nc', '')
# xx = sapply(all_ts_splt, function(x) x$site_code[1])
# all_ts_splt[[which(xx=='SYCA_extrapolate')]]
# zz = setdiff(yy, xx)
# zz[! grepl('^NHM_', zz)]
# setdiff(xx, yy)


ms_write_netcdf(df_list = all_ts_splt,
                path = 'CAMELS_macrosheds_combined/time_series')

## also for NHM pseudo-basins (process model output, rather than data)

nhm_all_ts_splt <- append(nhm_cmls_ts_splt, nhm_ms_ts_splt)

nhm_all_ts_splt <- lapply(nhm_all_ts_splt,
                          function(x){
                              q_src_int <- q_source_map[[x$q_source[1]]]
                              x$q_source <- q_src_int
                              x
                          })

# vapply(nhm_all_ts_splt, function(x) length(unique(rle(diff(x$date))$lengths)), 1)
# vapply(nhm_all_ts_splt, function(x) x$site_code[1], 'a')

ms_write_netcdf(df_list = nhm_all_ts_splt,
                path = 'CAMELS_macrosheds_combined/time_series')
### 6b MODIFY NETCDFs without re-running everything above ####

library(zoo)

ncdfs = list.files('CAMELS_macrosheds_combined/time_series/', full.names = TRUE)
ncdfs = grep('_GAPPED.nc$', ncdfs, invert = TRUE, value = TRUE)
ncdfs = grep('_MANUALQ.nc$', ncdfs, invert = TRUE, value = TRUE)

#remove modified netcdfs and replace them with a new set
# unlink('CAMELS_macrosheds_combined/time_series', recursive = TRUE)
# wd = getwd()
# system(paste0('cd ', wd, '; unzip CAMELS_macrosheds_combined/time_series.zip'))
# file.rename('~/git/macrosheds/qa_experimentation/data/time_series/',
#             file.path('~/git/macrosheds/qa_experimentation/data/CAMELS_macrosheds_combined/', 'time_series'))

#Qlags and Qlag binary indicators
for(i in seq_along(ncdfs)){

    f = ncdfs[i]

    xx = ncdf4::nc_open(f)
    dd = ncdf4::ncvar_get(xx, 'date')
    qq = ncdf4::ncvar_get(xx, 'discharge')
    Qlag0 = qq
    Qlag0[is.na(Qlag0)] = -999
    Qlag1 = lead(Qlag0, default = -999)
    Qlag2 = lead(Qlag0, 2, default = -999)
    Qlag3 = lead(Qlag0, 3, default = -999)
    Qlag5 = lead(Qlag0, 5, default = -999)
    Qlag7 = lead(Qlag0, 7, default = -999)
    Qlag10 = lead(Qlag0, 10, default = -999)
    Qbin0 = as.numeric(Qlag0 == -999)
    Qbin1 = as.numeric(Qlag1 == -999)
    Qbin2 = as.numeric(Qlag2 == -999)
    Qbin3 = as.numeric(Qlag3 == -999)
    Qbin5 = as.numeric(Qlag5 == -999)
    Qbin7 = as.numeric(Qlag7 == -999)
    Qbin10 = as.numeric(Qlag10 == -999)
    nc_close(xx)

    dim_time <- ncdim_def('date',
                          units = 'days since 1970-01-01 00:00',
                          calendar = 'standard',
                          vals = dd)

    xx = ncdf4::nc_open(f, write = TRUE)

    for(newvar_name in paste0('Qlag', c(0,1,2,3,5,7,10))){

        newvar = get(newvar_name)
        newvar_nc <- ncvar_def(name = newvar_name,
                               units = 'cms', #if you change this, change ms_write_netcdf and other instances
                               dim = list(dim_time),
                               missval = NULL)
        xx2 = ncvar_add(xx, newvar_nc)
        ncvar_put(xx2, newvar_name, newvar)
    }

    for(newvar_name in paste0('Qbin', c(0,1,2,3,5,7,10))){

        newvar = get(newvar_name)
        newvar_nc <- ncvar_def(name = newvar_name,
                               units = '',
                               dim = list(dim_time),
                               missval = NULL)
        xx2 = ncvar_add(xx, newvar_nc)
        ncvar_put(xx2, newvar_name, newvar)
    }

    nc_close(xx)
}

#thinned Q series
for(i in seq_along(ncdfs)){

    f = ncdfs[i]

    xx = ncdf4::nc_open(f)
    dd = ncdf4::ncvar_get(xx, 'date')
    qq = ncdf4::ncvar_get(xx, 'discharge')
    qq[is.na(qq)] = -999

    QalternatingNA = QrandomNA = qq

    nrecs = length(qq)
    alt_drop_locs = rep(c(F, T), length.out = nrecs)
    QalternatingNA[alt_drop_locs] = -999
    QalternatingNAind = QalternatingNA == -999

    random_drop_locs = sample(1:nrecs, nrecs %/% 2, replace = FALSE)
    QrandomNA[random_drop_locs] = -999
    QrandomNAind = QrandomNA == -999

    nc_close(xx)

    dim_time <- ncdim_def('date',
                          units = 'days since 1970-01-01 00:00',
                          calendar = 'standard',
                          vals = dd)

    xx = ncdf4::nc_open(f, write = TRUE)

    for(newvar_name in c('QalternatingNA', 'QrandomNA')){

        newvar = get(newvar_name)
        newvar_nc <- ncvar_def(name = newvar_name,
                               units = 'cms', #if you change this, change ms_write_netcdf and other instances
                               dim = list(dim_time),
                               missval = NULL)
        xx2 = ncvar_add(xx, newvar_nc)
        ncvar_put(xx2, newvar_name, newvar)
    }

    for(newvar_name in c('QalternatingNAind', 'QrandomNAind')){

        newvar = get(newvar_name)
        newvar_nc <- ncvar_def(name = newvar_name,
                               units = '',
                               dim = list(dim_time),
                               missval = NULL)
        xx2 = ncvar_add(xx, newvar_nc)
        ncvar_put(xx2, newvar_name, newvar)
    }

    nc_close(xx)
}

#moving average Q series
for(i in seq_along(ncdfs)){

    f = ncdfs[i]

    xx = ncdf4::nc_open(f)
    dd = ncdf4::ncvar_get(xx, 'date')
    qq = c(ncdf4::ncvar_get(xx, 'discharge'))

    Qma = qq
    na_inds = is.na(Qma)
    if(any(! na_inds)){
        Qma = imputeTS::na_locf(Qma)
        Qma = zoo::rollmean(Qma, 7, fill = 'extend')
    }
    Qma[na_inds] = -999
    QmaNAind = Qma == -999

    QmaExtraNA = Qma
    make_these_NA_too = sample(which(QmaExtraNA != -999), sum(! QmaNAind) * 0.1)
    QmaExtraNA[make_these_NA_too] = -999
    QmaExtraNAind = QmaExtraNA == -999

    nc_close(xx)

    dim_time <- ncdim_def('date',
                          units = 'days since 1970-01-01 00:00',
                          calendar = 'standard',
                          vals = dd)

    xx = ncdf4::nc_open(f, write = TRUE)

    for(newvar_name in c('Qma', 'QmaExtraNA')){

        newvar = get(newvar_name)
        newvar_nc <- ncvar_def(name = newvar_name,
                               units = 'cms', #if you change this, change ms_write_netcdf and other instances
                               dim = list(dim_time),
                               missval = NULL)
        xx2 = ncvar_add(xx, newvar_nc)
        ncvar_put(xx2, newvar_name, newvar)
    }

    for(newvar_name in c('QmaNAind', 'QmaExtraNAind')){

        newvar = get(newvar_name)
        newvar_nc <- ncvar_def(name = newvar_name,
                               units = '',
                               dim = list(dim_time),
                               missval = NULL)
        xx2 = ncvar_add(xx, newvar_nc)
        ncvar_put(xx2, newvar_name, newvar)
    }

    nc_close(xx)
}


### 6c CREATE NEW NETCDFs for experimentation (GAPPED), without re-running everything above ####

library(data.tree)

ncdfs = list.files('CAMELS_macrosheds_combined/time_series', full.names = TRUE)

segment_map = Node$new('gaps')
for(i in seq_len(nrow(neon_neighbors))){

    s = pull(neon_neighbors[i, 'nearest_neighbor'])
    f = grep(paste0('(?<!NHM_)', s, '\\.nc'), ncdfs, value = TRUE, perl = TRUE)
    f_newname = gsub('\\.nc$', '_GAPPED.nc', f)
    file.copy(f, f_newname, overwrite = TRUE)

    segment_map_site = segment_map$AddChild(s)

    xx = ncdf4::nc_open(f_newname, write = TRUE)
    qq = ncdf4::ncvar_get(xx, 'discharge')

    gap_lengths = c(64, 32, 16, 8, 4, 2)
    # if(sum(gap_lengths) > sum(! is.na(qq)) * 0.2) message(i)
    for(gl in gap_lengths){

        segment_map_g = segment_map_site$AddChild(gl)

        gap_placed = FALSE
        while(! gap_placed){

            ind = sample(gl:(length(qq) - gl - 1), 1)
            segment = ind:(ind + gl - 1)
            if(any(is.na(c(ind - 1,
                           qq[segment],
                           segment[length(segment)] + 1)))) next

            qq[segment] = NaN
            segment_map_g$AddChild(ind)
            gap_placed = TRUE
        }
    }

    ncvar_put(xx, 'discharge', qq)
    nc_close(xx)
}

segment_map = ToListSimple(segment_map)
saveRDS(segment_map, 'informative_stuff/introduced_testing_gaps.rds')

# gap_map = readRDS('informative_stuff/introduced_testing_gaps.rds')
# gap_map = as.Node(gap_map)

### 6d CREATE NEW NETCDFs for experimentation (neon manual Q as test Q), without re-running everything above ####

ncdfs = list.files('CAMELS_macrosheds_combined/time_series', full.names = TRUE)

#create _MANUALQ netcdfs (see next for back-extrapolatable versions of these
uniq_fls = c()
for(i in seq_along(neon_sites)){

    s = neon_sites[i]
    print(s)
    f = grep(paste0('(?<!NHM_)', s, '\\.nc'), ncdfs, value = TRUE, perl = TRUE)
    if(! length(f)) next
    f_newname = gsub('\\.nc$', '_MANUALQ.nc', f)
    file.copy(f, f_newname, overwrite = TRUE)

    # xx = ncdf4::nc_open(f, write = TRUE)
    # qq = ncdf4::ncvar_get(xx, 'discharge')
    # dd = ncdf4::ncvar_get(xx, 'date')
    # nc_close(xx)

    # xx = ncdf4::nc_open(f_newname)
    # qq = ncdf4::ncvar_get(xx, 'Qlag0')
    # xx = ncdf4::nc_open(f, write = TRUE)
    xx = ncdf4::nc_open(f_newname, write = TRUE)
    qq = ncdf4::ncvar_get(xx, 'discharge')
    # tt = ncdf4::ncvar_get(xx, 'tmax')

    dd = ncdf4::ncvar_get(xx, 'date')

    #field discharge measurements
    qd = neonUtilities::loadByProduct('DP1.20048.001', site = s, check.size = FALSE)
    uniq_fls = union(uniq_fls, grep('fieldData', names(qd), value = TRUE))

    q1 = q2 = tibble()
    try({
        q1 = select(qd$dsc_fieldDataADCP, discharge = totalDischarge, date = startDate, totalDischargeUnits)
    }, silent = TRUE)
    try({
        q2 = select(qd$dsc_fieldData, discharge = totalDischarge, date = startDate, totalDischargeUnits)
    }, silent = TRUE)
    if(nrow(q1) && nrow(q2)){
        q = bind_rows(q1, q2)
    } else if(nrow(q1)){
        q = q1
    } else {
        q = q2
    }

    if(any(! q$totalDischargeUnits %in% c('cubicMetersPerSecond', 'litersPerSecond'))) stop('!')

    q = q %>%
        mutate(discharge = ifelse(totalDischargeUnits == 'cubicMetersPerSecond', discharge * 1000, discharge),
               date = as.Date(date),
               site_code = s) %>%
        left_join(ms_areas, by = 'site_code') %>%
        filter(! is.na(ws_area_ha)) %>%
        #      mm/d        L/s         m^3/L   ha/m^2 mm/m   s/d     ha
        mutate(discharge = discharge * 0.001 * 1e-4 * 1000 * 86400 / ws_area_ha) %>%
        select(-ws_area_ha, -site_code, -totalDischargeUnits) %>%
        distinct(date, .keep_all = TRUE) %>%
        # mutate(fieldmeas = TRUE) %>%
        as_tibble()

    # #salt-based discharge measurements for reaeration experiments
    # qs = neonUtilities::loadByProduct('DP1.20193.001', site = s, check.size = FALSE)

    qq2 = left_join(tibble(date = dd), mutate(q, date = as.numeric(date)), by = 'date')
    qq2$discharge[is.na(qq2$discharge)] = NaN

    trycatch = try({
    # plot(as.Date(dd, origin = '1970-01-01'), qq, type = 'l', main = s,
    #      xlim=as.Date(c('2015-01-01', '2021-01-01')), ylim = c(0, 10))
    plot(as.Date(dd, origin = '1970-01-01'), qq, type = 'l', main = s)
    # plot(as.Date(dd, origin = '1970-01-01'), qq, type = 'l', main = s, xlim=as.Date(c('2014-01-01', '2022-01-01')))
    points(q$date, q$discharge, col = 'red')
    # points(qq2$date, qq2$discharge, col='blue', pch = 3)
    }, silent = TRUE)
    if(inherits(trycatch, 'try-error')) plot(q$date, q$discharge, col = 'red')

    # readLines(con = stdin(), 1)
    if(nrow(qq2) != length(qq)) stop('!!')

    ncvar_put(xx, 'discharge', qq2$discharge)
    nc_close(xx)
}

#create _MANQ_EXTRAP netcdfs, for predicting the full time range but evaluating on NEON field Q
uniq_fls = c()
for(i in seq_along(neon_sites)){

    s = neon_sites[i]
    print(s)
    f = grep(paste0(s, '_extrapolate\\.nc'), ncdfs, value = TRUE, perl = TRUE)
    if(! length(f)) next
    f_newname = gsub('_extrapolate\\.nc$', '_MANQ_EXTRAP.nc', f)
    file.copy(f, f_newname, overwrite = TRUE)

    xx = ncdf4::nc_open(f_newname, write = TRUE)
    qq = ncdf4::ncvar_get(xx, 'discharge')

    dd = ncdf4::ncvar_get(xx, 'date')

    #field discharge measurements
    qd = neonUtilities::loadByProduct('DP1.20048.001', site = s, check.size = FALSE)
    uniq_fls = union(uniq_fls, grep('fieldData', names(qd), value = TRUE))

    q1 = q2 = tibble()
    try({
        q1 = select(qd$dsc_fieldDataADCP, discharge = totalDischarge, date = startDate, totalDischargeUnits)
    }, silent = TRUE)
    try({
        q2 = select(qd$dsc_fieldData, discharge = totalDischarge, date = startDate, totalDischargeUnits)
    }, silent = TRUE)
    if(nrow(q1) && nrow(q2)){
        q = bind_rows(q1, q2)
    } else if(nrow(q1)){
        q = q1
    } else {
        q = q2
    }

    if(any(! q$totalDischargeUnits %in% c('cubicMetersPerSecond', 'litersPerSecond'))) stop('!')

    q = q %>%
        mutate(discharge = ifelse(totalDischargeUnits == 'cubicMetersPerSecond', discharge * 1000, discharge),
               date = as.Date(date),
               site_code = s) %>%
        left_join(ms_areas, by = 'site_code') %>%
        filter(! is.na(ws_area_ha)) %>%
        #      mm/d        L/s         m^3/L   ha/m^2 mm/m   s/d     ha
        mutate(discharge = discharge * 0.001 * 1e-4 * 1000 * 86400 / ws_area_ha) %>%
        select(-ws_area_ha, -site_code, -totalDischargeUnits) %>%
        distinct(date, .keep_all = TRUE) %>%
        # mutate(fieldmeas = TRUE) %>%
        as_tibble()

    # #salt-based discharge measurements for reaeration experiments
    # qs = neonUtilities::loadByProduct('DP1.20193.001', site = s, check.size = FALSE)

    qq2 = left_join(tibble(date = dd), mutate(q, date = as.numeric(date)), by = 'date')
    qq2$discharge[is.na(qq2$discharge)] = NaN

    trycatch = try({
    # plot(as.Date(dd, origin = '1970-01-01'), qq, type = 'l', main = s,
    #      xlim=as.Date(c('2015-01-01', '2021-01-01')), ylim = c(0, 10))
    plot(as.Date(dd, origin = '1970-01-01'), qq, type = 'l', main = s)
    # plot(as.Date(dd, origin = '1970-01-01'), qq, type = 'l', main = s, xlim=as.Date(c('2014-01-01', '2022-01-01')))
    points(q$date, q$discharge, col = 'red')
    # points(qq2$date, qq2$discharge, col='blue', pch = 3)
    }, silent = TRUE)
    if(inherits(trycatch, 'try-error')) plot(q$date, q$discharge, col = 'red')

    # readLines(con = stdin(), 1)
    if(nrow(qq2) != length(qq)) stop('!!')

    ncvar_put(xx, 'discharge', qq2$discharge)
    nc_close(xx)
}

### 7 write CSVs for static attributes of basins ####

sites_includedA <- sapply(ms_ts_splt, function(x) x$site_code[1])
sites_includedB <- sapply(nhm_ms_ts_splt, function(x) x$site_code[1])
sites_included <- union(sites_includedA, sites_includedB)

ms_all %>%
    filter(site_code %in% sites_included) %>%
    write_csv('CAMELS_macrosheds_combined/attributes/ms_attributes.csv')

write_csv(camels_attrs,
          'CAMELS_macrosheds_combined/attributes/camels_attributes.csv')

### 7b update static attributes and date ranges (jk, that happens in combined_ms_prep) of GAPPED sites ####

camels_attrs = read_csv('CAMELS_macrosheds_combined/attributes/camels_attributes.csv')
ms_attrs = read_csv('CAMELS_macrosheds_combined/attributes/ms_attributes.csv')

nghbs = neon_neighbors$nearest_neighbor

ms_attrs %>%
    filter(site_code %in% nghbs) %>%
    mutate(site_code = paste0(site_code, '_GAPPED')) %>%
    bind_rows(ms_attrs) %>%
    write_csv('CAMELS_macrosheds_combined/attributes/ms_attributes.csv')

camels_attrs %>%
    filter(site_code %in% nghbs) %>%
    mutate(site_code = paste0(site_code, '_GAPPED')) %>%
    bind_rows(camels_attrs) %>%
    write_csv('CAMELS_macrosheds_combined/attributes/camels_attributes.csv')


# ms_rngs = read_csv('../imputation/data/nh_methods/dateranges_ms.csv')
camels_rngs = read_csv('../imputation/data/nh_methods/dateranges_cmls.csv')

nghbs = neon_neighbors$nearest_neighbor

# ms_rngs %>%
#     filter(basin_id %in% nghbs) %>%
#     mutate(basin_id = paste0(basin_id, '_GAPPED')) %>%
#     bind_rows(ms_rngs) %>%
#     write_csv('../imputation/data/nh_methods/dateranges_ms.csv')

camels_rngs %>%
    filter(basin_id %in% nghbs) %>%
    mutate(basin_id = paste0(basin_id, '_GAPPED')) %>%
    bind_rows(camels_rngs) %>%
    write_csv('../imputation/data/nh_methods/dateranges_cmls.csv')

### 7c update static attributes and date ranges (jk, that happens in combined_ms_prep) of MANUALQ sites ####

ms_attrs = read_csv('CAMELS_macrosheds_combined/attributes/ms_attributes.csv')

ms_attrs %>%
    filter(site_code %in% neon_sites) %>%
    mutate(site_code = paste0(site_code, '_MANUALQ')) %>%
    bind_rows(ms_attrs) %>%
    write_csv('CAMELS_macrosheds_combined/attributes/ms_attributes.csv')

# ms_rngs = read_csv('../imputation/data/nh_methods/dateranges_ms.csv')
#
# ms_rngs %>%
#     filter(basin_id %in% neon_sites) %>%
#     mutate(basin_id = paste0(basin_id, '_MANUALQ')) %>%
#     bind_rows(ms_rngs) %>%
#     write_csv('../imputation/data/nh_methods/dateranges_ms.csv')

### 8 write list of MS sites and domains with netcdf files ####

s <- read_csv('informative_stuff/site_data.csv') %>%
    filter(site_type == 'stream_gauge') %>%
    select(domain, site_code)

sites_with_attrs <- read_csv('CAMELS_macrosheds_combined/attributes/ms_attributes.csv') %>%
    filter(! substr(site_code, 1, 3) == 'NHM') %>%
    pull(site_code)

ms_basins_with_netcdfs <-
    list.files('CAMELS_macrosheds_combined/time_series/') %>%
    str_subset('NHM', negate = TRUE) %>%
    str_sub(1, nchar(.) - 3) %>%
    str_subset('^[0-9]+', negate = TRUE) %>%
    tibble(site_code = .) %>%
    left_join(s, by = 'site_code') %>%
    select(domain, site_code) %>%
    filter(site_code %in% sites_with_attrs)

write_csv(ms_basins_with_netcdfs, 'informative_stuff/ms_basins_with_netcdfs.csv')

### --- end of normal sequence --- 9a for troubleshooting netcdf objects ####

library(dygraphs)
library(xts)

xx = ncdf4::nc_open('in/lstm_data/time_series/NHM_PH.nc')
xx = ncdf4::nc_open('/home/mike/git/macrosheds/qa_experimentation/data/CAMELS_macrosheds_combined/time_series/NHM_PH.nc')
xx = ncdf4::nc_open('CAMELS_macrosheds_combined/time_series/BIGC_MANUALQ.nc')
xx = ncdf4::nc_open('CAMELS_macrosheds_combined/time_series/TECR_MANQ_EXTRAP.nc')
xx = ncdf4::nc_open('CAMELS_macrosheds_combined/time_series/east_fork_GAPPED.nc')
xx = ncdf4::nc_open('CAMELS_macrosheds_combined/time_series/UNHC.nc')

# for(v in names(xx$var)){
#     vv = ncdf4::ncvar_get(xx, v)
#     cat(paste(v, '\t\t\t', round(max(vv, na.rm = TRUE)), '\n'))
# }
# for(v in names(xx$var)){
#     vv = ncdf4::ncvar_get(xx, v)
#     cat(paste(v, '\t\t\t', any(is.na(vv)), '\n'))
# }

qq = ncdf4::ncvar_get(xx, 'discharge')
qq[qq == -999] = NA
mean(qq, na.rm=T)
dd = ncdf4::ncvar_get(xx, 'date')
zz = bind_cols(dd,
               qq,
               ncdf4::ncvar_get(xx, 'prcp'),
               ncdf4::ncvar_get(xx, 'dayl'),
               ncdf4::ncvar_get(xx, 'srad'),
               ncdf4::ncvar_get(xx, 'swe'),
               ncdf4::ncvar_get(xx, 'tmax'),
               ncdf4::ncvar_get(xx, 'tmin'),
               ncdf4::ncvar_get(xx, 'vp'),
               # ncdf4::ncvar_get(xx, 'Qlag0'),
               ncdf4::ncvar_get(xx, 'pet'))
               # ncdf4::ncvar_get(xx, 'QrandomNA'))
colnames(zz) = c('date', 'discharge', 'prcp', 'dayl', 'srad', 'swe', 'tmax', 'tmin', 'vp', 'pet')
# colnames(zz) = c('date', 'discharge', 'prcp', 'dayl', 'srad', 'swe', 'tmax', 'tmin', 'vp', 'pet', 'QrandomNA')
zz = mutate(zz, date = as.Date(date))#origin = '1970-01-01'))
apply(zz, 2, function(x) sum(is.na(x)) / length(x))
# dygraphs::dygraph(xts::xts(select(zz, -date),
dygraphs::dygraph(xts::xts(select(zz, discharge),
            order.by = zz$date)) %>%
    dygraphs::dyRangeSelector()
ncdf4::nc_close(xx)
# plot(zz$date, zz$prcp, type='l', xlim = as.Date(c('2019-01-01', '2019-03-01')))
# plot(zz$date, zz$dayl, type='l', xlim = as.Date(c('2019-01-01', '2019-03-01')))
# plot(zz$date, zz$srad, type='l', xlim = as.Date(c('2019-01-01', '2019-03-01')))
# plot(zz$date, zz$swe, type='l', xlim = as.Date(c('2019-01-01', '2019-03-01')))
# plot(zz$date, zz$tmax, type='l', xlim = as.Date(c('2019-01-01', '2019-03-01')))
# plot(zz$date, zz$tmin, type='l', xlim = as.Date(c('2019-01-01', '2019-03-01')))
# plot(zz$date, zz$vp, type='l', xlim = as.Date(c('2019-01-01', '2019-03-01')))
# plot(zz$date, zz$discharge, type='l', xlim = as.Date(c('2019-01-01', '2019-03-01')))

### 9b for finding usgs basins that are included in kratzert set, and tallying -999s ####

# PRECEDED BY: grep -rc '\-999' . > chili

ss = read_lines('CAMELS/basin_dataset_public_v1p2/usgs_streamflow/chili')
ss = grep(':0$', ss, value = TRUE, invert = TRUE)
ss = str_match(ss, './([0-9]+)/([0-9]+)_streamflow_qc.txt:([0-9]+)$')[, 2:4] %>%
    as_tibble() %>%
    rename(group = 1,
           basin = 2,
           n = 3) %>%
    mutate(n = as.numeric(n))

sss = ss %>%
    filter(n > 400,
           basin %in% k_list) %>%
    print(n = 50)

### 9c for troubleshooting usgs discharge files (CAMELS) ####

library(dygraphs)
library(xts)

qq <- read.fwf(
    # 'CAMELS/basin_dataset_public_v1p2/usgs_streamflow/08/07376000_streamflow_qc.txt',
    'CAMELS/basin_dataset_public_v1p2/usgs_streamflow/03/02196000_streamflow_qc.txt',
    widths = c(9, 4, 3, 3, 9, 2),
    header = FALSE,
    strip.white = TRUE,
    colClasses = 'character',
    col.names = c('basin_id', 'year', 'month', 'day', 'Q', 'flag')
) %>%
    as_tibble() %>%
    mutate(date = ymd(paste(year, month, day)),
           Q = as.numeric(Q)) %>%
    select(date, discharge = Q)

dygraphs::dygraph(xts(x = select(qq, -date),
                      order.by = qq$date)) %>%
    dyRangeSelector()
options(max.print = 2000)
which(qq$discharge == -999)
options(max.print = 1000)
nrow(qq)
sum(qq$discharge == -999) / nrow(qq)

### 9d exhaustive version of 9c, using results from 9b ####

library(dygraphs)
library(xts)

defpar = par(mar = c(2,5,0,0), oma = c(0,0,0,0))
plot(1, 1, type = 'n', xlim = c(0, 100), ylim = c(0, nrow(sss)), ylab = '', yaxt = 'n')
abline(v = c(85, 92))
for(i in seq_len(nrow(sss))){

    sssr = sss[i, ]

    qq <- read.fwf(
        glue('CAMELS/basin_dataset_public_v1p2/usgs_streamflow/{g}/{b}_streamflow_qc.txt',
             g = sssr$group,
             b = sssr$basin),
        widths = c(9, 4, 3, 3, 9, 2),
        header = FALSE,
        strip.white = TRUE,
        colClasses = 'character',
        col.names = c('basin_id', 'year', 'month', 'day', 'Q', 'flag')
    ) %>%
        as_tibble() %>%
        mutate(date = ymd(paste(year, month, day)),
               Q = as.numeric(Q)) %>%
        select(date, discharge = Q)

    missing_vals = which(qq$discharge == -999)
    scaled_missing = scales::rescale(missing_vals, to = c(0, 100), from = c(1, nrow(qq)))
    points(scaled_missing, rep(i, length(scaled_missing)), pch = 15)
    text(x = -8, y = i, labels = sssr$basin, xpd = NA)
}
par(defpar)



### 9e for exploring the omitted CAMELS basins ####

#SEE kratzert_newman_basin_exploration.R

# k_list = read_lines('/home/mike/git/macrosheds/qa_experimentation/outside_resources/lstm_for_pub/data/basin_list.txt')
# k_list = union(k_list, c('06879650', '14158790', '09510200', '01187300'))
#
# full_list = list.files('CAMELS/basin_dataset_public_v1p2/usgs_streamflow/',
#                        recursive = TRUE) %>%
#     str_match('^[0-9]{2}/([0-9]+)') %>%
#     {.[, 2]}
#
# o_list <- setdiff(full_list, k_list)
#
# x = read_csv2('/home/mike/git/macrosheds/qa_experimentation/data/CAMELS/camels_attributes_v2.0/camels_topo.txt') %>%
#     select(gauge_id, a1 = area_gages2, a2 = area_geospa_fabric)
#
# x %>%
#     mutate(discrep = abs(a2 - a1) / mean(a2, a1)) %>%
#     filter(a1 > 2000 & a2 > 2000)

### 9f for comparing NHM, true, Rodeo, etc. ####

library(dygraphs)
library(xts)

load_ms_ncdf <- function(site_code){
    xx = ncdf4::nc_open(glue('CAMELS_macrosheds_combined/time_series/{z}.nc',
                             z = site_code))
    qq = ncdf4::ncvar_get(xx, 'discharge')
    qq[qq == -999] = NA
    mean(qq, na.rm=T)
    dd = ncdf4::ncvar_get(xx, 'date')
    zz = bind_cols(dd,
                   qq,
                   ncdf4::ncvar_get(xx, 'prcp'),
                   ncdf4::ncvar_get(xx, 'dayl'),
                   ncdf4::ncvar_get(xx, 'srad'),
                   ncdf4::ncvar_get(xx, 'swe'),
                   ncdf4::ncvar_get(xx, 'tmax'),
                   ncdf4::ncvar_get(xx, 'tmin'),
                   ncdf4::ncvar_get(xx, 'vp'))
    colnames(zz) = c('date', 'discharge', 'prcp', 'dayl', 'srad', 'swe', 'tmax', 'tmin', 'vp')
    zz = mutate(zz, date = as.Date(date))#origin = '1970-01-01'))
    return(zz)
}
dy_q_cmp <- function(x){

    dygraphs::dygraph(xts::xts(select(x, discharge),
                               order.by = x$date)) %>%
        dygraphs::dyRangeSelector()
}


#HERE: can't join because of something weird about date columns
zz_nhm = load_ms_ncdf('NHM_BLUE')
zz_true = load_ms_ncdf('BLUE')
zz = full_join(select(zz_nhm, date, discharge),
               select(zz_true, date, discharge),
               by = 'date', suffix = c('_NHM', '_MS'))
dy_q(zz_nhm)
dy_q(zz_true)



plot(zz_nhm$date, zz_nhm$discharge, type='l', xlim = as.Date(c('2012-01-01', '2021-09-01')))
lines(zz_true$date, zz_true$discharge, col = 'red')

### 9g do any CAMELS sites have SWE data? ####

daymet_files <- list.files('CAMELS/basin_dataset_public_v1p2/basin_mean_forcing/daymet',
                           full.names = TRUE,
                           recursive = TRUE)
daymet_files <- grep('011230\\*|208111310_forcing|344894205_forcing', #rm wonky files
                     daymet_files,
                     value = TRUE,
                     invert = TRUE)
for(i in seq_along(daymet_files)){
    x = read.fwf(
        daymet_files[i],
        skip = 4,
        widths = c(5, 3, 3, 100),
        header = FALSE,
        strip.white = TRUE,
        colClasses = 'numeric',
        col.names = c('year', 'month', 'day', 'hour', 'dayl', 'prcp', 'srad',
                      'swe', 'tmax', 'tmin', 'vp')
    ) %>%
        pull(swe)
    if(any(x != 0)) stop()
}

maurer_files <- list.files('CAMELS/basin_dataset_public_v1p2/basin_mean_forcing/maurer/',
                           full.names = TRUE,
                           recursive = TRUE)
maurer_files <- grep('011230\\*|208111310_forcing|344894205_forcing', #rm wonky files
                     daymet_files,
                     value = TRUE,
                     invert = TRUE)
for(i in seq_along(maurer_files)){
    x = read.fwf(
        maurer_files[i],
        skip = 4,
        widths = c(5, 3, 3, 100),
        header = FALSE,
        strip.white = TRUE,
        colClasses = 'numeric',
        col.names = c('year', 'month', 'day', 'hour', 'dayl', 'prcp', 'srad',
                      'swe', 'tmax', 'tmin', 'vp')
    ) %>%
        pull(swe)
    if(any(x != 0)) stop()
}

nldas_files <- list.files('CAMELS/basin_dataset_public_v1p2/basin_mean_forcing/nldas/',
                           full.names = TRUE,
                           recursive = TRUE)
nldas_files <- grep('011230\\*|208111310_forcing|344894205_forcing', #rm wonky files
                     daymet_files,
                     value = TRUE,
                     invert = TRUE)
for(i in seq_along(nldas_files)){
    x = read.fwf(
        nldas_files[i],
        skip = 4,
        widths = c(5, 3, 3, 100),
        header = FALSE,
        strip.white = TRUE,
        colClasses = 'numeric',
        col.names = c('year', 'month', 'day', 'hour', 'dayl', 'prcp', 'srad',
                      'swe', 'tmax', 'tmin', 'vp')
    ) %>%
        pull(swe)
    if(any(x != 0)) stop()
}

#wow. nope.

### 9h visually inspect all MS sites with netcdfs ####

get_response_1char <- function(msg,
                               possible_chars,
                               subsequent_prompt = FALSE){

    #msg: character. a message that will be used to prompt the user
    #possible_chars: character vector of acceptable single-character responses
    #subsequent prompt: not to be set directly. This is handled by
    #   get_response_mchar during recursion.

    if(subsequent_prompt){
        cat(paste('Please choose one of:',
                  paste(possible_chars,
                        collapse = ', '),
                  '\n> '))
    } else {
        cat(msg)
    }

    ch <- as.character(readLines(con = stdin(), 1))

    if(length(ch) == 1 && ch %in% possible_chars){
        return(ch)
    } else {
        get_response_1char(msg = msg,
                           possible_chars = possible_chars,
                           subsequent_prompt = TRUE)
    }
}

mbwn <- read_csv('/home/mike/git/macrosheds/qa_experimentation/data/informative_stuff/ms_basins_with_netcdfs.csv') %>%
    filter(! domain %in% c('santa_barbara', 'neon', 'bear'))
ms_basins_with_netcdfs <- as.list(mbwn$site_code)
names(ms_basins_with_netcdfs) <- mbwn$domain

wonky_sites = legit_sites = questionable_sites = needs_work_sites = short_sites = c()
for(i in seq_along(ms_basins_with_netcdfs)){

    # dmn = names(ms_basins_with_netcdfs)[i]
    sit = ms_basins_with_netcdfs[[i]]

    xx = ncdf4::nc_open(glue('CAMELS_macrosheds_combined/time_series/{s}.nc',
                             s = sit))
    qq = ncdf4::ncvar_get(xx, 'discharge')
    qq[qq == -999] = NA
    mean(qq, na.rm=T)
    dd = ncdf4::ncvar_get(xx, 'date')
    zz = suppressMessages(bind_cols(dd,
                   qq,
                   ncdf4::ncvar_get(xx, 'prcp'),
                   ncdf4::ncvar_get(xx, 'dayl'),
                   ncdf4::ncvar_get(xx, 'srad'),
                   ncdf4::ncvar_get(xx, 'swe'),
                   ncdf4::ncvar_get(xx, 'tmax'),
                   ncdf4::ncvar_get(xx, 'tmin'),
                   ncdf4::ncvar_get(xx, 'vp'),
                   ncdf4::ncvar_get(xx, 'pet')))
    colnames(zz) = c('date', 'discharge', 'prcp', 'dayl', 'srad', 'swe', 'tmax', 'tmin', 'vp', 'pet')
    zz = mutate(zz, date = as.Date(date))# = '1970-01-01'))
    # dygraphs::dygraph(xts::xts(select(zz, -date),
    plt = dygraphs::dygraph(xts::xts(select(zz, discharge),
                               order.by = zz$date)) %>%
        dygraphs::dyRangeSelector()
    # plot(zz$date, zz$discharge, type='l', xlim = as.Date(c('2019-01-01', '2019-03-01')))

    na_proportions <- apply(select(zz, -date),
                            MARGIN = 2,
                            FUN = function(x) sum(is.na(x)) / length(x))

    # if(any(na_proportions > 0.9)){
        print(na_proportions)
    # }
    readLines(con=stdin(), 1)

    print(sit)
    print(plt)
    #
    # resp <- get_response_1char('[y/g] if good, [n/b] if bad, [q] if questionable, [s] if short, [w] if needs work> ',
    #                            c('y', 'g', 'n', 'b', 'q', 's', 'w'))
    #
    # if(resp %in% c('y', 'g')){
    #     legit_sites[length(legit_sites) + 1] = sit
    # } else if(resp %in% c('n', 'b')){
    #     wonky_sites[length(wonky_sites) + 1] = sit
    # } else if(resp == 'q'){
    #     questionable_sites[length(questionable_sites) + 1] = sit
    # } else if(resp == 's'){
    #     short_sites[length(short_sites) + 1] = sit
    # } else if(resp == 'w'){
    #     needs_work_sites[length(needs_work_sites) + 1] = sit
    # }

    # print(q_status_df$q_status[q_status_df$site_code == sit])
    # resp = readline('[x] if short, [q] if questionable, anything else to continue> ')
    # if(resp == 'x'){
    #     q_status_df$q_status[q_status_df$site_code == sit] = 'short'
    # } else if(resp == 'q'){
    #     q_status_df$q_status[q_status_df$site_code == sit] = 'questionable'
    # }
}

q_status_df = tibble(site_code = legit_sites,
       q_status = 'okay') %>%
    bind_rows(tibble(site_code = wonky_sites,
                     q_status = 'fail')) %>%
    bind_rows(tibble(site_code = questionable_sites,
                     q_status = 'questionable')) %>%
    bind_rows(tibble(site_code = short_sites,
                     q_status = 'short')) %>%
    bind_rows(tibble(site_code = needs_work_sites,
                     q_status = 'needs work')) %>%
    left_join(mbwn) %>%
    select(domain, site_code, q_status) %>%
    arrange(domain)

print(q_status_df, n=100)

write_csv(q_status_df, 'informative_stuff/ms_q_status.csv')
# q_status_df <- read_csv('informative_stuff/ms_q_status.csv')
# q_status_df$q_status[q_status_df$site_code == 'C4'] = 'questionable'
# q_status_df$q_status[q_status_df$site_code %in% c('GSCC01', 'GSCC02', 'GSCC03', 'GSCC04', 'WS79', 'QS', 'GFVN')] = 'okay'

#check the list of wonky sites
filter(q_status_df, q_status != 'okay') %>% print(n=100)

#take a tour of sites by status

q_status_filtered = filter(q_status_df, q_status == 'short')
for(i in 1:nrow(q_status_filtered)){

    sit = q_status_filtered$site_code[i]

    xx = ncdf4::nc_open(glue('CAMELS_macrosheds_combined/time_series/{s}.nc',
                             s = sit))
    qq = ncdf4::ncvar_get(xx, 'discharge')
    qq[qq == -999] = NA
    mean(qq, na.rm=T)
    dd = ncdf4::ncvar_get(xx, 'date')
    zz = bind_cols(dd,
                   qq,
                   ncdf4::ncvar_get(xx, 'prcp'),
                   ncdf4::ncvar_get(xx, 'dayl'),
                   ncdf4::ncvar_get(xx, 'srad'),
                   ncdf4::ncvar_get(xx, 'swe'),
                   ncdf4::ncvar_get(xx, 'tmax'),
                   ncdf4::ncvar_get(xx, 'tmin'),
                   ncdf4::ncvar_get(xx, 'vp'))
    colnames(zz) = c('date', 'discharge', 'prcp', 'dayl', 'srad', 'swe', 'tmax', 'tmin', 'vp')
    zz = mutate(zz, date = as.Date(date))#origin = '1970-01-01'))
    # dygraphs::dygraph(xts::xts(select(zz, -date),
    plt = dygraphs::dygraph(xts::xts(select(zz, discharge),
                                     order.by = zz$date)) %>%
        dygraphs::dyRangeSelector()

    print(sit)
    print(plt)
    resp = readline('any key>')
}

### 9i inspect neon netcdfs ####

nrmlz <- function(x){
    rngx <- range(x, na.rm = TRUE)
    nrm <- (x - rngx[1]) / (rngx[2] - rngx[1])
    return(nrm)
}

mbwnn <- read_csv('/home/mike/git/macrosheds/qa_experimentation/data/informative_stuff/ms_basins_with_netcdfs.csv') %>%
    filter(domain == 'neon')

for(i in 1:nrow(mbwnn)){

    sit = mbwnn$site_code[i]

    # xx = ncdf4::nc_open(glue('CAMELS_macrosheds_combined/time_series/{s}.nc',
    xx = ncdf4::nc_open(glue('~/git/macrosheds/papers/q_sim/in/lstm_data/time_series/{s}.nc',
                             s = sit))
    qq = ncdf4::ncvar_get(xx, 'discharge')
    qq[qq == -999] = NA
    # print(paste(sum(is.na(qq))))
    mean(qq, na.rm=T)
    dd = ncdf4::ncvar_get(xx, 'date')
    zz = bind_cols(dd,
                   qq,
                   ncdf4::ncvar_get(xx, 'prcp'),
                   ncdf4::ncvar_get(xx, 'dayl'),
                   ncdf4::ncvar_get(xx, 'srad'),
                   ncdf4::ncvar_get(xx, 'swe'),
                   ncdf4::ncvar_get(xx, 'tmax'),
                   ncdf4::ncvar_get(xx, 'tmin'),
                   ncdf4::ncvar_get(xx, 'vp'))
    colnames(zz) = c('date', 'discharge', 'prcp', 'dayl', 'srad', 'swe', 'tmax', 'tmin', 'vp')
    zz = mutate(zz, date = as.Date(date))#origin = '1970-01-01'))
    # dygraphs::dygraph(xts::xts(select(zz, -date),
    # plt = dygraphs::dygraph(xts::xts(select(zz, discharge),
    nrmzz <- select(zz, discharge, prcp) %>% mutate(across(everything(), nrmlz))
    plt = dygraphs::dygraph(xts::xts(nrmzz,
                                     order.by = zz$date)) %>%
        dygraphs::dyRangeSelector()
    print(paste(sit, 'P and Q'))
    print(plt)
    resp = readline('any key>')

    plt = dygraphs::dygraph(xts::xts(select(zz, discharge),
                                     order.by = zz$date)) %>%
        dygraphs::dyRangeSelector()
    print(paste(sit, 'Q'))
    print(plt)
    resp = readline('any key>')
}
### 9i inspect NHM netcdfs ####

setwd('~/git/macrosheds/papers/q_sim/')

nrmlz <- function(x){
    rngx <- range(x, na.rm = TRUE)
    nrm <- (x - rngx[1]) / (rngx[2] - rngx[1])
    return(nrm)
}
qqq = list.files('in/lstm_data/time_series/', pattern = 'NHM_')
qqq = grep('^NHM_[0-9]+\\.nc$', qqq, value = TRUE, invert = TRUE)
for(q in qqq){
    print(q)
    xx = ncdf4::nc_open(paste0('in/lstm_data/time_series/', q))
    tryCatch(ncdf4::ncvar_get(xx, 'pet'), error = function(e) print('shit'))
}


### 9j inspect ms attributes ####

ddd = read_csv('CAMELS_macrosheds_combined/attributes/ms_attributes.csv')
apply(ddd, 2, function(x) sum(is.na(x)) / length(x))

### 9k evaluate artificial gaps from _GAPPED sites ####

reticulate::use_condaenv('nh2')

gap_map = readRDS('informative_stuff/introduced_testing_gaps.rds')

xr <- reticulate::import("xarray")
pd <- reticulate::import("pandas")
np <- reticulate::import("numpy")

# site = names(gap_map)[2]
for(site in names(gap_map)){

    if(site == 'name') next
    a = gap_map[[site]]
    print(site)

    # gaplen = 4
    for(gaplen in c(2, 4, 8, 16, 32, 64)){

        st = as.numeric(grep('name', names(a[[as.character(gaplen)]]), invert = TRUE, value = TRUE))
        gapinds_plus = seq(st - 10, st + gaplen + 9, 1)

        xx = ncdf4::nc_open(glue('CAMELS_macrosheds_combined/time_series/{s}.nc',
                                 s = paste0(site)))
        qq = ncdf4::ncvar_get(xx, 'discharge')
        qq[qq == -999] = NA
        dd = as.Date(ncdf4::ncvar_get(xx, 'date'))
        ncdf4::nc_close(xx)

        # xx = reticulate::py_load_object('../imputation/src/nh_methods/runs/run1360_1907_011320/test/model_epoch030/test_results.p')
        xx = reticulate::py_load_object('../imputation/src/nh_methods/runs/run1361_1907_032422/test/model_epoch030/test_results.p')

        pred = xx[[paste0(site, '_GAPPED')]]$`1D`$xr$discharge_sim$to_pandas()
        pred = tibble(date = as.Date(rownames(pred)), Q = pred$`0`)

        ylm = try({range(c(qq[gapinds_plus], pred$Q[gapinds_plus]))})
        if(inherits(ylm, 'try-error') | all(is.na(ylm))) next
        plot(dd[gapinds_plus], qq[gapinds_plus], type = 'l', ylim = ylm)
        abline(v=c(dd[st], dd[st + gaplen - 1]), col = 'gray')
        lines(pred, col = 'orange')

        readLines(con=stdin(), 1)
    }
}

### 9l non-graphically explore all netcdfs ####

for(i in seq_along(ms_basins_with_netcdfs$site_code)){

    sit = ms_basins_with_netcdfs$site_code[i]
    xx = ncdf4::nc_open(glue('CAMELS_macrosheds_combined/time_series/{s}.nc',
                             s = sit))
    qq = ncdf4::ncvar_get(xx, 'discharge')
    if(any(is.na(qq) & ! is.nan(qq))) stop()
    ncdf4::nc_close(xx)
}

### 9m identify best runs whose IDs have been forgotten ####
gg = read_csv('../imputation/src/accumulated_results.csv')
filter(gg, round(fine_NSE, 2) == 0.75, site_code == 'BLWA')
filter(gg, site_code == 'MCRA') %>% select(experiment_name, base_NSE, fine_NSE) %>% print(n=500)
...

### 9n insert pet into all macrosheds NHM netcdfs ####

fs <- list.files('CAMELS_macrosheds_combined/time_series', pattern = '^NHM_')
nhm_ms_sites <- grep('^NHM_[0-9]+\\.nc$', fs, value = TRUE, invert = TRUE)
ms_sites_ <- sub('NHM_', '', nhm_ms_sites)
# ms_sites_ <- sub('\\.nc', '_extrapolate.nc', ms_sites_)

for(i in seq_along(nhm_ms_sites)){

    f_nhm <- nhm_ms_sites[i]
    f <- ms_sites_[i]

    xx = try({
        ncdf4::nc_open(glue('CAMELS_macrosheds_combined/time_series/{f}'))
    }, silent = TRUE)

    if(inherits(xx, 'try-error')){
        print(f)
        next
    }

    pet = tibble(date = ncdf4::ncvar_get(xx, 'date'),
                 pet = ncdf4::ncvar_get(xx, 'pet'))

    nc_close(xx)

    xx2 = ncdf4::nc_open(glue('CAMELS_macrosheds_combined/time_series/{f_nhm}'), write = TRUE)
    # ncdf4::ncvar_get(xx2, 'pet')

    dd = ncdf4::ncvar_get(xx2, 'date')
    new_ <- tibble(date = dd) %>%
        left_join(pet)

    dim_time <- ncdim_def('date',
                          units = 'days since 1970-01-01 00:00',
                          calendar = 'standard',
                          vals = dd)

    newvar_nc <- ncvar_def(name = 'pet',
                           units = 'mm',
                           dim = list(dim_time),
                           missval = NULL)

    # xx2_ = ncvar_add(xx2, newvar_nc)
    # ncvar_put(xx2_, 'pet', new_$pet)
    zz <- ncvar_get(xx2, 'pet')
    # ncvar_change_missval(xx2, 'pet', NA_real_)
    # zz[is.na(zz)] = -999
    ncvar_put(xx2, 'pet', new_$pet)

    nc_close(xx2)
}

# xx=ncdf4::nc_open(glue('CAMELS_macrosheds_combined/time_series/FLNT_extrapolate.nc'))
# pet = tibble(date = as.Date(ncdf4::ncvar_get(xx, 'date')),
#              pet = ncdf4::ncvar_get(xx, 'pet'))
# range(as.Date(pet$date))
# tail(pet)
# nc_close(xx)

### 9o create ms_basins_with_netcdfs again ####

sites_with_attrs <- read_csv('in/lstm_data/attributes/ms_attributes.csv') %>%
    filter(! substr(site_code, 1, 3) == 'NHM') %>%
    pull(site_code)

ms_basins_with_netcdfs <- list.files('in/lstm_data/time_series/') %>%
    str_subset('NHM', negate = TRUE) %>%
    str_sub(1, nchar(.) - 3) %>%
    str_subset('^[0-9]+', negate = TRUE) %>%
    tibble(site_code = .) %>%
    left_join(ms_site_data, by = 'site_code') %>%
    select(domain, site_code) %>%
    filter(site_code %in% sites_with_attrs)

write_csv(ms_basins_with_netcdfs, '~/git/macrosheds/qa_experimentation/data/informative_stuff/ms_basins_with_netcdfs.csv')

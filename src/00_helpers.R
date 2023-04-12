#helper functions for the other R scripts. No need to source these directly.
#start with 01_linear_regression.R

#general (regression)

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

# setup (regression)

build_dir_structure <- function(){

    dir.create('out', showWarnings = FALSE)
    dir.create('figs', showWarnings = FALSE)

    dir.create('in/usgs_Q', showWarnings = FALSE)

    dir.create('figs/lm_plots', showWarnings = FALSE)
    dir.create('figs/lm_plots/diag', showWarnings = FALSE)
    dir.create('figs/lm_plots/pred', showWarnings = FALSE)
    dir.create('figs/lm_plots/fit', showWarnings = FALSE)
    dir.create('figs/lm_plots/val', showWarnings = FALSE)

    dir.create('out/lm_out', showWarnings = FALSE)
    dir.create('out/lm_out/predictions', showWarnings = FALSE)
    dir.create('out/lm_out/fit', showWarnings = FALSE)
    dir.create('out/lm_out/summary', showWarnings = FALSE)
}

rename_dir_structure <- function(){

    file.rename('figs/lm_plots', 'figs/lm_plots_specQ')
    file.rename('out/lm_out', 'out/lm_out_specQ')
}

# retrieval

get_neon_field_discharge <- function(neon_sites){

    dir.create('in/NEON', showWarnings = FALSE)
    dir.create('in/neon_field_Q', showWarnings = FALSE)

    for(i in seq_along(neon_sites)){

        s <- neon_sites[i]
        print(s)

        #field discharge measurements
        # ff <- paste0('field_q_rds_files/', s, '.rds')
        if(file.exists(ff)){
            qd <- readRDS(ff)
        } else {
            qd <- neonUtilities::loadByProduct('DP1.20048.001', site = s, check.size = FALSE)
            # saveRDS(qd, paste0('field_q_rds_files/', s, '.rds'))
        }

        q1 <- q2 <- tibble()
        try({

            flag1 <- qd$dsc_fieldDataADCP[['totalDischargeCalcQF']]
            flag2 <- qd$dsc_fieldDataADCP[['dataQF']]
            if('finalDischarge' %in% colnames(qd$dsc_fieldDataADCP)) stop()
            q1 <- select(qd$dsc_fieldDataADCP, discharge1 = totalDischarge,
                         date = startDate, totalDischargeUnits) %>%
                mutate(discharge2 = NA_real_)
            q1$flag1 <- if(! is.null(flag1)) as.numeric(flag1) else rep(NA_real_, nrow(q1))
            q1$flag2 <- if(! is.null(flag2)) as.numeric(flag2) else rep(NA_real_, nrow(q1))

        }, silent = TRUE)

        try({

            flag1 <- qd$dsc_fieldData[['totalDischargeCalcQF']]
            flag2 <- qd$dsc_fieldData[['dataQF']]
            q2 <- select(qd$dsc_fieldData, discharge1 = totalDischarge, discharge2 = finalDischarge,
                         date = startDate, totalDischargeUnits)
            q2$flag1 <- if(! is.null(flag1)) as.numeric(flag2) else rep(NA_real_, nrow(q2))
            q2$flag2 <- if(! is.null(flag2)) as.numeric(flag2) else rep(NA_real_, nrow(q2))

        }, silent = TRUE)

        if(nrow(q1) && nrow(q2)){
            q <- bind_rows(q1, q2)
        } else if(nrow(q1)){
            q <- q1
        } else {
            q <- q2
        }

        if(any(! q$totalDischargeUnits %in% c('cubicMetersPerSecond', 'litersPerSecond'))) stop('new unit detected. account for this.')

        q <- mutate(q, discharge1 = ifelse(totalDischargeUnits == 'cubicMetersPerSecond', discharge1 * 1000, discharge1))
        # q[is.na(q$discharge1) + is.na(q$discharge2) == 1, ]
        q$discharge <- q$discharge2
        missing_or_bogus <- is.na(q$discharge) | q$discharge < 0
        q$discharge[missing_or_bogus] <- q$discharge1[missing_or_bogus]

        print(paste('removing', sum(! is.na(q$flag2) & q$flag2 == 1), 'rows'))

        q <- q %>%
            filter(discharge >= 0,
                   is.na(flag2) | flag2 <= 0) %>%
            select(datetime = date, discharge) %>%
            mutate(site_code = s) %>%
            as_tibble()

        write_csv(q, glue('in/neon_field_Q/{s}.csv'))
    }
}

get_neon_inst_discharge <- function(neon_sites){

    dir.create('in/NEON', showWarnings = FALSE)
    dir.create('in/NEON/neon_continuous_Q')

    for(i in seq_along(neon_sites)){

        s <- neon_sites[i]
        # print(s)

        #continuous discharge measurements
        # if(file.exists(paste0('debris/neon_q_dl/', s, '.rds'))) next #

        qd <- neonUtilities::loadByProduct('DP4.00130.001', site = s,
                                           check.size = FALSE,
                                           release = 'RELEASE-2023')

        # write_lines(paste0(s, '\n'), 'log/neonq_qc.txt', append = T) #
        # saveRDS(qd, paste0('debris/neon_q_dl', s, '.rds')) #

        q1 <- q2 <- tibble()

        if('csd_continuousDischarge' %in% names(qd)){

            q1 <- qd$csd_continuousDischarge

            # xx = table(q1$dischargeFinalQF, useNA = 'always')
            # write_lines('FinalQF', 'log/neonq_qc.txt', append = T)
            # write_lines(paste(names(xx), collapse = ', '), 'log/neonq_qc.txt', append = T)
            # write_lines(paste(xx, collapse = ', '), 'log/neonq_qc.txt', append = T)
            # xx = table(q1$dischargeFinalQFSciRvw, useNA = 'always')
            # write_lines('FinalQFSciRvm', 'log/neonq_qc.txt', append = T)
            # write_lines(paste(names(xx), collapse = ', '), 'log/neonq_qc.txt', append = T)
            # write_lines(paste(xx, collapse = ', '), 'log/neonq_qc.txt', append = T)
            # write_lines('', 'log/neonq_qc.txt', append = T)

            q1 <- q1 %>%
                filter(is.na(dischargeFinalQF) | dischargeFinalQF == 0,
                       is.na(dischargeFinalQFSciRvw) | dischargeFinalQFSciRvw == 0) %>%
                select(discharge = maxpostDischarge,
                       datetime = endDate, discharge_lower = withRemnUncQLower2Std,
                       discharge_upper = withRemnUncQUpper2Std)
        }

        if('csd_continuousDischargeUSGS' %in% names(qd)){

            q2 <- qd$csd_continuousDischargeUSGS

            # xx = table(q2$dischargeFinalQF, useNA = 'always')
            # write_lines('FinalQF', 'log/neonq_qc.txt', append = T)
            # write_lines(paste(names(xx), collapse = ', '), 'log/neonq_qc.txt', append = T)
            # write_lines(paste(xx, collapse = ', '), 'log/neonq_qc.txt', append = T)
            # xx = table(q2$dischargeFinalQFSciRvw, useNA = 'always')
            # write_lines('FinalQFSciRvm', 'log/neonq_qc.txt', append = T)
            # write_lines(paste(names(xx), collapse = ', '), 'log/neonq_qc.txt', append = T)
            # write_lines(paste(xx, collapse = ', '), 'log/neonq_qc.txt', append = T)
            # write_lines('', 'log/neonq_qc.txt', append = T)

            q2 <- q2 %>%
                filter(is.na(dischargeFinalQFSciRvw) | dischargeFinalQFSciRvw == 0) %>%
                select(discharge = usgsDischarge,
                       datetime = endDate, discharge_lower = withRegressionUncQLower2Std,
                       discharge_upper = withRegressionUncQUpper2Std)
        }

        if(nrow(q1) && nrow(q2)){
            q <- bind_rows(q1, q2)
        } else if(nrow(q1)){
            q <- q1
        } else {
            q <- q2
        }

        # write_lines('---------------\n\n', 'log/neonq_qc.txt', append = T)

        q <- as_tibble(q) %>%
            mutate(site_code = s)

        write_csv(q, glue('in/NEON/neon_continuous_Q/{s}.csv'))
    }
}

# data prep

assemble_q_df <- function(neon_site, nearby_usgs_gages = NULL, ms_Q_data = NULL,
                          datetime_snapdist_hrs = 12, overwrite = TRUE,
                          scale_q_by_area = TRUE){

    #ms_Q_data: data.frame with posixct datetime column, and Q columns. Each Q column must
    #have its site_code as header. transformed columns will be created

    neon_q_manual = read_csv(glue('in/neon_field_Q/{neon_site}.csv')) %>%
        mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
        rename(discharge_manual = discharge) %>%
        distinct(datetime, .keep_all = TRUE)
    earliest_date = as.character(date(min(neon_q_manual$datetime)))

    if(scale_q_by_area){
        wsa = filter(ms_areas, site_code == neon_site) %>% pull(ws_area_ha)
        neon_q_manual$discharge_manual = neon_q_manual$discharge_manual / wsa * 1000
    }

    if(! file.exists(glue('in/usgs_Q/{neon_site}.csv')) || overwrite){

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

        write_csv(site_nearby, glue('in/usgs_Q/{neon_site}.csv'))

    } else {
        site_nearby = read_csv(glue('in/usgs_Q/{neon_site}.csv'))
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

assemble_q_df_daily <- function(neon_site, ms_Q_data = NULL, scale_q_by_area = TRUE,
                                overwrite = TRUE){

    #ms_Q_data: data.frame with date column, and Q columns. Each Q column must
    #have its site_code as header. transformed columns will be created

    neon_q_manual = read_csv(glue('in/neon_field_Q/{neon_site}.csv')) %>%
        mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
        rename(discharge_manual = discharge) %>%
        distinct(datetime, .keep_all = TRUE) %>%
        mutate(date = as.Date(datetime)) %>%
        select(-datetime)

    if(scale_q_by_area){
        wsa = filter(ms_areas, site_code == neon_site) %>% pull(ws_area_ha)
        neon_q_manual$discharge_manual = neon_q_manual$discharge_manual / wsa * 1000
    }

    if(! file.exists(glue('in/usgs_Q/{neon_site}.csv')) || overwrite){

        ms_sites = grep('date', colnames(ms_Q_data), invert = TRUE, value = TRUE)

        if(scale_q_by_area){

            for(g in ms_sites){

                wsa = filter(ms_areas, site_code == g) %>% pull(ws_area_ha)
                ms_Q_data[[g]] = ms_Q_data[[g]] / wsa * 1000
            }
        }

        ms_Q_data = ms_Q_data %>%
            mutate(across(-date, neglog, .names = '{.col}_log'))

        write_csv(ms_Q_data, glue('in/usgs_Q/{neon_site}.csv'))

    } else {
        ms_Q_data = read_csv(glue('in/usgs_Q/{neon_site}.csv'))
    }

    joined = left_join(neon_q_manual, ms_Q_data, by = 'date') %>%
        arrange(date) %>%
        rename(discharge = discharge_manual) %>%
        mutate(season = factor(lubridate::quarter(date)),
               discharge_log = neglog(discharge)) %>%
        select(site_code, date, discharge, discharge_log, everything())

    return(joined)
}

# linear regression

generate_nested_formulae <- function(full_spec, d,
                                     interactions = TRUE, through_origin = TRUE,
                                     min_points_per_param = 15,
                                     max_interaction = Inf){

    #full_spec: a formula representing the full model specification (all possible predictors)
    #d: a data frame. must include dependent variable from full_spec as a column name
    #interactions: are interaction terms allowed? see max_interaction
    #through_origin: should 0-intercept models be generated?
    #min_points_per_param: at least this many data points must be present for each
    #   parameter included in resulting formulae.
    #max_interaction: the maximum number of terms that may interact, e.g. if 3,
    #   then formulae with x:y:z:w will be filtered from the results.

    #returns a list of formulae

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

    if(! length(indeps)) stop('not enough data')

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

    mods <- Filter(function(x){
        trms <- attr(terms(x), 'term.labels')
        length(trms) > 1 || trms != 'season'
    }, mods)

    return(mods)
}

eval_model_set <- function(data, model_list, metric,
                           unscale_q_by_area = TRUE, k_folds = 10){

    #data: a data frame with a date column, and all columns referenced in model_list
    #model_list: a list of lm formulae, as generated by generate_nested_formulae
    #metric: [IGNORED - uses NSE] a function for evaluating the models. the first argument to this function
    #   will be the predictions, and the second will be the observations for the dependent variable

    log = 'xy'
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
            # ndata_per_fold = sum(! is.na(dd[[dep]])) %/% 30
            # ndata_per_fold = ifelse(ndata_per_fold == 0, 1, ndata_per_fold)

            tryCatch({
                dd = CVlm(data = dd, model_list[[i]], m = k_folds,
                          plotit = FALSE, printit = FALSE)
            }, warning = function(w) stop('rank-deficiency detected. moving on.'))

            if(log == 'xy'){

                neon_site = data$site_code[1]
                if(unscale_q_by_area){
                    wsa = filter(ms_areas, site_code == neon_site) %>% pull(ws_area_ha)
                    sim = inv_neglog(dd$cvpred) * wsa / 1000
                    nsenum = sum((inv_neglog(dd$cvpred) * wsa / 1000 - inv_neglog(dd[[dep]]) * wsa / 1000)^2) #just checking
                } else {
                    sim = inv_neglog(dd$cvpred)
                    nsenum = sum((inv_neglog(dd$cvpred) - inv_neglog(dd[[dep]]))^2)
                }

            # } else if(log == 'y'){
            #     stop('do not use log = "y"')
            #     nsenum = sum((inv_neglog(dd$cvpred) - dd[[dep]])^2)
            # } else if(log == FALSE){
            #     stop('do not use log = FALSE')
            #     nsenum = sum((dd$cvpred - dd[[dep]])^2)
            }

            if(unscale_q_by_area){
                obs = pull(dd[dep_transformed]) * wsa / 1000
            } else {
                obs = pull(dd[dep_transformed])
            }

            nsedenom = sum((obs - mean(obs))^2)
            nse_ = 1 - (nsenum / nsedenom)

            nse = hydroGOF::NSE(sim, obs)

            if(! all.equal(nse, nse_)) stop('oi')

            if(nse > best_score){
                best_score = nse
                best_mod = i

                metr_cv = tibble(
                    kge = hydroGOF::KGE(sim, obs),
                    rmse = hydroGOF::rmse(sim, obs),
                    mae = hydroGOF::mae(sim, obs),
                    pbias = hydroGOF::pbias(sim, obs)
                )
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

    nse_out = hydroGOF::NSE(data$lm, data$discharge)

    metr = tibble(
        kge = hydroGOF::KGE(data$lm, data$discharge),
        rmse = hydroGOF::rmse(data$lm, data$discharge),
        mae = hydroGOF::mae(data$lm, data$discharge),
        pbias = hydroGOF::pbias(data$lm, data$discharge)
    )

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
    #     # xx = reticulate::py_load_object('src/nh_methods/runs/run1360_1907_011320/test/model_epoch030/test_results.p')
    #     xx = reticulate::py_load_object('src/nh_methods/runs/run1361_1907_032422/test/model_epoch030/test_results.p')
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
               other_metrics = metr,
               other_metrics_crossval = metr_cv,
               # plot = dg,
               fits = gd)

    return(out)
}

plots_and_results <- function(neon_site, best, lm_df, results, return_plot = FALSE,
                              unscale_q_by_area = TRUE){

    if(length(best$prediction) != nrow(lm_df)) stop('oi')

    #load corroborating usgs/ms site data

    sites_nearby = read_csv(glue('in/usgs_Q/{neon_site}.csv'))

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

    neon_q_auto = read_csv(glue('in/NEON/neon_continuous_Q/{neon_site}.csv')) %>%
        filter(! is.na(discharge)) %>%
        rename(discharge_auto = discharge)

    q_eval = read_csv('in/neon_q_eval.csv') %>%
        filter(site == neon_site)

    q_eval = q_eval %>%
        group_by(site, year, month) %>%
        summarize(keep = all(final_qual %in% c('Tier1', 'Tier2')), #if issues occur at any point in a month, reject that whole month
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
               # Q_neon_continuous_filtered = discharge_auto_qc,
               Q_neon_continuous = discharge_auto) %>%
        filter(if_any(-datetime, ~cumsum(! is.na(.)) != 0)) %>%  #remove leading NA rows
        arrange(datetime)

    dg = dygraphs::dygraph(xts(x = select(out_data, -datetime, -Q_pred_int_2.5, -Q_pred_int_97.5) %>% tail(5e5),
                               order.by = tail(out_data$datetime, 5e5))) %>%
        dyRangeSelector()

    saveWidget(dg, glue('figs/lm_plots/pred/{neon_site}_log.html'))

    if(inherits(best$fits, 'grob')){

        #make diagnostic and fit plots
        png(glue('figs/lm_plots/diag/{neon_site}_diag_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
        defpar = par(mfrow=c(2,2))
        plot(best$best_model_object)
        par(defpar)
        dev.off()

        png(glue('figs/lm_plots/fit/{neon_site}_fit_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
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

    png(glue('figs/lm_plots/val/{neon_site}_obs_v_pred.png'), 6, 6, 'in', type = 'cairo', res = 300)
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
    write_csv(out_data, glue('out/lm_out/predictions/{neon_site}.csv'))

    #save fit data as CSV
    write_csv(best$lm_data, glue('out/lm_out/fit/{neon_site}.csv'))

    #save model summary
    sink(glue('out/lm_out/summary/{neon_site}.txt'))
    print(summary(best$best_model_object))
    sink()

    #return results frame, updated with this site's results
    results$bestmod[results$site_code == neon_site] = as.character(best$best_model)[3]
    results$nse[results$site_code == neon_site] = best$score
    results$nse_cv[results$site_code == neon_site] = best$score_crossval
    results$kge[results$site_code == neon_site] = best$other_metrics$kge
    results$kge_cv[results$site_code == neon_site] = best$other_metrics_crossval$kge
    results$pbias[results$site_code == neon_site] = best$other_metrics$pbias
    results$pbias_cv[results$site_code == neon_site] = best$other_metrics_crossval$pbias
    results$adj_r_squared[results$site_code == neon_site] = summary(best$best_model_object)$adj.r.squared

    if(unscale_q_by_area){
        write_csv(results, 'out/lm_out/results_specificq.csv')
    } else {
        write_csv(results, 'out/lm_out/results.csv')
    }

    if(return_plot){
        return(list(plot = dg, results = results))
    } else {
        return(results)
    }
}

plots_and_results_daily_composite <- function(neon_site, best1, best2, lm_df1,
                                              lm_df2, results,
                                              unscale_q_by_area = TRUE){

    #best2 and lm_df2 should represent the model with more terms included

    if(length(best1$prediction) != nrow(lm_df1)) stop('oi')
    if(length(best2$prediction) != nrow(lm_df2)) stop('oi')

    #load corroborating usgs/ms site data

    sites_nearby = read_csv(glue('in/usgs_Q/{neon_site}.csv'))

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

    neon_q_daily = read_csv(glue('in/NEON/neon_continuous_Q/{neon_site}.csv')) %>%
        group_by(date = date(datetime)) %>%
        summarize(discharge = mean(discharge, na.rm = TRUE)) %>%
        filter(! is.na(discharge)) %>%
        rename(discharge_daily = discharge)

    q_eval = read_csv('in/neon_q_eval.csv') %>%
        filter(site == neon_site)

    q_eval = q_eval %>%
        group_by(site, year, month) %>%
        summarize(keep = all(final_qual %in% c('Tier1', 'Tier2')), #if issues occur at any point in a month, reject that whole month
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

    saveWidget(dg, glue('figs/lm_plots/pred/{neon_site}_log.html'))

    #make diagnostic and fit plots
    png(glue('figs/lm_plots/diag/{neon_site}_diag_log_modA.png'), 6, 6, 'in', type = 'cairo', res = 300)
    defpar = par(mfrow=c(2,2))
    plot(best1$best_model_object)
    par(defpar)
    dev.off()

    png(glue('figs/lm_plots/diag/{neon_site}_diag_log_modB.png'), 6, 6, 'in', type = 'cairo', res = 300)
    defpar = par(mfrow=c(2,2))
    plot(best2$best_model_object)
    par(defpar)
    dev.off()

    png(glue('figs/lm_plots/fit/{neon_site}_fit_log_modA.png'), 6, 6, 'in', type = 'cairo', res = 300)
    plot(best1$fits)
    dev.off()

    png(glue('figs/lm_plots/fit/{neon_site}_fit_log_modB.png'), 6, 6, 'in', type = 'cairo', res = 300)
    plot(best2$fits)
    dev.off()

    #plot predictions versus field measurements. need to round field meas to Q interval
    nse = hydroGOF::NSE(out_data$Q_predicted, out_data$Q_neon_field)
    kge = hydroGOF::KGE(out_data$Q_predicted, out_data$Q_neon_field)
    pbias = hydroGOF::pbias(out_data$Q_predicted, out_data$Q_neon_field)
    axlim = c(0, max(c(out_data$Q_predicted, out_data$Q_neon_field), na.rm = TRUE))

    png(glue('figs/lm_plots/val/{neon_site}_obs_v_pred.png'), 6, 6, 'in', type = 'cairo', res = 300)
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
    write_csv(out_data, glue('out/lm_out/predictions/{neon_site}.csv'))

    #save fit data as CSV
    write_csv(best2$lm_data, glue('out/lm_out/fit/{neon_site}.csv'))

    #save model summary
    sink(glue('out/lm_out/summary/{neon_site}.txt'))
    print(summary(best1$best_model_object))
    print('---')
    print(summary(best2$best_model_object))
    sink()

    #return results frame, updated with this site's results
    results$bestmod[results$site_code == neon_site] = paste(as.character(best1$best_model)[3],
                                                            as.character(best1$best_model)[3],
                                                            sep = ' & ')
    results$nse[results$site_code == neon_site] = nse
    results$nse_cv[results$site_code == neon_site] = mean(best1$score_crossval, best2$score_crossval)
    results$adj_r_squared[results$site_code == neon_site] = mean(summary(best1$best_model_object)$adj.r.squared,
                                                                 summary(best2$best_model_object)$adj.r.squared)

    results$kge[results$site_code == neon_site] = kge
    results$kge_cv[results$site_code == neon_site] = mean(best1$other_metrics$kge, best2$other_metrics$kge)
    results$pbias[results$site_code == neon_site] = pbias
    results$pbias_cv[results$site_code == neon_site] = mean(best1$other_metrics$pbias, best2$other_metrics$pbias)

    if(unscale_q_by_area){
        write_csv(results, 'out/lm_out/results_specificq.csv')
    } else {
        write_csv(results, 'out/lm_out/results.csv')
    }

    return(results)
}

# random forest regression

eval_model_set_rf <- function(data){

    #data: a data frame with a date column, and all columns referenced in model_list

    # getModelInfo('ranger')
    #establish param search grid
    maxtry <- max(1, floor(sqrt(ncol(data) / 2 - 1)))
    tunegrid <- expand.grid(mtry = 1:maxtry,
                            splitrule = c('variance', 'extratrees', 'maxstat'),
                            min.node.size = 1:10)

    #train
    tc <- trainControl(method = 'cv', number = 10, search = 'grid')
    model <- train(as.formula(glue(formulaA)),
                   data = data[complete.cases(data), ],
                   trControl = tc, method = 'ranger', metric = 'RMSE', #using RMSE because only care about prediction accuracy, not variance explained
                   tuneGrid = tunegrid)
                   # respect.unordered.factors = 'order')
                   # importance = 'permutation')

    results <- filter(model$results, RMSE == min(RMSE))

    #predict on train set
    data$rf <- NA_real_
    data$rf[complete.cases(select(data, -rf))] <- ranger::predictions(model$finalModel)
    data$rf[data$rf < 0] = 0

    results$NSE <- hydroGOF::NSE(data$rf, data$discharge)
    results$KGE <- hydroGOF::KGE(data$rf, data$discharge)

    #clean and return
    first_non_na = Position(function(x) ! is.na(x), data$rf)
    last_non_na = nrow(data) - Position(function(x) ! is.na(x), rev(data$rf)) + 1
    plot_data = data[first_non_na:last_non_na, ]

    plot_data = select(plot_data, -ends_with('_log')) %>%
        select(site_code, any_of(c('date', 'datetime')), Q_neon_field = discharge, Q_predicted = rf,
               everything())

    out = list(best_model_object = model$finalModel,
               prediction = unname(data$rf),
               rf_data = plot_data,
               results = results,
               # importance = ranger::importance(model$finalModel),
               search_plot = ggplot(model))

    return(out)
}

plots_and_results_rf <- function(neon_site, best, rf_df, results,
                                 return_plot = FALSE){

    if(length(best$prediction) != nrow(rf_df)) stop('oi')

    #load corroborating usgs/ms site data

    sites_nearby = read_csv(glue('in/usgs_Q/{neon_site}.csv')) %>%
        select(-ends_with('_log')) %>%
        rename_with(~paste0('x', .), matches('^[0-9]+$'))

    #filter usgs/ms site data that didn't end up in the model
    # trms = rownames(attributes(terms(best$best_model))$factors)
    # dep = gsub('`', '', trms[1])
    # indeps = gsub('`', '', trms[-1])
    # indeps = gsub('_log', '', indeps)
    # if(length(indeps) == 1 && indeps == 'season') stop('season is the only indep. rubbish site(s)')
    # site_indeps = grep('season', indeps, invert = TRUE, value = TRUE)
    # site_indeps_log = paste0(site_indeps, '_log')

    # sites_nearby = sites_nearby %>%
    #     select(datetime, all_of(site_indeps), all_of(site_indeps_log))

    # if('season' %in% indeps){
    sites_nearby$season = factor(lubridate::quarter(sites_nearby$datetime))
    # }

    #assemble neon sensor data, filtered neon sensor data, neon field data
    #into one frame and plot it

    neon_q_auto = read_csv(glue('in/NEON/neon_continuous_Q/{neon_site}.csv')) %>%
        filter(! is.na(discharge)) %>%
        rename(discharge_auto = discharge)

    q_eval = read_csv('in/neon_q_eval.csv') %>%
        # filter(site == 'BLDE')
        filter(site == neon_site)

    # check1 = is.na(q_eval$regression_status) | q_eval$regression_status %in% c('good')
    # check2 = is.na(q_eval$drift_status) | q_eval$drift_status %in% c('likely_no_drift', 'not_assessed')
    # check3 = is.na(q_eval$rating_curve_status) | q_eval$rating_curve_status %in% c('Tier1')
    #
    # q_eval$keep = check1 & check2 & check3
    q_eval = q_eval %>%
        group_by(site, year, month) %>%
        summarize(keep = all(final_qual %in% c('Tier1', 'Tier2')), #if issues occur at any point in a month, reject that whole month
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

    dummies <- model.matrix(~qall$season)[, -1]
    colnames(dummies) <- paste0('season', 2:(ncol(dummies) + 1))
    qall_pred <- bind_cols(qall, dummies) %>%
        select(-season, -discharge_auto, -site_code,
               -discharge_auto_qc) %>%
        filter(if_all(everything(), ~! is.na(.)))

    qall_pred$predictions <- predict(best$best_model_object, data = qall_pred)$predictions
    qall_pred$predictions[qall_pred$predictions < 0] <- 0
    # qall_pred$predictions <- predict(best$best_model_object, data = qall_pred)$predictions

    qall = left_join(qall, select(qall_pred, datetime, predictions), by = 'datetime')

    out_data = qall %>%
        select(-ends_with('_log'), -any_of('season')) %>%
        full_join(select(best$rf_data, datetime, Q_neon_field), by = 'datetime') %>%
        # select(datetime, Q_predicted = fit,
        #        Q_pred_int_2.5 = lwr, Q_pred_int_97.5 = upr, Q_neon_field,
        #        Q_neon_continuous_filtered = discharge_auto_qc,
        #        Q_neon_continuous_raw = discharge_auto) %>%
        filter(if_any(-datetime, ~cumsum(! is.na(.)) != 0)) %>%  #remove leading NA rows
        arrange(datetime)

    dg = dygraphs::dygraph(xts(x = select(out_data, predictions, Q_neon_field) %>% tail(5e5),
                               order.by = tail(out_data$datetime, 5e5))) %>%
        dyRangeSelector()

    # saveWidget(dg, glue('figs/lm_plots/pred/{neon_site}_log.html'))

    return(dg)
}

# LSTM

run_lstm <- function(strategy, runset){

    #strategy: one of 'generalist', 'specialist', 'pdgl' (aka process-guided specialist)
    #runset: a list of run IDs, possibly split up into batches of runs as defined in
    #   in/lstm_configs

    strtgy <- ifelse(strategy == 'pgdl', 'specialist', strategy)

    #crude way to pass variables into python script from R
    r2pyenv_template <- new.env()
    r2pyenv_template$wdir <- getwd()
    r2pyenv_template$confdir <- confdir
    r2pyenv_template$rundir <- rundir
    r2pyenv_template$strategy <- strtgy
    r2pyenv_template$runset <- paste0('runs_', paste(range(runset), collapse = '-'))
    r2pyenv <- as.list(r2pyenv_template)

    py_run_file('src/lstm_dungeon/run_lstms.py')
}

# potential for future use

plots_and_results_daily <- function(neon_site, best, lm_df, results){

    stop('not updated')
    if(length(best$prediction) != nrow(lm_df)) stop('oi')

    d = lm_df %>%
        select(-ends_with('_log'), -season, -matches('[0-9]+'))

    #load corroborating usgs/ms site data

    sites_nearby = read_csv(glue('in/usgs_Q/{neon_site}.csv'))

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

    neon_q_daily = read_csv(glue('in/NEON/neon_continuous_Q/{neon_site}.csv')) %>%
        group_by(date = date(datetime)) %>%
        summarize(discharge = mean(discharge, na.rm = TRUE)) %>%
        filter(! is.na(discharge)) %>%
        rename(discharge_daily = discharge)

    q_eval = read_csv('in/neon_q_eval.csv') %>%
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

    saveWidget(dg, glue('figs/lm_plots/pred/{neon_site}_log.html'))

    if(inherits(best$fits, 'grob')){

        #make diagnostic and fit plots
        png(glue('figs/lm_plots/diag/{neon_site}_diag_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
        defpar = par(mfrow=c(2,2))
        plot(best$best_model_object)
        par(defpar)
        dev.off()

        png(glue('figs/lm_plots/fit/{neon_site}_fit_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
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

    png(glue('figs/lm_plots/val/{neon_site}_obs_v_pred.png'), 6, 6, 'in', type = 'cairo', res = 300)
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
    write_csv(out_data, glue('out/lm_out/by_site/{neon_site}.csv'))

    #return results frame, updated with this site's results
    results$bestmod_logq[results$site_code == neon_site] = as.character(best$best_model)[3]
    results$nse_logq[results$site_code == neon_site] = best$score
    results$nse_cv_logq[results$site_code == neon_site] = best$score_crossval

    return(results)
}

# misc helpers required to generate figures/datasets

dms_to_decdeg = function(x){

    #x: an integer or character vector of latlongs in dms format, e.g. "123456" or 1234567
    #                                                                   DDMMSS     DDDMMSS

    decdeg = c()
    for(i in seq_along(x)){

        xx = x[i]

        if(is.na(xx)){
            decdeg[i] = NA
            next
        }

        if(! nchar(xx) %in% 6:7){
            warning(paste(nchar(xx), 'characters in x. need 6 (or 7 for some longitudes)'))
            decdeg[i] = NA
        }

        deginc = if(nchar(xx) == 7) 1 else 0

        degs = as.numeric(substr(xx, 1, 2 + deginc))
        mins = as.numeric(substr(xx, 3 + deginc, 4 + deginc))
        secs = as.numeric(substr(xx, 5 + deginc, 6 + deginc))

        decdeg[i] = degs + mins / 60 + secs / 3600
    }

    return(decdeg)
}

gaugeid_to_location <- function(id){

    out <- dataRetrieval::readNWISsite(id)
    lat <- dms_to_decdeg(as.integer(out$lat_va))
    long <- -dms_to_decdeg(as.integer(out$long_va))

    sfobj <- tibble(lat = lat, long = long) %>%
        st_as_sf(coords = c('long', 'lat'), crs = 4269) %>%
        st_transform(crs = 4326)

    return(sfobj)
}

get_osm_streams <- function(extent_raster, outfile = NULL){

    #extent_raster: either a terra spatRaster or a rasterLayer. The output
    #   streams will have the same crs, and roughly the same extent, as this raster.
    #outfile: string. If supplied, output shapefile will be written to this
    #   location. If not supplied, the output will be returned.

    message('Downloading streams layer from OpenStreetMap')

    extent_raster <- terra::rast(extent_raster)
    # rast_crs <- as.character(extent_raster@crs)
    rast_crs <- terra::crs(extent_raster,
                           proj = TRUE)

    extent_raster_wgs84 <- terra::project(extent_raster,
                                          y = 'epsg:4326')

    dem_bounds <- terra::ext(extent_raster_wgs84)[c(1, 3, 2, 4)]

    streams_query <- osmdata::opq(dem_bounds) %>%
        osmdata::add_osm_feature(key = 'waterway',
                                 value = c('river', 'stream'))

    streams_query$prefix <- sub('timeout:25', 'timeout:180', streams_query$prefix)

    streams <- osmdata::osmdata_sf(streams_query)
    streams <- streams$osm_lines$geometry

    streams_proj <- streams %>%
        sf::st_transform(crs = rast_crs) %>%
        sf::st_union() %>%
        # sf::st_transform(crs = WGS84) %>%
        sf::st_as_sf() %>%
        rename(geometry = x) %>%
        mutate(FID = 0:(n() - 1)) %>%
        dplyr::select(FID, geometry)

    if(! is.null(outfile)){

        sf::st_write(streams_proj,
                     dsn = outfile,
                     layer = 'streams',
                     driver = 'ESRI Shapefile',
                     delete_layer = TRUE,
                     quiet = TRUE)

        message(paste('OSM streams layer written to', outfile))

    } else {
        return(streams_proj)
    }
}

locate_test_results <- function(nh_dir, runid){

    run_dirs <- list.files(nh_dir)

    rundir <- grep(paste0('^finetune', runid, '_'), run_dirs, value = TRUE)
    if(length(rundir) > 1){

        warning('multiple finetune dirs')
        return()

    } else if(length(rundir) == 0){

        rundir <- grep(paste0('^run', runid, '_'), run_dirs, value = TRUE)
        if(length(rundir) != 1){

            warning('0 or multiple run dirs')
            return()

        }
    }

    testdir <- file.path('test', list.files(file.path(nh_dir, rundir, 'test')))
    if(length(testdir) != 1){
        warning('0 or multiple test dirs')
        return()
    }

    res <- file.path(nh_dir, rundir, testdir, 'test_results.p')

    return(res)
}

retrieve_test_results <- function(runids){

    nse_out <- kge_out <-
        matrix(NA_real_,
               nrow = length(neon_sites),
               ncol = length(runids),
               dimnames = list(neon_sites, NULL))

    for(i in seq_along(runids)){

        td <- try(locate_test_results(nh_dir, runids[i]), silent = TRUE)
        if(inherits(td, 'try-error') || is.null(td)) next
        # write_lines(td, '/tmp/rundirs.txt', append = T)

        xx = reticulate::py_load_object(td)

        for(s in neon_sites){
            try({
                nse_out[rownames(nse_out) == s, i] <- xx[[paste0(s, '_MANUALQ')]]$`1D`$discharge_NSE
            }, silent = TRUE)
            try({
                kge_out[rownames(kge_out) == s, i] <- xx[[paste0(s, '_MANUALQ')]]$`1D`$discharge_KGE
            }, silent = TRUE)
        }
        # if(all(is.na(nse_out[, i]))) stop('no nse', i)
        # if(all(is.na(kge_out[, i]))) stop('no kge', i)
    }

    return(list(nse = nse_out, kge = kge_out))
}

polygon_with_gaps <- function(df){

    rl = rle(is.na(df$Q_pred_int_2.5))
    vv = ! rl$values
    chunkfac = rep(cumsum(vv), rl$lengths)
    chunkfac[chunkfac == 0] = 1
    chunks = split(df, chunkfac)
    noNAchunks = lapply(chunks, function(x) x[! is.na(x$Q_pred_int_2.5), ])

    for(i in 1:length(noNAchunks)){
        polygon(x=c(noNAchunks[[i]]$datetime, rev(noNAchunks[[i]]$datetime)),
                y=c(noNAchunks[[i]]$Q_pred_int_2.5, rev(noNAchunks[[i]]$Q_pred_int_97.5)),
                col=adjustcolor('red', alpha.f=0.2), border=NA)
    }
}

polygon_with_gaps2 <- function(df, gapcol, lowval, highval, col = 'blue',
                               invert = FALSE, border = col){

    if(invert){
        rl = rle(is.na(df[[gapcol]]))
    } else {
        rl = rle(! is.na(df[[gapcol]]))
    }

    vv = ! rl$values
    chunkfac = rep(cumsum(vv), rl$lengths)
    chunkfac[chunkfac == 0] = 1
    chunks = split(df, chunkfac)

    if(invert){
        NAchunks = lapply(chunks, function(x) x[! is.na(x[[gapcol]]), ])
    } else {
        NAchunks = lapply(chunks, function(x) x[is.na(x[[gapcol]]), ])
    }

    NAchunks = Filter(function(x) difftime(x$datetime[nrow(x)], x$datetime[1], units = 'hours') >= 6,
                      NAchunks)

    if(length(NAchunks)){

        for(i in 1:length(NAchunks)){

            d <- NAchunks[[i]]

            # if(nrow(d) == 1 || difftime(d$datetime[nrow(d)], d$datetime[1], units = 'mins') < 60){
            #     d <- bind_rows(d,
            #                    mutate(d, datetime = datetime[1] + 60 * 60))
            # }

            polygon(x=c(d$datetime, rev(d$datetime)),
                    y=c(rep(lowval, nrow(d)),
                        rep(highval, nrow(d))),
                    col=col, border=border)
                    # col=adjustcolor(col, alpha.f=0.2), border=NA)
        }
    }
}

ts_plot <- function(site, yr, boldgray = FALSE, ymax = Inf, scaled = TRUE){

    if(scaled){
        pred <- read_csv(paste0('out/lm_out_specQ/predictions/', site, '.csv')) %>%
            select(datetime, starts_with('Q_pred'), Q_neon_continuous_raw)
        fit <- read_csv(paste0('out/lm_out_specQ/fit/', site, '.csv')) %>%
            select(datetime, Q_neon_field)
    } else {
        pred <- read_csv(paste0('out/lm_out/predictions/', site, '.csv')) %>%
            select(datetime, starts_with('Q_pred'), Q_neon_continuous_raw)
        fit <- read_csv(paste0('out/lm_out/fit/', site, '.csv')) %>%
            select(datetime, Q_neon_field)
    }

    plotd <- full_join(pred, fit, by = 'datetime') %>%
        filter(datetime >= as.POSIXct(paste0(yr, '-01-01')),
               datetime <= as.POSIXct(paste0(yr, '-12-31'))) %>%
        arrange(datetime)

    ymax_ <- max(c(plotd$Q_neon_continuous_raw, plotd$Q_pred_int_97.5, plotd$Q_neon_field), na.rm = TRUE)
    rle_ <- rle2(as.numeric(month(plotd$datetime)))

    if(nrow(rle_) < 12){
        rle_ <- rle_[-1, ]
    }

    xaxis_labs <- substr(month.abb, 1, 1)[rle_$values]

    plot(plotd$datetime, plotd$Q_neon_continuous_raw, type = 'n',
         ylim = c(0, min(ymax_, ymax)),
         xlab = '', ylab = '', xaxt = 'n')
    axis(1, plotd$datetime[rle_$starts], xaxis_labs)
    axis(1, plotd$datetime[rle_$stops[nrow(rle_)]], yr, tick = F, font = 2)
    polygon_with_gaps(plotd)
    if(boldgray){
        lines(plotd$datetime, plotd$Q_neon_continuous_raw, col = 'gray50', lwd = 3)
    } else {
        lines(plotd$datetime, plotd$Q_neon_continuous_raw, col = 'gray50')
    }
    lines(plotd$datetime, plotd$Q_predicted, col = 'red')
    points(plotd$datetime, plotd$Q_neon_field, col = 'black', pch = 1)
    mtext(site, 3, -1, adj = 0.01, padj = 0.1)
}

rle2 <- function(x){

    r <- rle(x)
    ends <- cumsum(r$lengths)

    r <- tibble(values = r$values,
                starts = c(1, ends[-length(ends)] + 1),
                stops = ends,
                lengths = r$lengths)

    return(r)
}

reduce_results <- function(res, name, f, ...){

    r_ <- suppressWarnings(apply(res, 1, f, ...))

    r_[is.infinite(r_)] <- NA

    out <- tibble(site = names(r_),
                  x = unname(r_)) %>%
        rename(!!sym(name) := x)

    return(out)
}

approxjoin_datetime <- function(x,
                                y,
                                rollmax = '7:30',
                                keep_datetimes_from = 'x',
                                indices_only = FALSE){
    #direction = 'forward'){

    #x and y: macrosheds standard tibbles with only one site_code,
    #   which must be the same in x and y. Nonstandard tibbles may also work,
    #   so long as they have datetime columns, but the only case where we need
    #   this for other tibbles is inside precip_pchem_pflux_idw, in which case
    #   indices_only == TRUE, so it's not really set up for general-purpose joining
    #rollmax: the maximum snap time for matching elements of x and y.
    #   either '7:30' for continuous data or '12:00:00' for grab data
    #direction [REMOVED]: either 'forward', meaning elements of x will be rolled forward
    #   in time to match the next y, or 'backward', meaning elements of
    #   x will be rolled back in time to reach the previous y
    #keep_datetimes_from: string. either 'x' or 'y'. the datetime column from
    #   the corresponding tibble will be kept, and the other will be dropped
    #indices_only: logical. if TRUE, a join is not performed. rather,
    #   the matching indices from each tibble are returned as a named list of vectors..

    #good datasets for testing this function:
    # x <- tribble(
    #     ~datetime, ~site_code, ~var, ~val, ~ms_status, ~ms_interp,
    #     '1968-10-09 04:42:00', 'GSWS10', 'GN_alk', set_errors(27.75, 1), 0, 0,
    #     '1968-10-09 04:44:00', 'GSWS10', 'GN_alk', set_errors(21.29, 1), 0, 0,
    #     '1968-10-09 04:47:00', 'GSWS10', 'GN_alk', set_errors(21.29, 1), 0, 0,
    #     '1968-10-09 04:59:59', 'GSWS10', 'GN_alk', set_errors(16.04, 1), 0, 0,
    #     '1968-10-09 05:15:01', 'GSWS10', 'GN_alk', set_errors(17.21, 1), 1, 0,
    #     '1968-10-09 05:30:59', 'GSWS10', 'GN_alk', set_errors(16.50, 1), 0, 0) %>%
    # mutate(datetime = as.POSIXct(datetime, tz = 'UTC'))
    # y <- tribble(
    #     ~datetime, ~site_code, ~var, ~val, ~ms_status, ~ms_interp,
    #     '1968-10-09 04:00:00', 'GSWS10', 'GN_alk', set_errors(1.009, 1), 1, 0,
    #     '1968-10-09 04:15:00', 'GSWS10', 'GN_alk', set_errors(2.009, 1), 1, 1,
    #     '1968-10-09 04:30:00', 'GSWS10', 'GN_alk', set_errors(3.009, 1), 1, 1,
    #     '1968-10-09 04:45:00', 'GSWS10', 'GN_alk', set_errors(4.009, 1), 1, 1,
    #     '1968-10-09 05:00:00', 'GSWS10', 'GN_alk', set_errors(5.009, 1), 1, 1,
    #     '1968-10-09 05:15:00', 'GSWS10', 'GN_alk', set_errors(6.009, 1), 1, 1) %>%
    #     mutate(datetime = as.POSIXct(datetime, tz = 'UTC'))

    #tests
    if('site_code' %in% colnames(x) && length(unique(x$site_code)) > 1){
        stop('Only one site_code allowed in x at the moment')
    }
    if('var' %in% colnames(x) && length(unique(drop_var_prefix(x$var))) > 1){
        stop('Only one var allowed in x at the moment (not including prefix)')
    }
    if('site_code' %in% colnames(y) && length(unique(y$site_code)) > 1){
        stop('Only one site_code allowed in y at the moment')
    }
    if('var' %in% colnames(y) && length(unique(drop_var_prefix(y$var))) > 1){
        stop('Only one var allowed in y at the moment (not including prefix)')
    }
    if('site_code' %in% colnames(x) &&
       'site_code' %in% colnames(y) &&
       x$site_code[1] != y$site_code[1]) stop('x and y site_code must be the same')
    if(! rollmax %in% c('7:30', '12:00:00')) stop('rollmax must be "7:30" or "12:00:00"')
    # if(! direction %in% c('forward', 'backward')) stop('direction must be "forward" or "backward"')
    if(! keep_datetimes_from %in% c('x', 'y')) stop('keep_datetimes_from must be "x" or "y"')
    if(! 'datetime' %in% colnames(x) || ! 'datetime' %in% colnames(y)){
        stop('both x and y must have "datetime" columns containing POSIXct values')
    }
    if(! is.logical(indices_only)) stop('indices_only must be a logical')

    #deal with the case of x or y being a specialized "flow" tibble
    # x_is_flowtibble <- y_is_flowtibble <- FALSE
    # if('flow' %in% colnames(x)) x_is_flowtibble <- TRUE
    # if('flow' %in% colnames(y)) y_is_flowtibble <- TRUE
    # if(x_is_flowtibble && ! y_is_flowtibble){
    #     varname <- y$var[1]
    #     y$var = NULL
    # } else if(y_is_flowtibble && ! x_is_flowtibble){
    #     varname <- x$var[1]
    #     x$var = NULL
    # } else if(! x_is_flowtibble && ! y_is_flowtibble){
    #     varname <- x$var[1]
    #     x$var = NULL
    #     y$var = NULL
    # } else {
    #     stop('x and y are both "flow" tibbles. There should be no need for this')
    # }
    # if(x_is_flowtibble) x <- rename(x, val = flow)
    # if(y_is_flowtibble) y <- rename(y, val = flow)

    #data.table doesn't work with the errors package, so error needs
    #to be separated into its own column. also give same-name columns suffixes

    if('val' %in% colnames(x)){

        x <- x %>%
            # mutate(err = errors::errors(val),
            #        val = errors::drop_errors(val)) %>%
            rename_with(.fn = ~paste0(., '_x'),
                        .cols = everything()) %>%
            data.table::as.data.table()

        y <- y %>%
            # mutate(err = errors::errors(val),
            #        val = errors::drop_errors(val)) %>%
            rename_with(.fn = ~paste0(., '_y'),
                        .cols = everything()) %>%
            data.table::as.data.table()

    } else {

        if(indices_only){
            x <- rename(x, datetime_x = datetime) %>%
                # mutate(across(where(~inherits(., 'errors')),
                #               ~errors::drop_errors(.))) %>%
                data.table::as.data.table()

            y <- rename(y, datetime_y = datetime) %>%
                # mutate(across(where(~inherits(., 'errors')),
                #               ~errors::drop_errors(.))) %>%
                data.table::as.data.table()
        } else {
            stop('this case not yet handled')
        }

    }

    #alternative implementation of the "on" argument in data.table joins...
    #probably more flexible, so leaving it here in case we need to do something crazy
    # data.table::setkeyv(x, 'datetime')
    # data.table::setkeyv(y, 'datetime')

    #convert the desired maximum roll distance from string to integer seconds
    rollmax <- ifelse(test = rollmax == '7:30',
                      yes = 7 * 60 + 30,
                      no = 12 * 60 * 60)

    #leaving this here in case the nearest neighbor join implemented below is too
    #slow. then we can fall back to a basic rolling join with a maximum distance
    # rollmax <- ifelse(test = direction == 'forward',
    #                   yes = -rollmax,
    #                   no = rollmax)
    #rollends will move the first/last value of x in the opposite `direction` if necessary
    # joined <- y[x, on = 'datetime', roll = rollmax, rollends = c(TRUE, TRUE)]

    #create columns in x that represent the snapping window around each datetime
    x[, `:=` (datetime_min = datetime_x - rollmax,
              datetime_max = datetime_x + rollmax)]
    y[, `:=` (datetime_y_orig = datetime_y)] #datetime col will be dropped from y

    # if(indices_only){
    #     y_indices <- y[x,
    #                    on = .(datetime_y <= datetime_max,
    #                           datetime_y >= datetime_min),
    #                    which = TRUE]
    #     return(y_indices)
    # }

    #join x rows to y if y's datetime falls within the x range
    joined <- y[x, on = .(datetime_y <= datetime_max,
                          datetime_y >= datetime_min)]
    joined <- na.omit(joined, cols = 'datetime_y_orig') #drop rows without matches

    #for any datetimes in x or y that were matched more than once, keep only
    #the nearest match
    joined[, `:=` (datetime_match_diff = abs(datetime_x - datetime_y_orig))]
    joined <- joined[, .SD[which.min(datetime_match_diff)], by = datetime_x]
    joined <- joined[, .SD[which.min(datetime_match_diff)], by = datetime_y_orig]

    if(indices_only){
        y_indices <- which(y$datetime_y %in% joined$datetime_y_orig)
        x_indices <- which(x$datetime_x %in% joined$datetime_x)
        return(list(x = x_indices, y = y_indices))
    }

    #drop and rename columns (data.table makes weird name modifications)
    if(keep_datetimes_from == 'x'){
        joined[, c('datetime_y', 'datetime_y.1', 'datetime_y_orig', 'datetime_match_diff') := NULL]
        data.table::setnames(joined, 'datetime_x', 'datetime')
    } else {
        joined[, c('datetime_x', 'datetime_y.1', 'datetime_y', 'datetime_match_diff') := NULL]
        data.table::setnames(joined, 'datetime_y_orig', 'datetime')
    }

    #restore error objects, var column, original column names (with suffixes).
    #original column order
    # joined <- tibble::as_tibble(joined) %>%
    #     mutate(val_x = errors::set_errors(val_x, err_x),
    #            val_y = errors::set_errors(val_y, err_y)) %>%
    #     dplyr::select(-err_x, -err_y)
    # mutate(var = !!varname)

    # if(x_is_flowtibble) joined <- rename(joined,
    #                                      flow = val_x,
    #                                      ms_status_flow = ms_status_x,
    #                                      ms_interp_flow = ms_interp_x)
    # if(y_is_flowtibble) joined <- rename(joined,
    #                                      flow = val_y,
    #                                      ms_status_flow = ms_status_y,
    #                                      ms_interp_flow = ms_interp_y)

    # if(! sum(grepl('^val_[xy]$', colnames(joined))) > 1){
    #     joined <- rename(joined, val = matches('^val_[xy]$'))
    # }

    joined <- dplyr::select(joined,
                            datetime,
                            # matches('^val_?[xy]?$'),
                            # any_of('flow'),
                            starts_with('site_code'),
                            any_of(c(starts_with('var_'), matches('^var$'))),
                            any_of(c(starts_with('val_'), matches('^val$'))),
                            starts_with('ms_status_'),
                            starts_with('ms_interp_'))

    return(joined)
}

drop_var_prefix <- function(x){

    unprefixed <- substr(x, 4, nchar(x))

    return(unprefixed)
}

mode_interval_dt <- function(x){
    #mode interval in minutes for datetime vector
    Mode(diff(as.numeric(x)) / 60)
}

ms_write_netcdf <- function(df_list, path){

    vars_units0 <- c(discharge = 'cms', dayl = 's', prcp = 'mm/day',
                     srad = 'W/m2', swe = 'mm', tmax = 'C', tmin = 'C', vp = 'Pa',
                     pet = 'mm')

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

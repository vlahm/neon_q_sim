#general

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

# setup

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

# retrieval

get_neon_field_discharge <- function(neon_sites){

    dir.create('in/neon_field_Q')

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

        write_csv(q, glue('in/neon_field_Q/{s}.csv'))
    }
}

get_neon_inst_discharge <- function(neon_sites){

    dir.create('in/neon_continuous_Q')

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

        write_csv(q, glue('in/neon_continuous_Q/{s}.csv'))
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

    #this one allows all interactions if interactions == TRUE

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

    neon_q_auto = read_csv(glue('in/neon_continuous_Q/{neon_site}.csv')) %>%
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
               Q_neon_continuous_filtered = discharge_auto_qc,
               Q_neon_continuous_raw = discharge_auto) %>%
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

    neon_q_daily = read_csv(glue('in/neon_continuous_Q/{neon_site}.csv')) %>%
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
    results$bestmod_logq[results$site_code == neon_site] = paste(as.character(best1$best_model)[3],
                                                                 as.character(best1$best_model)[3],
                                                                 sep = ' & ')
    results$nse_logq[results$site_code == neon_site] = nse
    results$nse_cv_logq[results$site_code == neon_site] = NA_real_
    results$adj_r_squared[results$site_code == neon_site] = mean(summary(best1$best_model_object)$adj.r.squared,
                                                                 summary(best2$best_model_object)$adj.r.squared)

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

    neon_q_auto = read_csv(glue('in/neon_continuous_Q/{neon_site}.csv')) %>%
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

    neon_q_daily = read_csv(glue('in/neon_continuous_Q/{neon_site}.csv')) %>%
        group_by(date = date(datetime)) %>%
        summarize(discharge = mean(discharge, na.rm = TRUE)) %>%
        filter(! is.na(discharge)) %>%
        rename(discharge_daily = discharge)

    q_eval = read_csv('/in/neon_q_eval.csv') %>%
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

get_substance_by_cas = function(casrn_vec,
                                exclude_synonyms = TRUE,
                                include_former_and_incorrect_cas = FALSE,
                                include_related = FALSE){

    query_base = 'https://cdxnodengn.epa.gov/cdx-srs-rest/substances/cas?casList='

    casrn_str = casrn_vec %>%
        as.character() %>%
        str_replace('-', '') %>%
        paste(collapse='%7C')

    excl_syn = ifelse(exclude_synonyms, '&excludeSynonyms=true', '')
    incl_other = ifelse(include_former_and_incorrect_cas, '&includeOtherCas=true', '')
    incl_rel = ifelse(include_related, '&includeOtherTsn=true', '')

    query = glue('{qb}{cc}{es}{io}{ir}&qualifier=exact',
                 qb = query_base,
                 cc = casrn_str,
                 es = excl_syn,
                 io = incl_other,
                 ir = incl_rel)

    r = httr::GET(query)
    d = httr::content(r, as='text', encoding='UTF-8')
    d = jsonlite::fromJSON(d)

    return(d)
}

get_ef_tablechunk = function(table_name,
                             column_name=NULL,
                             operator=NULL,
                             column_value=NULL,
                             rows=NULL,
                             rtn_fmt,
                             timeout_,
                             timeout_action,
                             debug_){

    # table_name: char string; the name of an envirofacts table
    # column_name: char: optional; use to filter by column value (see operator)
    # operator: =, !=, <, >, BEGINNING, CONTAINING; if filtering, the filter expression's operator (with column_value)
    # column_value: if filtering, the value on the RHS of the filtering expression
    # rows: 'x:y'; max of 100,000
    # see query_ef_table

    query_base = 'https://data.epa.gov/efservice/'
    # query_base = 'http://iaspub.epa.gov/enviro/efservice/' #old? still in one of their examples

    table_names_ext = ''
    if(length(table_name) > 1) stop('table_name max length is 1 (maybe raise an issue with EPA about this? should be allowed length 3)')
    # if(length(table_names) > 1){ #deprecated for now
    #     table_names_ext = paste0('/',
    #                              paste(table_names[2:length(table_names)],
    #                                    collapse = '/'))
    #     table_names = table_names[1]
    # }

    if(! is.null(column_name)){

        cl = length(column_name)
        if(length(operator) != cl && length(column_value) != cl){
            stop('lengths of column_name, operator, and column_value must all be equal')
        }

        filter_str = ''
        for(i in seq_len(cl)){
            filter_str = paste(filter_str,
                               paste(column_name[i], operator[i], column_value[i],
                                     sep = '/'),
                               sep = '/')
        }
    } else {
        filter_str = ''
    }

    # query = 'https://data.epa.gov/efservice/ICIS_LIMIT/REF_POLLUTANT/SRS_ID/439216/rows/1:5/JSON'
    # query='https://data.epa.gov/efservice/ICIS_LIMIT/SRS_ID/=/5520/csv' ##***
    # query = glue('{qb}{tn}{cn}{op}{cv}{tnx}{ro}/csv',
    query = glue('{qb}{tn}{fs}{ro}/{fmt}',
                 qb=query_base,
                 tn=paste(table_name, collapse='/'),
                 # cn=ifelse(is.null(column_name), '', paste0('/', column_name)),
                 # op=ifelse(is.null(operator), '', paste0('/', operator)),
                 # cv=ifelse(is.null(column_value), '', paste0('/', column_value)),
                 # tnx=table_names_ext,
                 fs=filter_str,
                 ro=ifelse(! is.null(rows), paste0('/rows/', rows), ''),
                 fmt=rtn_fmt)

    if(debug_) browser()

    # query = 'https://data.epa.gov/efservice/FACILITIES/rows/40000:49999/csv'
    # query = 'https://data.epa.gov/efservice/FACILITIES/rows/50000:59999/csv'
    # query = 'https://data.epa.gov/efservice/FACILITIES/rows/18000:18010/csv'
    r = try(httr::GET(query, timeout(timeout_)),
            silent = TRUE)
    if(inherits(r, 'try-error')){
        if(timeout_action == 'skip' && grepl('Timeout was reached', r[[1]])){
            message('timout reached. skipping this chunk. see timeout_ parameter.')
            return(tibble())
        }
        stop('timout reached. try increasing the value passed to the timeout_ parameter')
    }
    d = httr::content(r, as='text', encoding='UTF-8')
    # sw(sm(read_csv(d, col_types = cols(.default = 'c')))) %>%
    #     select(-starts_with('...')) %>%
    #     rename_with(function(x) str_match(x, '\\.([^\\.]+)$')[, 2])

    #sometimes SQL returns a server-side error when csv or json is requested.
    #so, trying csv first, and json if that fails.

    if(rtn_fmt == 'csv'){
        dout <- try({
            sw(sm(read_csv(d, col_types = cols(.default = 'c')))) %>%
                select(-starts_with('...')) %>%
                rename_with(function(x) str_match(x, '\\.([^\\.]+)$')[, 2])
        }, silent = TRUE)
    }

    if(rtn_fmt == 'json' || inherits(dout, 'try-error')){

        query = sub('csv$', 'json', query)
        r = try(httr::GET(query, timeout(timeout_)),
                silent = TRUE)
        if(inherits(r, 'try-error')){
            if(timeout_action == 'skip' && grepl('Timeout was reached', r[[1]])){
                message('timout reached. skipping this chunk. see timeout_ parameter.')
                return(tibble())
            }
            stop('timout reached. try increasing the value passed to the timeout_ parameter')
        }
        d = httr::content(r, as='text', encoding='UTF-8')
        dout = as_tibble(jsonlite::fromJSON(d)) %>%
            mutate(across(where(~! is.character(.)), as.character))
    }

    return(dout)
}

get_response_1char = function(msg,
                              possible_chars,
                              subsequent_prompt = FALSE){

    #msg: character. a message that will be used to prompt the user
    #possible_chars: character vector of acceptable single-character responses
    #subsequent prompt: not to be set directly. This is handled by
    #   get_response_1char during recursion.

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

# nrows=1; maxrows=1e5
# nrows=259000; maxrows=1e4
get_chunksets = function(nrows, maxrows){

    # nrows: the number of rows in an envirofacts table. probably queried via query_ef_rows()
    # maxrows: the maximum number of rows that can be returned via the envirofacts REST-API

    options(scipen = 100)

    nchunks = nrows %/% maxrows
    rem = nrows %% maxrows
    if(rem > 0) nchunks = nchunks + 1

    chunk_starts = seq(0, nrows, maxrows)
    chunk_ends = 1:nchunks * maxrows - 1
    if(rem != 0){
        chunk_ends[nchunks] = chunk_starts[nchunks] + rem
    }

    chunk_sets = paste(chunk_starts, chunk_ends,
                       sep=':')

    options(scipen = 0)

    return(chunk_sets)
}

query_ef_rows = function(table_name,
                         column_name=NULL,
                         operator=NULL,
                         column_value=NULL){

    # table_name: char string; the name of an envirofacts table
    # column_name: char: optional; use to filter by column value (see operator)
    # operator: =, !=, <, >, BEGINNING, CONTAINING; if filtering, the filter expression's operator (with column_value)
    # column_value: if filtering, the value on the RHS of the filtering expression

    query_base = 'https://data.epa.gov/efservice/'

    if(length(table_name) > 1) stop('this function is for querying rows of a single table')

    if(! is.null(column_name)){

        cl = length(column_name)
        if(length(operator) != cl && length(column_value) != cl){
            stop('lengths of column_name, operator, and column_value must all be equal')
        }

        filter_str = ''
        for(i in seq_len(cl)){
            filter_str = paste(filter_str,
                               paste(column_name[i], operator[i], column_value[i],
                                     sep = '/'),
                               sep = '/')
        }
    } else {
        filter_str = ''
    }

    # query = glue('{qb}{tn}{cn}{op}{cv}/count/json',
    query = glue('{qb}{tn}{fs}/count/json',
                 qb=query_base,
                 tn=paste(table_name, collapse='/'),
                 fs = filter_str)
    # cn=ifelse(is.null(column_name), '', paste0('/', column_name)),
    # op=ifelse(is.null(operator), '', paste0('/', operator)),
    # cv=ifelse(is.null(column_value), '', paste0('/', column_value)))

    r = httr::GET(query)
    d = httr::content(r, as='text', encoding='UTF-8')
    d = jsonlite::fromJSON(d)[[1]]

    return(d)
}

query_ef_table = function(table_name,
                          column_name=NULL,
                          operator=NULL,
                          column_value=NULL,
                          specify_rows=TRUE,
                          chunk_size=1e4,
                          custom_chunks=NULL,
                          rtn_fmt='csv',
                          timeout_=120,
                          timeout_action='fail',
                          warn=TRUE,
                          verbose=FALSE,
                          debug_=FALSE){

    # table_name: char string; the name of an envirofacts table
    # column_name: char: optional; use to filter by column value (see operator)
    #   if filtering by multiple columns, include them all here as a character vector.
    # operator: =, !=, <, >, BEGINNING, CONTAINING; if filtering, the filter
    #   expression's operator (with column_value). If filtering on multiple columns,
    #   must supply a separate operator for each column as a character vector.
    #   Make sure elements of column_name, operator, and column_value line
    #   up correctly.
    # column_value: if filtering, the value on the RHS of the filtering expression.
    #   If filtering on multiple columns, must supply a value for each column as
    #   a vector. Make sure elements of column_name, operator, and column_value line
    #   up correctly.
    # specify_rows: logical. if TRUE, the API request will include "/rows/x:y".
    #   There's no reason to change the default unless you experience hanging requests
    #   or unexplainable request failures. Something might be misspecified server-
    #   side, and this could help. Note though that a request can't be larger
    #   than 100,000 rows. If your query would return more than this,
    #   specify_rows must be TRUE.
    # chunk_size: the maximum number of rows to ask for at a time. Envirofacts
    #   won't let you ask for more than 100,000 (the default).
    # custom_chunks: if all else fails, you can figure out exactly which chunks are
    #   causing trouble and tell get_ef_tablechunk exactly how to split
    #   up its requests. Amazingly, there are some chunks that can't be
    #   retrieved whole, but can be retrieved as two separate requests, notably
    #   https://data.epa.gov/efservice/FACILITIES/rows/31500:31600/csv will fail,
    #   but 31500:31550 and 31550:31600 will both work. or, 31500:31600/json will
    #   work. if specified, this argument overrides chunk_size.
    # rtn_fmt: character. semi-deprecated. This function will now try to return
    #   CSV results first, and if that fails, it'll try again with JSON results.
    #   You may also specify "json" to try that first, but then it won't fall back to CSV.
    #   see details
    # timeout_: integer. the amount of time to wait for a response from the server.
    # timeout_action: character. either "skip" to move on to the next chunk,
    #   or "fail" to raise an error.
    # warn: logical. if TRUE and more than 1 million rows are queried, this function will
    #   ask for confirmation before proceeding.
    # verbose: logical. if TRUE, you'll receive more information.
    # debug_: logical; only for interactive use. If TRUE, the debugger will be
    #   entered before the first request is executed.

    #DETAILS
    #Some envirofacts tables generate errors when returning one format or the
    #   other (JSON or CSV). The output of this function will be a
    #   data.frame (tibble) in either case. This parameter is just a way of fiddling with
    #   an imperfect data retrieval system.

    #RETURN VALUE
    #an all-character tibble of evirofacts results

    if(! rtn_fmt %in% c('csv', 'json')){
        stop('rtn_fmt must be "csv" or "json"')
    }

    if(! timeout_action %in% c('skip', 'fail')){
        stop('timeout_action must be "skip" or "fail"')
    }

    nrows = query_ef_rows(table_name = table_name,
                          column_name = column_name,
                          operator = operator,
                          column_value = column_value)

    if(warn && nrows > 1e6){
        resp = get_response_1char(msg = paste('This table has',
                                              nrows,
                                              'rows. still load into memory? (y/n) >'),
                                  c('y', 'n'))
        if(resp == 'n'){
            message('aborted query')
            rtn_abrt = 'user aborted'
            class(rtn_abrt) = 'ef_err'
            return(rtn_abrt)
        }
    } else {
        if(verbose) print(paste('This table has', nrows, 'row(s).'))
    }

    if(nrows == 0){
        if(verbose) message('Returning empty tibble')
        return(tibble())
    }

    if(is.null(custom_chunks)){
        chunksets = get_chunksets(nrows=nrows,
                                  maxrows=chunk_size)
    } else {
        chunksets <- custom_chunks
    }

    full_table = tibble()
    task_start = proc.time()
    if(specify_rows){
        for(i in seq_along(chunksets)){

            if(verbose) print(paste0('Retrieving chunk ', i, ' of ', length(chunksets),
                                     '; (', chunksets[i], ')'))

            chnk = get_ef_tablechunk(table_name=table_name,
                                     column_name=column_name,
                                     operator=operator,
                                     column_value=column_value,
                                     rows=chunksets[i],
                                     rtn_fmt=rtn_fmt,
                                     timeout_=timeout_,
                                     timeout_action=timeout_action,
                                     debug_=debug_)

            full_table = bind_rows(full_table, chnk)
        }
    } else {
        if(nrows >= 1e5) stop('envirofacts does not allow queries of more than 100,000 rows.')

        if(verbose) print(paste('attempting to retrieve all rows (specify_rows is FALSE).'))

        full_table = get_ef_tablechunk(table_name=table_name,
                                       column_name=column_name,
                                       operator=operator,
                                       column_value=column_value,
                                       rows=NULL,
                                       rtn_fmt=rtn_fmt)
    }

    task_time = unname(round((proc.time() - task_start)[3] / 60, 2))
    if(verbose) print(paste0('Got table ', table_name, ' (', nrows, ' rows) in ', task_time, ' minutes'))

    return(full_table)
}

# table_name = c('ICIS_LIMIT'); rows = '1:10'
# table_names = c('tri_facility', 'tri_reporting_form', 'tri_chem_info'); rows = '1:10'
# column_name = 'state_abbr'; operator = '='; column_value = 'VA'
# return_count = FALSE
# table_name = c('ICIS_LIMIT'); rows = '1:10'

# query_ef_rows('XREF_AIR_FAC_INT_PROG_SUBPART')
# xx = query_ef_table('ref_sic')

query_envirofacts = function(table_names,
                             filter_column=NULL,
                             filter_operator=NULL,
                             filter_value=NULL,
                             join_column,
                             save_intermediates=FALSE,
                             warn=TRUE){

    # table_names: char vector of envirofacts table names
    # filter_column: char: optional; use to filter by column value (see filter_operator). must be shared among all tables in table_names.
    # filter_operator: =, !=, <, >, BEGINNING, CONTAINING; if filtering, the filter expression's filter_operator (with filter_value)
    # filter_value: if filtering, the value on the RHS of the filtering expression
    # join_column: the column by which to join tables in table_names. must be shared by all.
    # save_intermediates: logical. save each table's results to global env, in case something goes wrong?
    #   DEPRECATED? (not yet): using arrow to manipulate files on disk, so saving intermediates to disk now required.

    #combines results from several envirofacts tables.

    combined_results = tibble()
    for(i in seq_along(table_names)){

        tname = table_names[i]

        print(paste('Working on', tname))

        table_result = query_ef_table(tname,
                                      column_name=filter_column,
                                      operator=filter_operator,
                                      column_value=filter_value,
                                      warn=warn)

        if(inherits(table_result, 'ef_err')) return()

        if(save_intermediates){
            print(paste('saving table', tname, 'to global env'))
            assign(tname, table_result, pos=.GlobalEnv)
        }

        if(nrow(table_result)){

            shared_cols = intersect(colnames(combined_result),
                                    colnames(table_result))
            shared_cols = shared_cols[shared_cols != join_column]

            if(length(shared_cols)){
                message(paste('dropping redundant columns from :',
                              paste(shared_cols, collapse = ', ')))
            }

            table_result = select(table_result, -any_of(shared_cols))
            combined_results = full_join(combined_result, table_result,
                                         by = toupper(join_column))
        } else {

            if(! is.null(filter_column)){
                warning(paste('no rows in', tname, 'for', filter_column,
                              operator, filter_value))
            } else {
                warning(paste('no rows in', tname))
            }
        }

    }

    return(combined_results)
}

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
    long <- dms_to_decdeg(as.integer(out$long_va))

    sfobj <- tibble(lat = lat, long = long) %>%
        st_as_sf(coords = c('long', 'lat'), crs = 4269) %>%
        st_transform(crs = 4326)

    return(sfobj)
}

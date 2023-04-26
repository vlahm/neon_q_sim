#helper functions for the other R scripts. No need to source these directly.
#start with 01_linear_regression.R

#general (regression) ####

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

# setup (regression) ####

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

# retrieval ####

get_neon_field_discharge <- function(neon_sites){

    dir.create('in/NEON', showWarnings = FALSE)
    dir.create('in/NEON/neon_field_Q', showWarnings = FALSE)

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

        write_csv(q, glue('in/NEON/neon_field_Q/{s}.csv'))
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

# data prep ####

assemble_q_df <- function(neon_site, nearby_usgs_gages = NULL, ms_Q_data = NULL,
                          datetime_snapdist_hrs = 12, overwrite = FALSE,
                          scale_q_by_area = TRUE){

    #ms_Q_data: data.frame with posixct datetime column, and Q columns. Each Q column must
    #have its MacroSheds site_code as header. transformed columns will be created

    #can't remember if i built this to handle nearby_usgs_gages and ms_Q_data ,
    #but in any case that wasn't relevant for this study.

    neon_q_manual = read_csv(glue('in/NEON/neon_field_Q/{neon_site}.csv')) %>%
        mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
        rename(discharge_manual = discharge) %>%
        distinct(datetime, .keep_all = TRUE)
    earliest_date = as.character(date(min(neon_q_manual$datetime)))

    if(scale_q_by_area){
        wsa = filter(ms_areas, site_code == neon_site) %>% pull(ws_area_ha)
        neon_q_manual$discharge_manual = neon_q_manual$discharge_manual / wsa * 1000 #scale by 1000 so that neglog skew is negligible
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
                                overwrite = FALSE){

    #ms_Q_data: data.frame with date column, and Q columns. Each Q column must
    #have its site_code as header. transformed columns will be created

    neon_q_manual = read_csv(glue('in/NEON/neon_field_Q/{neon_site}.csv')) %>%
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

# regression ####

glmnet_wrap <- function(data, full_spec, unscale_q_by_area = TRUE,
                        k_folds = 10, cv_ensemble = 10, plot_marg = TRUE,
                        alpha = 0, ...){

    #data: a data frame with a date column, and all columns referenced in model_list
    #full_spec: formula; the fully specified model to regularize via L2 (ridge regression)
    #cv_ensemble: the number of cross-validation procedures to carry out. optimal
    #   lambda varies according to how observations are randomly allocated to CV folds.
    #   so, select `cv_ensemble` optimal lambdas and average them.
    #plot_marg: logical. plot marginal relationships
    #alpha: 0 = ridge; 1 = lasso; see ?glmnet
    #...: additional arguments passed to glmnet AND cv.glmnet

    ## run regression with 10-fold crossval

    x <- model.matrix(full_spec, data = data)
    mod_inds <- complete.cases(data)
    y <- data$discharge_log[mod_inds]
    data_filt <- data[mod_inds, ]

    ## get cross-validated predictions and metrics

    best_mse <- Inf
    neon_site <- data$site_code[1]
    wsa <- filter(ms_areas, site_code == neon_site) %>% pull(ws_area_ha)
    for(i in 0:1){

        n_samp <- length(y)
        # fold_ids <- sample(rep(1:10, length.out = n_samp), size = n_samp)

        sim_mat <- obs_mat <- matrix(NA_real_, nrow = n_samp, ncol = cv_ensemble)
        lambda_vec <- rep(NA_real_, cv_ensemble)
        for(j in 1:cv_ensemble){

            cv_model_ <- cv.glmnet(x, y, alpha = alpha, n_folds = k_folds, #foldid = fold_ids
                                   intercept = as.logical(i), type.measure = 'mse', ...)
            lambda_vec[j] <- cv_model_$lambda.min
            cvpred <- c(predict(cv_model_, s = 'lambda.min', newx = x))

            if(unscale_q_by_area){
                sim_mat[, j] <- inv_neglog(cvpred) * wsa / 1000 #un-scale by 1000. see assemble_q_df
                obs_mat[, j] <- data_filt$discharge * wsa / 1000
            } else {
                sim_mat[, j] <- inv_neglog(cvpred)
                obs_mat[, j] <- data_filt$discharge
            }
        }

        best_lambda_ <- mean(lambda_vec)
        sim_ <- apply(sim_mat, 1, mean, na.rm = TRUE)
        obs_ <- apply(obs_mat, 1, mean, na.rm = TRUE)

        # print(tibble(nse = hydroGOF::NSE(sim_, obs_),
        #              kge = hydroGOF::KGE(sim_, obs_)))

        new_mse <- hydroGOF::mse(sim_, obs_)
        if(new_mse < best_mse){
            best_mse <- new_mse
            # cv_model <- cv_model_
            sim <- sim_
            obs <- obs_
            best_lambda <- best_lambda_
            intcpt <- as.logical(i)
        }
    }

    metr_cv <- tibble(
        nse = hydroGOF::NSE(sim, obs),
        kge = hydroGOF::KGE(sim, obs),
        mse = hydroGOF::mse(sim, obs),
        mae = hydroGOF::mae(sim, obs),
        pbias = hydroGOF::pbias(sim, obs)
    )

    ## plot marginal associations

    trms <- rownames(attributes(terms(full_spec))$factors)
    dep <- gsub('`', '', trms[1])
    indeps <- gsub('`', '', trms[-1])
    site_indeps <- grep('season', indeps, invert = TRUE, value = TRUE)

    ggps <- list()
    for(i in seq_along(site_indeps)){
        d2 <- filter(data, ! is.na(discharge_log) & ! is.na(!!sym(site_indeps[i])))
        ggps[[i]] <- ggplot(d2, aes(x = !!sym(site_indeps[i]), y = discharge_log)) +
            geom_point()
    }
    gd <- do.call("grid.arrange", c(ggps))
    if(plot_marg) print(gd)

    ## handle novel seasonal factor values

    if('season' %in% indeps){
        modeled_seasons <- sort(as.character(unique(data_filt$season)))
        modseas_inds <- as.character(data_filt$season) %in% modeled_seasons
    } else {
        modeled_seasons <- as.character(1:4)
        modseas_inds <- rep(TRUE, nrow(data_filt))
    }

    ## run regression with optimal lambda and generate predictions + metrics

    best_model <- glmnet(x, y, alpha = alpha, lambda = best_lambda, intercept = intcpt, ...)

    data$lm <- NA_real_
    bestpred <- predict(best_model, s = best_lambda, newx = x[modseas_inds, ])

    if(unscale_q_by_area){
        data$lm[mod_inds][modseas_inds] <- inv_neglog(bestpred) * wsa / 1000
        data$discharge <- data$discharge * wsa / 1000
    } else {
        data$lm[mod_inds][modseas_inds] <- inv_neglog(bestpred)
    }

    data$lm[data$lm < 0] <- 0

    metr <- tibble(
        nse = hydroGOF::NSE(data$lm, data$discharge),
        kge = hydroGOF::KGE(data$lm, data$discharge),
        mse = hydroGOF::mse(data$lm, data$discharge),
        mae = hydroGOF::mae(data$lm, data$discharge),
        pbias = hydroGOF::pbias(data$lm, data$discharge)
    )

    ## trim leading/trailing NAs, clean up data for plotting
    first_non_na <- Position(function(x) ! is.na(x), data$lm)
    last_non_na <- nrow(data) - Position(function(x) ! is.na(x), rev(data$lm)) + 1
    plot_data <- data[first_non_na:last_non_na, ]

    site_indeps <- c(site_indeps, sub('_log', '', site_indeps))
    drop_cols <- grep('^[0-9]', colnames(plot_data), value = TRUE)
    drop_cols <- drop_cols[! drop_cols %in% site_indeps]
    plot_data <- select(plot_data, -any_of(drop_cols), -ends_with('_log')) %>%
        select(site_code, any_of(c('date', 'datetime')), Q_neon_field = discharge, Q_predicted = lm,
               everything())

    ## unscale usgs Q

    if(unscale_q_by_area){

        nearby_usgs_gages <- grep('^[0-9]+$', colnames(plot_data), value = TRUE)
        for(g in nearby_usgs_gages){

            if(g == '06190540'){
                plot_data[[g]] <- plot_data[[g]] * 51089.73 #ha
                next
            }

            nwissite <- dataRetrieval::readNWISsite(g)

            if(is.null(nwissite)) stop('cannot find nwis watershed area')

            if(is.na(nwissite$contrib_drain_area_va)){
                wsa <- nwissite$drain_area_va * 258.999 #mi^2 -> ha
            } else {
                wsa <- nwissite$contrib_drain_area_va * 258.999 #mi^2 -> ha
            }

            plot_data[[g]] <- plot_data[[g]] * wsa / 1000
        }

        nearby_ms_gages <- select(plot_data,
                                  -site_code,
                                  -any_of(c('date', 'datetime')),
                                  -starts_with('Q_'),
                                  -season,
                                  -any_of(nearby_usgs_gages)) %>%
            colnames()

        for(g in nearby_ms_gages){
            wsa <- filter(ms_areas, site_code == g) %>% pull(ws_area_ha)
            plot_data[[g]] <- plot_data[[g]] * wsa / 1000
        }
    }

    full_spec_ignoreint <- full_spec
    if(! intcpt){
        full_spec <- update(full_spec, ~ . - 1)
    }

    out <- list(best_model = full_spec,
                best_model_ignoreint = full_spec_ignoreint,
                best_model_object = best_model,
                alpha = alpha,
                prediction = unname(data$lm),
                lm_data = plot_data,
                metrics = metr,
                metrics_crossval = metr_cv,
                fits = gd)

    return(out)
}

bootstrap_ci_glmnet <- function(ncores, in_df, frm, best, has_intcpt, newx_){

    clst_type <- ifelse(.Platform$OS.type == 'windows', 'PSOCK', 'FORK')
    # ncores <- parallel::detectCores() %/% 4.5
    clst <- makeCluster(spec = ncores, type = clst_type)
    registerDoParallel(clst)

    bootstrap_samps <- parallel::parSapply(clst, 1:1000, function(x){

        resamp <- stratified_resample(in_df, 'season')
        x <- model.matrix(frm, data = resamp)
        mod_inds <- complete.cases(resamp)
        y <- resamp$discharge_log[mod_inds]

        s <- cv.glmnet(x, y, alpha = best$alpha, intercept = has_intcpt)$lambda.min
        m <- glmnet(x, y, alpha = best$alpha, intercept = has_intcpt, lambda = s)

        predict(m, s = s, newx = model.matrix(update(frm, NULL ~ .), newx_))
    }, simplify = TRUE)

    ci <- parallel::parApply(clst, bootstrap_samps, 1, quantile, probs = c(0.025, 0.975))

    # ci_lwr <- do.call(function(...) Map(function(...) parallel::clusterMap(
    #     clst, quant025, ..., SIMPLIFY = TRUE, USE.NAMES = FALSE), ...),
    #     bootstrap_samps)
    # ci_upr <- do.call(function(...) Map(function(...) parallel::clusterMap(
    #     clst, quant975, ..., SIMPLIFY = TRUE, USE.NAMES = FALSE), ...),
    #     bootstrap_samps)

    stopCluster(clst)

    return(list(ci_lwr = ci[1, ], ci_upr = ci[2, ]))
}

# quant975 <- function(...) quantile(c(...), probs = 0.975)
# quant025 <- function(...) quantile(c(...), probs = 0.025)

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

    #if one of the terms is "dummy" it will not be considered for interaction

    #returns a list of formulae

    trms = rownames(attributes(terms(full_spec))$factors)

    dummy_present <- 'dummy' %in% trms
    if(dummy_present){
        trms <- grep('dummy', trms, value = TRUE, invert = TRUE)
    }

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

    if(dummy_present){
        mods <- c(mods, lapply(mods, function(x) update(x, . ~ . + dummy)))
    }

    param_counts = sapply(mods, function(x) sum(grepl('season', attr(terms(x), 'term.labels')) * 2 + 1)) + 1
    is_through_origin = sapply(mods, function(x) grepl('^0 +', as.character(x)[3]))
    param_counts = param_counts - is_through_origin
    mods = mods[param_counts <= max_params]

    mods <- Filter(function(x){
        trms <- attr(terms(x), 'term.labels')
        ! all(trms %in% c('season', 'dummy'))
    }, mods)

    return(mods)
}

lm_wrap <- function(data, model_list,
                    unscale_q_by_area = TRUE, k_folds = 10){

    #data: a data frame with a date column, and all columns referenced in model_list
    #model_list: a list of lm formulae, as generated by generate_nested_formulae

    #lots of vestigial code in this function. sincere apologies.

    log = 'xy'
    best_mod = NA
    best_score = Inf
    for(i in seq_along(model_list)){

        try({

            trms = rownames(attributes(terms(model_list[[i]]))$factors)
            indeps = gsub('`', '', trms[-1])
            dep = gsub('`', '', trms[1])
            dep_transformed = sub('_log', '', dep)

            dd = filter(data, if_all(all_of(c(dep, indeps)), ~ ! is.na(.)))

            tryCatch({
                dd = CVlm(data = dd, model_list[[i]], m = k_folds,
                          plotit = FALSE, printit = FALSE)
            }, warning = function(w) stop('rank-deficiency detected. moving on.'))

            if(log == 'xy'){

                neon_site = data$site_code[1]
                if(unscale_q_by_area){
                    wsa = filter(ms_areas, site_code == neon_site) %>% pull(ws_area_ha)
                    sim = inv_neglog(dd$cvpred) * wsa / 1000
                    # nsenum = sum((inv_neglog(dd$cvpred) * wsa / 1000 - inv_neglog(dd[[dep]]) * wsa / 1000)^2) #just checking
                } else {
                    sim = inv_neglog(dd$cvpred)
                    # nsenum = sum((inv_neglog(dd$cvpred) - inv_neglog(dd[[dep]]))^2)
                }
            }

            if(unscale_q_by_area){
                obs = pull(dd[dep_transformed]) * wsa / 1000
            } else {
                obs = pull(dd[dep_transformed])
            }

            # nsedenom = sum((obs - mean(obs))^2)
            # nse_ = 1 - (nsenum / nsedenom)

            # nse = hydroGOF::NSE(sim, obs)
            mserr = hydroGOF::mse(sim, obs)

            # if(! all.equal(nse, nse_)) stop('oi')

            if(mserr < best_score){
                best_score = mserr
                best_mod = i

                metr_cv = tibble(
                    nse = hydroGOF::NSE(sim, obs),
                    kge = hydroGOF::KGE(sim, obs),
                    mse = mserr,
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
    site_indeps = grep('season|dummy', indeps, invert = TRUE, value = TRUE)

    dd = filter(data, if_all(all_of(c(dep, indeps)), ~ ! is.na(.)))
    m = lm(model_list[[best_mod]], data = dd)

    ggps = list()
    yvar = ifelse(log %in% c('y', 'xy'), 'discharge_log', 'discharge')
    for(i in seq_along(site_indeps)){
        d2 = filter(data, ! is.na(!!sym(yvar)) & ! is.na(!!sym(site_indeps[i])))
        ggps[[i]] = ggplot(d2, aes(x = !!sym(site_indeps[i]), y = !!sym(yvar))) +
            geom_point()

        if(! 'dummy' %in% indeps){
            ggps[[i]] <- ggps[[i]] + stat_smooth(method = "lm", col = "red")
        }
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
        nse = nse_out,
        kge = hydroGOF::KGE(data$lm, data$discharge),
        mse = hydroGOF::mse(data$lm, data$discharge),
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

        nearby_ms_gages = select(plot_data, -site_code, -any_of(c('date', 'datetime', 'dummy')),
                                 -starts_with('Q_'), -season, -any_of(nearby_usgs_gages)) %>%
            colnames()

        for(g in nearby_ms_gages){
            wsa = filter(ms_areas, site_code == g) %>% pull(ws_area_ha)
            plot_data[[g]] = plot_data[[g]] * wsa / 1000
        }
    }

    out = list(best_model = model_list[[best_mod]],
               best_model_object = m,
               prediction = unname(data$lm),
               lm_data = plot_data,
               # score = nse_out,
               # score_crossval = best_score,
               metrics = metr,
               metrics_crossval = metr_cv,
               # plot = dg,
               fits = gd)

    return(out)
}

segmented_wrap <- function(data, model_list,
                           unscale_q_by_area = TRUE, k_folds = 10){

    #data: a data frame with a date column, and all columns referenced in model_list
    #model_list: a list of lm formulae, as generated by generate_nested_formulae

    best_mod <- NA
    best_score <- Inf
    for(i in seq_along(model_list)){

        try({

            trms <- rownames(attributes(terms(model_list[[i]]))$factors)
            indeps <- gsub('`', '', trms[-1])
            dep <- gsub('`', '', trms[1])
            dep_transformed <- sub('_log', '', dep)

            dd <- filter(data, if_all(all_of(c(dep, indeps)), ~ ! is.na(.)))

            tryCatch({
                dd <- CVlm(data = dd, model_list[[i]], m = k_folds,
                          plotit = FALSE, printit = FALSE)
            }, warning = function(w) stop('rank-deficiency detected. moving on.'))

            neon_site <- data$site_code[1]
            if(unscale_q_by_area){
                wsa <- filter(ms_areas, site_code == neon_site) %>% pull(ws_area_ha)
                sim <- inv_neglog(dd$cvpred) * wsa / 1000
            } else {
                sim <- inv_neglog(dd$cvpred)
            }

            if(unscale_q_by_area){
                obs <- pull(dd[dep_transformed]) * wsa / 1000
            } else {
                obs <- pull(dd[dep_transformed])
            }

            mserr <- hydroGOF::mse(sim, obs)

            if(mserr < best_score){
                best_score <- mserr
                best_mod <- i
                metr_cv <- tibble(nse = NA_real_, kge = NA_real_, mse = mserr,
                                  mae = NA_real_, pbias = NA_real_)
            }
        })
    }

    trms <- rownames(attributes(terms(model_list[[best_mod]]))$factors)
    dep <- gsub('`', '', trms[1])
    indeps <- gsub('`', '', trms[-1])
    if(length(indeps) == 1 && indeps == 'season') stop('season is the only indep. rubbish site(s)')
    site_indeps <- grep('season', indeps, invert = TRUE, value = TRUE)

    #segmented package can't handle variable names with backticks
    model_list[[best_mod]] <- model_list[[best_mod]] %>%
        as.character() %>%
        str_replace_all('`?([0-9]{2,}(?:_log)?)`?', 'x\\1') %>%
        {paste(.[c(2, 1, 3)], collapse = ' ')} %>%
        as.formula()

    data_filt <- filter(data, if_all(all_of(c(dep, indeps)), ~ ! is.na(.))) %>%
        rename_with(~sub('`?([0-9]{2,}_log)`?', 'x\\1', .))

    seg_vars <- update(model_list[[best_mod]], NULL ~ . - season)
    if(! attr(terms(model_list[[best_mod]]), 'intercept')){
        seg_vars <- update(seg_vars, NULL ~ . + 1)
    }

    m <- lm(model_list[[best_mod]], data = data_filt) %>%
        segmented(seg.Z = seg_vars, npsi = 1)

    ggps = list()
    for(i in seq_along(site_indeps)){
        d2 = filter(data, ! is.na(discharge_log) & ! is.na(!!sym(site_indeps[i])))
        ggps[[i]] = ggplot(d2, aes(x = !!sym(site_indeps[i]), y = discharge_log)) +
            geom_point()
    }
    gd = do.call("grid.arrange", c(ggps))
    print(gd)

    newdata = select(data, all_of(indeps)) %>%
        rename_with(~sub('`?([0-9]{2,}_log)`?', 'x\\1', .))

    if('season' %in% indeps){
        modeled_seasons = m$xlevels$season
        modseas_inds = as.character(newdata$season) %in% modeled_seasons
    } else {
        modeled_seasons = as.character(1:4)
        modseas_inds = rep(TRUE, nrow(newdata))
    }

    data$lm = NA_real_

    if(unscale_q_by_area){
        data$lm[modseas_inds] = inv_neglog(predict(m, newdata = newdata[modseas_inds, ])) * wsa / 1000
        data$discharge = data$discharge * wsa / 1000
    } else {
        data$lm[modseas_inds] = inv_neglog(predict(m, newdata = newdata[modseas_inds, ]))
    }

    data$lm[data$lm < 0] = 0

    metr = tibble(
        nse = hydroGOF::NSE(data$lm, data$discharge),
        kge = hydroGOF::KGE(data$lm, data$discharge),
        mse = hydroGOF::mse(data$lm, data$discharge),
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

    out = list(best_model = model_list[[best_mod]],
               best_model_object = m,
               prediction = unname(data$lm),
               lm_data = plot_data,
               metrics = metr,
               metrics_crossval = metr_cv,
               fits = gd)

    return(out)
}

regress <- function(neon_site, framework, ..., scale_q_by_area = TRUE,
                    precomputed_df = NULL, custom_formula = NULL,
                    dummy_break = NULL, custom_gauges = NULL,
                    no_write = FALSE, bootstrap_ci = TRUE,
                    ncores = max(parallel::detectCores() %/% 4, 2)){

    #neon_site: 4-letter NEON site code, e.g. "REDB"
    #framework: string. 'lm', 'glmnet', or 'segmented'
    #...: other arguments to lm_wrap, glmnet_wrap, or segmented_wrap
    #     depending on the value of `framework`
    #scale_q_by_area: logical. should discharge be scaled by watershed area?
    #precomputed_df: data.frame as produced by assemble_q_df or
    #   assemble_q_df_daily. only supply this if you need to make modifications.
    #   otherwise it will be generated automatically by this function. may
    #   include a factor column called "dummy" that contains a dummy variable
    #custom_formula: a `glue` specification . if not provided, the fully specified model
    #   is assumed. See function definition for examples. only set up to work
    #   with special cases of framework == 'lm'.
    #dummy_break: if a dummy variable is used, this specifies the date on which to
    #   shift from dummy = 0 to dummy = 1. ignored otherwise
    #custom_gauges: character vector. use if you need to override the defaults
    #   provided by cft/donor_gauges.yml, e.g. if for building composite models.
    #no_write: logical. if TRUE, do not write to the `results` data.frame in the
    #   global environment, and do not write data or plots via plots_and_results.
    #bootstrap_ci: logical. should 95% confidence intervals for glmnet models
    #   be computed? Ignored if framework != 'glmnet'. see ncores.
    #ncores: integer. number of cores to use for parallel computation of 95%
    #   confidence intervals

    #updates `results` in memory and out/lm_out/results_specificq.csv (unless no_write is TRUE)
    #   if scale_q_by_area is FALSE, writes instead to out/lm_out/results.csv
    #returns list including:
    #   best_model: chosen model, as a formula
    #   best_model_ignoreint: glmnet only. same thing as best_model, but if
    #       regression through the origin was successful, this won't indicate it.
    #       used internally.
    #   best_model_object: the glmnet, lm, or segmented fit object
    #   prediction: vector of best-model predictions for each observation
    #   lm_data: a tibble of inputs and outputs. misnomer in case of glmnet or segmented
    #   metrics: nse, kge, rmse, mae, pbias for best model
    #   metrics_crossval: the same, but computed on out-of-sample folds during
    #       cross-validation
    #   fits: misnomer. just a plot of the marginal associations between each
    #       predictor gauge and the response. regression line included if model
    #       framework == 'lm'

    if(! framework %in% c('lm', 'glmnet', 'segmented')){
        stop('framework must be "lm", "glmnet", or "segmented"')
    }

    if(is.null(custom_gauges)){
        gage_ids <- donor_gauges[[neon_site]]
    } else {
        gage_ids <- custom_gauges
    }

    if(is.null(precomputed_df)){
        in_df <- assemble_q_df(neon_site = neon_site,
                               nearby_usgs_gages = gage_ids,
                               scale_q_by_area = scale_q_by_area)
    } else {
        in_df <- precomputed_df
    }

    formulaA <- 'discharge_log ~ `{paste(paste0(gage_ids, "_log"), collapse = "`+`")}` + season'
    formulaB <- 'discharge_log ~ `{paste(paste0(gage_ids, "_log"), collapse = "`*`")}` * season'
    if(! is.null(custom_formula)){
        formulaA <- formulaB <- custom_formula
    }

    if(framework == 'lm'){

        mods <- generate_nested_formulae(
            full_spec = as.formula(glue(formulaA)),
            d = in_df,
            max_interaction = 3
        )

        best <- lm_wrap(
            data = in_df,
            model_list = mods,
            unscale_q_by_area = scale_q_by_area,
            ...
        )

    } else if(framework == 'segmented'){

        mods <- generate_nested_formulae(
            full_spec = as.formula(glue(formulaA)),
            d = in_df,
            max_interaction = 3
        )

        best <- segmented_wrap(
            data = in_df,
            model_list = mods,
            unscale_q_by_area = scale_q_by_area,
            ...
        )

    } else { #glmnet

        best <- glmnet_wrap(
            data = in_df,
            full_spec = as.formula(glue(formulaB)),
            unscale_q_by_area = scale_q_by_area,
            ...
        )
    }

    if(! no_write){

        results <- plots_and_results(
            neon_site, best, results, in_df,
            unscale_q_by_area = scale_q_by_area,
            dummy_break = dummy_break,
            bootstrap_ci = bootstrap_ci,
            ncores = ncores
        )
        results[results$site_code == neon_site, 'method'] <- framework
        results <<- results

        if(scale_q_by_area){
            write_csv(results, 'out/lm_out/results_specificq.csv')
        } else {
            write_csv(results, 'out/lm_out/results.csv')
        }
    }

    return(best)
}

plots_and_results <- function(neon_site, best, results, in_df,
                              return_plot = FALSE, unscale_q_by_area = TRUE,
                              dummy_break = NULL, bootstrap_ci, ncores){

    #load corroborating usgs/ms site data
    sites_nearby = read_csv(glue('in/usgs_Q/{neon_site}.csv'))

    if(inherits(best$best_model_object, 'segmented')){
        sites_nearby <- rename_with(sites_nearby, ~sub('`?([0-9]{2,}(?:_log)?)`?', 'x\\1', .))
    }

    trms = rownames(attributes(terms(best$best_model))$factors)

    dummy_present <- 'dummy' %in% trms
    if(dummy_present){
        trms <- grep('dummy', trms, value = TRUE, invert = TRUE)
        if(is.null(dummy_break)) stop('dummy_break must be supplied if dummy variable present')
    }
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

    # q_eval = read_csv('in/NEON/neon_q_eval.csv') %>%
    #     filter(site == neon_site)
    #
    # q_eval = q_eval %>%
    #     group_by(site, year, month) %>%
    #     summarize(keep = all(final_qual %in% c('Tier1', 'Tier2')), #if issues occur at any point in a month, reject that whole month
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

    #predict Q for all datetimes with predictor data

    qall = left_join(sites_nearby, neon_q_auto, by = 'datetime')
        # left_join(select(neon_q_auto_qc, -site_code), by = 'datetime')

    if(dummy_present){
        qall$dummy <- 0
        qall$dummy[qall$datetime > dummy_break] <- 1
        qall$dummy <- as.factor(qall$dummy)
    }

    if(inherits(best$best_model_object, 'segmented')){

        qall <- predict(best$best_model_object,
                        newdata = select(qall, all_of(site_indeps_log), any_of('season')),
                        interval = 'predict') %>%
            as_tibble() %>%
            mutate(across(everything(), ~ifelse(. < 0, 0, .))) %>%
            mutate(across(everything(), inv_neglog)) %>%
            bind_cols(qall)

    } else if(inherits(best$best_model_object, 'glmnet')){

        betas <- best$best_model_object$beta@Dimnames[[1]]
        seasons_present <- str_extract(betas, 'season[1-4]') %>%
            unique() %>%
            na.omit() %>%
            str_split_i('season', 2)

        seasons_present <- c(seasons_present,
                             min(as.numeric(seasons_present)) - 1)

        newx_ <- select(qall, all_of(site_indeps_log), any_of('season')) %>%
            mutate(season = ifelse(! season %in% seasons_present, NA, season),
                   season = as.factor(season))
        valid_obs <- complete.cases(newx_)
        frm <- best$best_model_ignoreint
        has_intcpt <- '(Intercept)' %in% betas

        pred <- predict(
            best$best_model_object,
            s = best$best_model_object$lambda,
            newx = model.matrix(update(frm, NULL ~ .), newx_),
            type = 'response'
        )

        if(bootstrap_ci){
            cat('Bootstrapping 95% confidence intervals on', ncores, 'cores. Watch your RAM!')
            ci <- bootstrap_ci_glmnet(ncores, in_df, frm, best, has_intcpt, newx_)
        } else {
            ci <- list(ci_lwr = rep(NA_real_, length(pred)),
                       ci_upr = rep(NA_real_, length(pred)))
        }

        pred[pred < 0] <- 0
        ci$ci_lwr[ci$ci_lwr < 0] <- 0
        ci$ci_upr[ci$ci_upr < 0] <- 0

        qall$upr <- qall$lwr <- qall$fit <- NA_real_
        qall[valid_obs, 'fit'] <- inv_neglog(pred)
        qall[valid_obs, 'lwr'] <- inv_neglog(ci$ci_lwr)
        qall[valid_obs, 'upr'] <- inv_neglog(ci$ci_upr)

    } else { #lm

        qall <- predict(best$best_model_object,
                        newdata = qall,
                        interval = 'predict') %>%
            as_tibble() %>%
            mutate(across(everything(), ~ifelse(. < 0, 0, .))) %>%
            mutate(across(everything(), inv_neglog)) %>%
            bind_cols(qall)
    }

    if(unscale_q_by_area){
        wsa = filter(ms_areas, site_code == neon_site) %>% pull(ws_area_ha)
        qall = mutate(qall, fit = fit * wsa / 1000, lwr = lwr * wsa / 1000, upr = upr * wsa / 1000)
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

        #make diagnostic and "fit" plots

        png(glue('figs/lm_plots/diag/{neon_site}_diag_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
        defpar = par(mfrow=c(2, 2))
        if(inherits(best$best_model_object, 'segmented')){
            m <- best$best_model_object
            plot(m, conf.level=.95, is=TRUE, isV=FALSE, col=2, shade = TRUE, res=TRUE, res.col=4, pch=3)
            plot(fitted(m), m$residuals, xlab = 'Fitted values', ylab = 'Residuals')
            qqnorm(m$residuals, ylab = 'non-standardized residuals')
            qqline(m$residuals, lty = 3, col = 'gray40', lwd = 2)
        } else if(inherits(best$best_model_object, 'glmnet')){
            plot(best$best_model_object, xvar = 'norm', label = TRUE)
            plot(best$best_model_object, xvar = 'lambda', label = TRUE)
            plot(best$best_model_object, xvar = 'dev', label = TRUE)
        } else {
            plot(best$best_model_object)
        }
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

    asterisk <- if(inherits(best$best_model_object, 'segmented')) '*' else ''
    png(glue('figs/lm_plots/val/{neon_site}_obs_v_pred.png'), 6, 6, 'in', type = 'cairo', res = 300)
    plot(zz$Q_neon_field, zz$Q_predicted, xlab = 'NEON Field Discharge (L/s)',
         ylab = 'Predicted Discharge (L/s)',
         main = glue('Site: {neon_site}; KGE: {kge1}; KGE crossval{a}: {kge2}',
                     kge1 = round(best$metrics$kge, 2),
                     kge2 = round(best$metrics_crossval$kge, 2),
                     a = asterisk),
         xlim = axlim, ylim = axlim, xaxs = 'i', yaxs = 'i')
    abline(a = 0, b = 1, col = 'blue')
    legend('topleft', legend = '1:1', lty = 1, col = 'blue', bty = 'n')
    dev.off()

    out_data = filter(out_data, ! is.na(Q_predicted)) %>% select(-Q_neon_field)

    #save predictions as CSV
    out_data = left_join(out_data,
                         select(sites_nearby, datetime, all_of(site_indeps), any_of('season')),
                         by = 'datetime')
    if(dummy_present){
        out_data$dummy <- 0
        out_data$dummy[out_data$datetime > dummy_break] <- 1
        out_data$dummy <- as.factor(out_data$dummy)
    }

    write_csv(out_data, glue('out/lm_out/predictions/{neon_site}.csv'))

    #save fit data as CSV
    write_csv(best$lm_data, glue('out/lm_out/fit/{neon_site}.csv'))

    #save model summary
    sink(glue('out/lm_out/summary/{neon_site}.txt'))
    if(inherits(best$best_model_object, 'glmnet')){
        print(best$best_model_object)
        cat('\n---\n')
        print(summary(best$best_model_object))
        cat('\n---\n')
        print(best$best_model_object$beta)
    } else {
        print(summary(best$best_model_object))
    }
    sink()

    #return results frame, updated with this site's results
    results$bestmod[results$site_code == neon_site] = as.character(best$best_model)[3]
    results$nse[results$site_code == neon_site] = best$metrics$nse
    results$nse_cv[results$site_code == neon_site] = best$metrics_crossval$nse
    results$kge[results$site_code == neon_site] = best$metrics$kge
    results$kge_cv[results$site_code == neon_site] = best$metrics_crossval$kge
    results$pbias[results$site_code == neon_site] = best$metrics$pbias
    results$pbias_cv[results$site_code == neon_site] = best$metrics_crossval$pbias
    if(! inherits(best$best_model_object, 'glmnet')){
        results$adj_r_squared[results$site_code == neon_site] = summary(best$best_model_object)$adj.r.squared
    }

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

plots_and_results_daily_composite <- function(neon_site, best1, best2, results,
                                              unscale_q_by_area = TRUE, in_df1,
                                              in_df2, bootstrap_ci, ncores){

    #best2 and lm_df2 should represent the model with more terms included

    if(! setequal(class(best1), class(best2))) stop('best1 and best2 must have the same class(es)')

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

    #predict Q for all datetimes with predictor data

    qall = left_join(sites_nearby, neon_q_daily, by = 'date')
        # left_join(neon_q_daily_qc, by = 'date')

    if(inherits(best1$best_model_object, 'lm')){

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

        stop('nothing below this point has been tested to work with lm since glmnet took over at COMO')

    } else if(inherits(best1$best_model_object, 'glmnet')){

        qall1 <- qall2 <- qall

        #mod1
        betas <- best1$best_model_object$beta@Dimnames[[1]]
        seasons_present <- str_extract(betas, 'season[1-4]') %>%
            unique() %>%
            na.omit() %>%
            str_split_i('season', 2)

        seasons_present <- c(seasons_present,
                             min(as.numeric(seasons_present)) - 1)

        newx_ <- select(qall, all_of(site_indeps_log1), any_of('season')) %>%
            mutate(season = ifelse(! season %in% seasons_present, NA, season),
                   season = as.factor(season))
        valid_obs <- complete.cases(newx_)
        frm <- best1$best_model_ignoreint
        has_intcpt <- '(Intercept)' %in% betas

        pred1 <- predict(
            best1$best_model_object,
            s = best1$best_model_object$lambda,
            newx = model.matrix(update(frm, NULL ~ .), newx_),
            type = 'response'
        )

        if(bootstrap_ci){
            cat('Bootstrapping 95% confidence intervals on', ncores, 'cores. Watch your RAM!')
            ci <- bootstrap_ci_glmnet(ncores, in_df1, frm, best1, has_intcpt, newx_)
        } else {
            ci <- list(ci_lwr = rep(NA_real_, length(pred1)),
                       ci_upr = rep(NA_real_, length(pred1)))
        }

        pred1[pred1 < 0] <- 0
        ci$ci_lwr[ci$ci_lwr < 0] <- 0
        ci$ci_upr[ci$ci_upr < 0] <- 0

        qall1$upr <- qall1$lwr <- qall1$fit <- NA_real_
        qall1[valid_obs, 'fit'] <- inv_neglog(pred1)
        qall1[valid_obs, 'lwr'] <- inv_neglog(ci$ci_lwr)
        qall1[valid_obs, 'upr'] <- inv_neglog(ci$ci_upr)

        #mod2
        betas <- best2$best_model_object$beta@Dimnames[[1]]
        seasons_present <- str_extract(betas, 'season[1-4]') %>%
            unique() %>%
            na.omit() %>%
            str_split_i('season', 2)

        seasons_present <- c(seasons_present,
                             min(as.numeric(seasons_present)) - 1)

        newx_ <- select(qall, all_of(site_indeps_log2), any_of('season')) %>%
            mutate(season = ifelse(! season %in% seasons_present, NA, season),
                   season = as.factor(season))
        valid_obs <- complete.cases(newx_)
        frm <- best2$best_model_ignoreint
        has_intcpt <- '(Intercept)' %in% betas

        pred2 <- predict(
            best2$best_model_object,
            s = best2$best_model_object$lambda,
            newx = model.matrix(update(frm, NULL ~ .), newx_),
            type = 'response'
        )

        if(bootstrap_ci){
            cat('Bootstrapping 95% confidence intervals on', ncores, 'cores. Watch your RAM!')
            ci <- bootstrap_ci_glmnet(ncores, in_df2, frm, best2, has_intcpt, newx_)
        } else {
            ci <- list(ci_lwr = rep(NA_real_, length(pred2)),
                       ci_upr = rep(NA_real_, length(pred2)))
        }

        pred2[pred2 < 0] <- 0
        ci$ci_lwr[ci$ci_lwr < 0] <- 0
        ci$ci_upr[ci$ci_upr < 0] <- 0

        qall2$upr <- qall2$lwr <- qall2$fit <- NA_real_
        qall2[valid_obs, 'fit'] <- inv_neglog(pred2)
        qall2[valid_obs, 'lwr'] <- inv_neglog(ci$ci_lwr)
        qall2[valid_obs, 'upr'] <- inv_neglog(ci$ci_upr)

    } else stop('this function only specified for lm and glmnet models')

    #qall[complete.cases(newx_), c('fit', 'lwr', 'upr')]
    #qall$upr <- qall$lwr <- qall$fit <- NA_real_

    pred1 <- select(qall1, fit, lwr, upr)
    pred <- select(qall2, fit, lwr, upr)

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
               # Q_neon_continuous_filtered = discharge_daily_qc,
               Q_neon_continuous = discharge_daily) %>%
        filter(if_any(-date, ~cumsum(! is.na(.)) != 0)) %>%  #remove leading NA rows
        arrange(date)

    dg = dygraphs::dygraph(xts(x = select(out_data, -date, Q_pred_int_2.5, -Q_pred_int_97.5) %>% tail(5e5),
                               order.by = tail(out_data$date, 5e5))) %>%
        dyRangeSelector()

    saveWidget(dg, glue('figs/lm_plots/pred/{neon_site}_log.html'))

    png(glue('figs/lm_plots/diag/{neon_site}_diag_log_modA.png'), 6, 6, 'in', type = 'cairo', res = 300)
    defpar = par(mfrow=c(2,2))
    par(defpar)
    dev.off()

    #make diagnostic and fit plots
    png(glue('figs/lm_plots/diag/{neon_site}_diag_log_modA.png'), 6, 6, 'in', type = 'cairo', res = 300)
    defpar = par(mfrow=c(2,2))
    if(inherits(best1$best_model_object, 'glmnet')){
        plot(best1$best_model_object, xvar = 'norm', label = TRUE)
        plot(best1$best_model_object, xvar = 'lambda', label = TRUE)
        plot(best1$best_model_object, xvar = 'dev', label = TRUE)
    } else {
        plot(best1$best_model_object)
    }
    par(defpar)
    dev.off()

    png(glue('figs/lm_plots/diag/{neon_site}_diag_log_modB.png'), 6, 6, 'in', type = 'cairo', res = 300)
    defpar = par(mfrow=c(2,2))
    if(inherits(best2$best_model_object, 'glmnet')){
        plot(best2$best_model_object, xvar = 'norm', label = TRUE)
        plot(best2$best_model_object, xvar = 'lambda', label = TRUE)
        plot(best2$best_model_object, xvar = 'dev', label = TRUE)
    } else {
        plot(best2$best_model_object)
    }
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
         ylab = 'Predicted Discharge (L/s)',
         main = glue('Site: {neon_site}; KGE: {kge1}; KGE crossval*: {kge2}',
                     kge1 = round(kge, 2),
                     kge2 = NA),
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
    if(inherits(best1$best_model_object, 'glmnet')){
        print(best1$best_model_object)
        cat('\n---\n')
        print(summary(best1$best_model_object))
        cat('\n---\n')
        print(best1$best_model_object$beta)

        cat('\n\n--- MODEL B (2) ---\n\n')

        print(best2$best_model_object)
        cat('\n---\n')
        print(summary(best2$best_model_object))
        cat('\n---\n')
        print(best2$best_model_object$beta)
    } else {
        print(summary(best1$best_model_object))
        cat('\n---\n')
        print(summary(best2$best_model_object))
    }
    sink()

    #return results frame, updated with this site's results
    results$bestmod[results$site_code == neon_site] = paste(as.character(best1$best_model)[3],
                                                            as.character(best1$best_model)[3],
                                                            sep = '  &  ')
    results$nse[results$site_code == neon_site] = nse
    results$nse_cv[results$site_code == neon_site] = mean(best1$metrics_crossval$nse, best2$metrics_crossval$nse)
    results$adj_r_squared[results$site_code == neon_site] = NA
    # results$adj_r_squared[results$site_code == neon_site] = mean(summary(best1$best_model_object)$adj.r.squared,
    #                                                              summary(best2$best_model_object)$adj.r.squared)

    results$kge[results$site_code == neon_site] = kge
    results$kge_cv[results$site_code == neon_site] = mean(best1$metrics$kge, best2$metrics$kge)
    results$pbias[results$site_code == neon_site] = pbias
    results$pbias_cv[results$site_code == neon_site] = mean(best1$metrics$pbias, best2$metrics$pbias)

    if(unscale_q_by_area){
        write_csv(results, 'out/lm_out/results_specificq.csv')
    } else {
        write_csv(results, 'out/lm_out/results.csv')
    }

    return(results)
}

stratified_resample <- function(df, factor){

    indices <- split(seq_len(nrow(df)), df[[factor]])
    resampled_indices <- lapply(indices, sample, replace = TRUE)
    resampled_indices <- unlist(resampled_indices)

    return(df[resampled_indices,])
}

# random forest regression ####

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

    q_eval = read_csv('in/NEON/neon_q_eval.csv') %>%
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

# LSTM ####

run_lstm <- function(strategy, runset, ensemble = FALSE){

    #strategy: either 'generalist' or 'specialist'
    #runset: a numeric vector of run IDs
    #ensemble: logical. if true, runs finetune1 and finetune2 in sequence for
    #   each iteration of runset (so continue* and finetune* filesets must both
    #   be present in the corresponding config directories)

    #crude way to pass variables into python script from R
    r2pyenv_template <- new.env()
    r2pyenv_template$wdir <- getwd()
    r2pyenv_template$confdir <- confdir
    r2pyenv_template$rundir <- rundir
    r2pyenv_template$strategy <- strategy
    r2pyenv_template$ensemble <- ensemble
    r2pyenv_template$runset <- paste0('runs_', paste(range(runset), collapse = '-'))
    r2pyenv <<- as.list(r2pyenv_template)

    py_run_file('src/lstm_dungeon/run_lstms_local.py')
}

eval_lstms <- function(strategy, runset){

    #strategy: either 'generalist' or 'specialist'
    #runset: a numeric vector of run IDs

    phase <- ifelse(strategy == 'generalist', 'run', 'finetune')
    phase_token <- ifelse(phase == 'run', 'continue', 'finetune')

    #adjust test periods
    run_range <- paste(range(runset), collapse = '-')
    runset_parent_dir <- paste0('in/lstm_configs/runs_', run_range)
    run_dirs <- list.files(runset_parent_dir, pattern = 'run', full.names = FALSE)
    # for(rd in run_dirs){
    #     testrng_path <- file.path(runset_parent_dir, rd, 'test_ranges.csv')
    #     read_csv(testrng_path) %>%
    #         mutate(start_dates = as.character(holdout[1]),
    #                end_dates = as.character(holdout[2])) %>%
    #         write_csv(sub('\\.csv$', '_holdout.csv', testrng_path))
    # }

    #pickle (serialize) holdout ranges
    r2pyenv_template <- new.env()
    r2pyenv_template$confdir <- confdir
    r2pyenv_template$runset <- paste0('runs_', run_range)
    # r2pyenv <<- as.list(r2pyenv_template)
    # py_run_file('src/lstm_dungeon/pickle_test_periods.py')

    #make new configs that point to */test_ranges_holdout.pkl
    # for(rd in run_dirs){
    #
    #     rundir_ <- list.files('out/lstm_runs', pattern = rd, full.names = TRUE)
    #     if(! length(rundir_)) next
    #     cfg <- file.path(rundir_, 'config.yml')
    #     file.rename(cfg, paste0(cfg, '.bak'))
    #
    #     read_lines(paste0(cfg, '.bak')) %>%
    #         str_replace('test_ranges\\.pkl', 'test_ranges_holdout.pkl') %>%
    #         write_lines(cfg)
    # }

    #re-evaluate
    r2pyenv_template$wdir <- getwd()
    r2pyenv_template$runrange <- as.integer(runset)
    r2pyenv_template$phase <- phase
    r2pyenv <<- as.list(r2pyenv_template)
    py_run_file('src/lstm_dungeon/re-evaluate_models.py')

    #compute metrics (NauralHydrology 1.3.0 doesn't write to test_metrics.csv?)
    metrics <- list()
    for(rd in run_dirs){

        sites <- read_lines(file.path(runset_parent_dir, rd, 'test.txt'))
        sites <- neon_areas %>%
            filter(site_code %in% substr(sites, 1, 4))

        ptn <- if(strategy == 'specialist') sub('run', 'finetune', rd) else rd
        rundir_ <- list.files('out/lstm_runs', pattern = ptn)

        if(! length(rundir_)){
            message('no results for ', ptn, '; moving on')
            next
        }

        epoch_dir <- list.files(glue('out/lstm_runs/{r}/test', r = rundir_),
                                full.names = TRUE)
        if(length(epoch_dir) > 1) stop('handle this')

        lstm_out <- reticulate::py_load_object(file.path(epoch_dir, 'test_results.p'))
        metrics[[rd]] <- mapply(function(x, site){

            area <- sites$ws_area_ha[sites$site_code == site]

            pred <- x$`1D`$xr$discharge_sim$to_pandas()
            pred <- tibble(date = as.Date(rownames(pred)),
                           discharge_sim = pred$`0`) %>%
                #      L/s            mm/d             L/m^3  m^2/ha mm/m   s/d     ha
                mutate(discharge_sim = discharge_sim * 1000 * 1e4 / 1000 / 86400 * area)

            obs <- x$`1D`$xr$discharge_obs$to_pandas()
            obs <- tibble(date = as.Date(rownames(obs)),
                          discharge_obs = obs$`0`) %>%
                mutate(discharge_obs = discharge_obs * 1000 * 1e4 / 1000 / 86400 * area)

            predobs <- left_join(pred, obs, by = 'date')

            tibble(#site_code = site,
                   nse = hydroGOF::NSE(predobs$discharge_sim, predobs$discharge_obs),
                   kge = hydroGOF::KGE(predobs$discharge_sim, predobs$discharge_obs),
                   pbias = hydroGOF::pbias(predobs$discharge_sim, predobs$discharge_obs))

        }, lstm_out, sub('_MANUALQ', '', names(lstm_out))) %>%
            t() %>%
            as.data.frame() %>%
            tibble::rownames_to_column('site_code') %>%
            mutate(site_code = sub('_MANUALQ', '', site_code))
    }

    # #set configs to point back to */test_ranges.pkl
    # for(rd in run_dirs){
    #
    #     rundir_ <- list.files('out/lstm_runs', pattern = rd, full.names = TRUE)
    #     if(! length(rundir_)) next
    #     cfg <- file.path(rundir_, 'config.yml')
    #     suppressWarnings(try(file.rename(paste0(cfg, '.bak'), cfg), silent = TRUE))
    # }

    return(metrics)
}

identify_best_models <- function(result_list, kge_thresh){

    #returns the best model for each site, so long as it's above kge_thresh

    filt <- Map(
        function(x, run){
            filter(x, kge >= kge_thresh) %>%
                as_tibble() %>%
                mutate(across(-site_code, as.numeric),
                       run = run) %>%
                relocate(run, .after = 'site_code')
        },
        result_list,
        names(result_list)
    )

    smry <- reduce(filt, function(x, y) bind_rows(x, y)) %>%
        # arrange(site_code)
        group_by(site_code) %>%
        filter(kge == max(kge)) %>%
        ungroup()

    return(smry)
}

build_ensemble_config <- function(sites, runid, param_search, ensemb_n = 30){

    r <- str_extract(runid, '[0-9]+')
    strategy <- names(which(sapply(param_search, function(x) r %in% x)))
    if(length(strategy) != 1) stop('!')

    ## build new config dir structure

    existing_runs <- list.files(confdir, pattern = '^run')

    last_runid <- sapply(existing_runs, function(x){
        as.numeric(str_match(x, 'runs_[0-9]+-([0-9]+)')[, 2])
    }, USE.NAMES = FALSE) %>%
        max()

    next_runid <- last_runid + 1
    run_seq <- seq(next_runid, next_runid + ensemb_n - 1)

    cfgnew_ <- paste0('runs_', paste(range(run_seq), collapse = '-'))
    cfgnew <- file.path(confdir, cfgnew_, paste0('run', run_seq))
    lapply(cfgnew, dir.create, recursive = TRUE)

    ## copy and modify config(s) from which this ensemble is derived...

    if(grepl('specialist', strategy)){

        runid2 <- runid
        r2 <- r

        cfg00_ <- paste0('runs_', paste(range(param_search[[strategy]]), collapse = '-'))
        cfg00 <- file.path(confdir, cfg00_, runid2)

        runid <- read_lines(file.path(cfg00, paste0('finetune', r, '.yml'))) %>%
            str_subset('base_run_dir') %>%
            str_extract('/(run[0-9]+)_[0-9_]+$', group = 1)
        r <- str_extract(runid, '[0-9]+')
    }

    ## handle finetune1 (continue) files and test files

    src_strat <- case_when(grepl('generalist', strategy) ~ strategy,
                           strategy == 'specialist' ~ 'generalist',
                           strategy == 'pgdl_specialist' ~ 'pgdl_generalist')

    cfg0_ <- paste0('runs_', paste(range(param_search[[src_strat]]), collapse = '-'))
    cfg0 <- file.path(confdir, cfg0_, runid)

    file.copy(file.path(cfg0, 'continue.txt'),
              file.path(cfgnew, 'continue.txt'))
    file.copy(file.path(cfg0, 'continue_train_ranges.csv'),
              file.path(cfgnew, 'continue_train_ranges.csv'))

    test_new <- read_lines(file.path(cfg0, 'test.txt')) %>%
        str_subset(paste(sites, collapse = '|'))
    lapply(cfgnew, function(x) write_lines(test_new, file.path(x, 'test.txt')))

    testrange_new <- read_csv(file.path(cfg0, 'test_ranges.csv')) %>%
        filter(basin_id %in% paste0(sites, '_MANUALQ'))
    lapply(cfgnew, function(x) write_csv(testrange_new, file.path(x, 'test_ranges.csv')))

    yml <- read_lines(file.path(cfg0, paste0('continue', r, '.yml')))
    yml <- sub('seed: [0-9]+', paste('seed:', sample(1:99999, 1)), yml)
    yml <- sub(cfg0_, cfgnew_, yml)
    yml <- sub('save_weights_every: 1', 'save_weights_every: 999', yml)
    lapply(run_seq, function(x){
        out <- sub(runid, paste0('run', x), yml)
        write_lines(out,
                    file.path(confdir, cfgnew_, paste0('run', x),
                              paste0('continue', x, '.yml')))
    })

    if(grepl('specialist', strategy)){

        ## handle finetune2 files

        file.copy(file.path(cfg00, 'finetune.txt'),
                  file.path(cfgnew, 'finetune.txt'))
        file.copy(file.path(cfg00, 'finetune_train_ranges.csv'),
                  file.path(cfgnew, 'finetune_train_ranges.csv'))

        yml <- read_lines(file.path(cfg00, paste0('finetune', r2, '.yml')))
        yml <- sub('seed: [0-9]+', paste('seed:', sample(1:99999, 1)), yml)
        yml <- sub(cfg00_, cfgnew_, yml)
        lapply(run_seq, function(x){
            out <- sub(runid2, paste0('run', x), yml)
            sub('experiment_name: finetune[0-9]+',
                paste0('experiment_name: finetune', x),
                out) %>%
                write_lines(file.path(confdir, cfgnew_, paste0('run', x),
                                      paste0('finetune', x, '.yml')))
        })
    }

    #pickle (serialize) train and test ranges
    r2pyenv_template <- new.env()
    r2pyenv_template$confdir <- confdir
    r2pyenv_template$runset <- cfgnew_
    r2pyenv <<- as.list(r2pyenv_template)
    py_run_file('src/lstm_dungeon/pickle_traintest_periods.py')
}

# misc helpers required to generate figures/datasets ####

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

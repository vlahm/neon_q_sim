library(feather)
library(xts)
library(tidyverse)
library(lubridate)
library(glue)
library(ncdf4)
library(purrr)
library(reticulate)
library(dygraphs)
library(htmlwidgets)

dir.create('in/lstm_out', showWarnings = FALSE)
dir.create('in/lstm_out/pred', showWarnings = FALSE)
dir.create('in/lstm_out/val', showWarnings = FALSE)
dir.create('in/lstm_out/predictions', showWarnings = FALSE)
dir.create('in/lstm_out/fit', showWarnings = FALSE)

### 1b. adjust the test periods to include all manual Q (run once) ####

# # note: gotta rerun pickle_train_periods.py first

# runset_parent_dir = '/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs/runs_1718-1747'
# rundirs0 = list.files(runset_parent_dir, pattern = 'run', full.names = FALSE)
# for(rd in rundirs0){
#     testrng_path = file.path(runset_parent_dir, rd, 'test_ranges.csv')
#     read_csv(testrng_path) %>%
#         mutate(start_dates = '2010-01-01') %>%
#         write_csv(testrng_path)
# }

### 2. get LSTM stuff together ####

neon_site = 'TECR'

neon_q_manual = read_csv(glue('../imputation/data/neon_field_Q/{neon_site}.csv')) %>%
    mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
    rename(discharge_manual = discharge) %>%
    distinct(datetime, .keep_all = TRUE) %>%
    mutate(date = date(datetime))

runlist = paste(1718:1747, collapse = '|')
rundirs = list.files('../imputation/src/nh_methods/runs', pattern = '^run') %>%
    str_match(glue('run(?:{runlist})_[0-9_]+')) %>%
    {.[, 1]} %>%
    na.omit() #still has attrs

ws_area_ha = na.omit(site_data$ws_area_ha[site_data$site_code == neon_site])
ensemble_runs = list()
nses = rep(NA, length(rundirs))
holdout_nses = rep(NA, length(rundirs))
for(i in seq_along(rundirs)){

    r = rundirs[i]
    lstm_out = reticulate::py_load_object(glue('../imputation/src/nh_methods/runs/{r}/test/model_epoch020/test_results.p'))

    pred = lstm_out[[paste0(neon_site, '_MANQ_EXTRAP')]]$`1D`$xr$discharge_sim$to_pandas()
    pred = tibble(date = as.Date(rownames(pred)), Q = pred$`0`) %>%
        rename(discharge_sim = Q) %>%
        #      L/s            mm/d             L/m^3  m^2/ha mm/m   s/d     ha
        mutate(discharge_sim = discharge_sim * 1000 * 1e4 / 1000 / 86400 * ws_area_ha)

    qall = full_join(neon_q_manual, pred, by = 'date') %>%
        filter(! is.na(datetime))
    nses[i] = hydroGOF::NSE(sim = qall$discharge_sim, obs = qall$discharge_manual)
    holdout_nses[i] = hydroGOF::NSE(sim = qall$discharge_sim[qall$date > as.Date('2019-06-01')],
                                    obs = qall$discharge_manual[qall$date > as.Date('2019-06-01')])
    ensemble_runs[[i]] = pred
    # plot(qall$date, qall$discharge_manual, ylim = c(0,100))
    # lines(qall$date, qall$discharge_sim)
    # lines(qall$date[qall$date > as.Date('2020-01-01')], qall$discharge_sim[qall$date > as.Date('2020-01-01')], col='red')
    # print(paste(round(nses[i], 2), round(holdout_nses[i], 2)))
    # readLines(n=1)
}

nse_mean_overall = mean(nses)
nse_median_overall = median(nses)
nse_mean_holdout = mean(holdout_nses)
nse_median_holdout = median(holdout_nses)

ensemble_m = map(ensemble_runs, ~.$discharge_sim) %>%
    reduce(cbind) %>%
    as.matrix()
npreds = nrow(ensemble_m)
ensemble_out = tibble(
    date = ensemble_runs[[1]]$date,
    lower_bound_2.5 = rep(NA, npreds),
    mean_run = rep(NA, npreds),
    upper_bound_97.5 = rep(NA, npreds))
for(i in seq_len(npreds)){

    ensemble_out$mean_run[i] = mean(ensemble_m[i, ])
    # ens_sd = sd(ensemble_m[i, -1])
    ensemble_out$lower_bound_2.5[i] = quantile(ensemble_m[i, ], probs = 0.025)
    ensemble_out$upper_bound_97.5[i] = quantile(ensemble_m[i, ], probs = 0.975)
    # ensemble_out$lower_bound_2.5[i] = ensemble_out$mean_run[i] - (ens_sd * 1.96) #dist should be centered on mean
    # ensemble_out$upper_bound_97.5[i] = ensemble_out$mean_run[i] +
}

qall_tecr = full_join(neon_q_manual, ensemble_out, by = 'date') %>%
    filter(! is.na(datetime))
qall_forreals_tecr = full_join(neon_q_manual, ensemble_out, by = 'date') %>%
    arrange(date)
# nses[i] = hydroGOF::NSE(sim = qall_tecr$discharge_sim, obs = qall_tecr$discharge_manual)
ensemble_nse_holdout = hydroGOF::NSE(sim = qall_tecr$mean_run[qall_tecr$date > as.Date('2019-06-01')],
                                     obs = qall_tecr$discharge_manual[qall_tecr$date > as.Date('2019-06-01')])
ensemble_nse_overall = hydroGOF::NSE(sim = qall_tecr$mean_run,
                                     obs = qall_tecr$discharge_manual)


### 3. get neon sensor Q ####

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
    select(-keep, -year, -month, -site_code) %>%
    rename(discharge_auto_qc = discharge_auto)

nadts = is.na(qall_forreals_tecr$datetime)
qall_forreals_tecr$datetime[nadts] = as_datetime(qall_forreals_tecr$date[nadts])
qall_tecr_plot = qall_forreals_tecr %>%
    full_join(neon_q_auto, by = 'datetime') %>%
    left_join(neon_q_auto_qc, by = 'datetime') %>%
    mutate(site_code = neon_site) %>%
    filter(! is.na(discharge_manual) | substr(datetime, 15, 16) %in% c('00', '15', '30', '45'))

### 4. make output plots ####

out_data = qall_tecr_plot %>%
    select(datetime, Q_predicted = mean_run,
           # Q_pred_int_2.5 = lwr, Q_pred_int_97.5 = upr,
           Q_neon_field = discharge_manual,
           Q_neon_continuous_filtered = discharge_auto_qc,
           Q_neon_continuous_raw = discharge_auto) %>%
    arrange(datetime)

dg = dygraphs::dygraph(xts(x = select(out_data, -datetime) %>% tail(5e5),
                           order.by = tail(out_data$datetime, 5e5))) %>%
    dyRangeSelector()

saveWidget(dg, 'in/lstm_out/pred/{neon_site}_log.html')

#plot predictions versus field measurements. need to round field meas to Q interval
zz = out_data
field_dts = filter(zz, ! is.na(Q_neon_field)) %>% pull(datetime)

zz = zz %>%
    filter(datetime %in% field_dts) %>%
    group_by(datetime) %>%
    summarize(across(everything(), ~mean(., na.rm = TRUE))) %>%
    ungroup()

axlim = c(0, max(c(zz$Q_predicted, zz$Q_neon_field), na.rm = TRUE))

png('in/lstm_out/val/{neon_site}_obs_v_pred.png', 6, 6, 'in', type = 'cairo', res = 300)
plot(zz$Q_neon_field, zz$Q_predicted, xlab = 'NEON Field Discharge (L/s)',
     ylab = 'Predicted Discharge (L/s)', main = glue('Site: {neon_site}; NSE overall: {nse1}; NSE holdout: {nse2}',
                                                     nse1 = round(ensemble_nse_overall, 2),
                                                     nse2 = round(ensemble_nse_holdout, 2)),
     xlim = axlim, ylim = axlim, xaxs = 'i', yaxs = 'i')
abline(a = 0, b = 1, col = 'blue')
legend('topleft', legend = '1:1', lty = 1, col = 'blue', bty = 'n')
dev.off()

### 5. write model performance file into output dir ####

## combine lm, scaled lm, and lstm results.

lm_results = read_csv('out/q_lm_outdata/results.csv') %>%
    mutate(method = 'regression on nearby gauge discharge')
scaled_lm_results = read_csv('out/specific_q_lm_outdata/results.csv') %>%
    mutate(method = 'regression on nearby gauge discharge (specific)')
# lstm_results = read_csv('out/simulated_discharge_export/simulated_discharge_2022-07-11/model_skill.csv')
lstm_results = read_csv('out/lstm_model_skill.csv') %>%
    mutate(NSE = !!ensemble_nse_overall,
           NSE_holdout = !!ensemble_nse_holdout)

lm_results$mean_nse = rowMeans(lm_results[, c('nse_logq', 'nse_cv_logq')], na.rm = TRUE)
scaled_lm_results$mean_nse = rowMeans(scaled_lm_results[, c('nse_logq', 'nse_cv_logq')], na.rm = TRUE)

#use mean of NSE and CV-NSE to choose better mod between scaled and unscaled (by watershed area)
results = full_join(lm_results, scaled_lm_results, by = 'site_code')
x_better = results$mean_nse.x > results$mean_nse.y
x_better[is.na(x_better) & ! is.na(results$mean_nse.x)] = TRUE
x_better[is.na(x_better)] = FALSE
results = results %>%
    mutate(nse_logq = ifelse(x_better, nse_logq.x, nse_logq.y),
           nse_cv_logq = ifelse(x_better, nse_cv_logq.x, nse_cv_logq.y),
           bestmod_logq = ifelse(x_better, bestmod_logq.x, bestmod_logq.y),
           adj_r_squared = ifelse(x_better, adj_r_squared.x, adj_r_squared.y),
           mean_nse = ifelse(x_better, mean_nse.x, mean_nse.y),
           method = ifelse(x_better, method.x, method.y)) %>%
    select(-ends_with(c('.x', '.y')))

#use mean of NSE and CV-NSE (lm) and mean NSE across 30 runs (LSTM) to choose the better
results = full_join(results, lstm_results, by = 'site_code')
results$NSE_bestmethod = NA_real_
x_better = results$mean_nse > results$NSE
x_better[is.na(x_better) & ! is.na(results$NSE)] = FALSE
x_better[is.na(x_better)] = TRUE
results = results %>%
    mutate(NSE_bestmethod = ifelse(x_better, mean_nse, NSE),
           nse_logq = ifelse(x_better, nse_logq, NA_real_),
           nse_cv_logq = ifelse(x_better, nse_cv_logq, NA_real_),
           bestmod_logq = ifelse(x_better, bestmod_logq, NA_character_),
           adj_r_squared = ifelse(x_better, adj_r_squared, NA_real_),
           method = ifelse(x_better, method.x, method.y)) %>%
    select(-ends_with(c('.x', '.y')), -mean_nse, -NSE)

## clean up, add notes, etc.

# results$method[results$nse_cv_logq > results$NSE | is.na(results$NSE)] = 'regression on nearby gages'
# piecewise_used = grepl('U1', results$bestmod_logq)
# results$method[piecewise_used] = 'piecewise regression on nearby gauges'
results$note = NA_character_
# results$note[piecewise_used] = 'NSE not yet cross-validated.'
# results$NSE = sw(map2_dbl(results$nse_cv_logq, results$NSE, ~max(.x, .y, na.rm = TRUE)))
# results$NSE[is.infinite(results$NSE)] = results$nse_logq[is.infinite(results$NSE)]
results$method = sub('cherrypicked from generalists.*', 'generalist LSTM', results$method)
# zz = ! grepl('^discharge_log', results$bestmod_logq)
# results$bestmod_logq[zz] = paste('discharge_log ~', results$bestmod_logq[zz])
# lstm_best = grepl('LSTM', results$method)
# results$nse_cv_logq[lstm_best] = NA
# results$nse_logq[lstm_best] = NA
# results$bestmod_logq[lstm_best] = NA

results = results %>%
    select(site_code, method, NSE_bestmethod,
           NSE_regression_crossval = nse_cv_logq, NSE_regression = nse_logq,
           NSE_holdout,
           selected_regression_formula = bestmod_logq, adj_r_squared, note) %>%
    arrange(desc(NSE_bestmethod)) %>%
    mutate(across(starts_with('NSE'), ~round(., 3))) %>%
    mutate(NSE_bestmethod = ifelse(is.na(NSE_holdout), NA, NSE_bestmethod)) %>%
    rename(NSE_overall = NSE_bestmethod)

# results$note[results$site_code  == 'KING'] = 'Performance will be improved by updated Konza LTER data.'
results$note[results$site_code %in% c('BIGC')] = 'Performance may be improved by updated KREW data.'
# results$note[results$site_code  == 'GUIL'] = 'Over 3000 linear models tested. Good fit could be spurious.'
results$note[results$site_code  == 'POSE'] = 'Over 2000 linear models tested. Good fit could be spurious despite cross-validation.'

results = rename(results,
       NSE_LSTM_ensemble_overall = NSE_overall,
       NSE_LSTM_ensemble_holdout = NSE_holdout)
results = relocate(results, NSE_LSTM_ensemble_holdout, NSE_LSTM_ensemble_overall, .after = NSE_regression)

write_csv(results, file.path(outdir, 'model_performance.csv'))

# copy diagnostic, fit, and prediction plots to output dir for regression sites ####

#ended up using two methods, and this isn't going to change very often. just manually copying.

# by_regression = filter(results, grepl('regression', method)) %>% pull(site_code)
#
# fs = list.files('out/q_lm_plots/diag', full.names = TRUE)
# fs_sites = str_match(fs, '([A-Z]{4})_diag_log.png')[, 2]
# file.copy(fs[fs_sites %in% by_regression],
#           file.path(outdir, 'plots/diag/'),
#           recursive = TRUE)
#
# fs = list.files('out/q_lm_plots/fit', full.names = TRUE)
# fs_sites = str_match(fs, '([A-Z]{4})_fit_log.png')[, 2]
# file.copy(fs[fs_sites %in% by_regression],
#           file.path(outdir, 'plots/fit/'),
#           recursive = TRUE)
#
# fs = list.files('out/q_lm_plots/pred', full.names = TRUE)
# fs = grep('_files$', fs, invert = TRUE, value = TRUE)
# fs_sites = str_match(fs, '([A-Z]{4})_log')[, 2]
# file.copy(fs[fs_sites %in% by_regression],
#           file.path(outdir, 'plots/pred/'),
#           recursive = TRUE)
#
# fs = list.files('out/q_lm_plots/val', full.names = TRUE)
# fs_sites = str_match(fs, '([A-Z]{4})_obs_v_pred')[, 2]
# file.copy(fs[fs_sites %in% by_regression],
#           file.path(outdir, 'plots/val/'),
#           recursive = TRUE)

### 6. save outdata for lstm ####

qall_tecr %>%
    select(site_code, datetime, Q_neon_field = discharge_manual,
           Q_predicted = mean_run) %>%
    write_csv('in/lstm_out/fit/TECR.csv')

qall_tecr_plot %>%
    mutate(date = as.Date(datetime)) %>%
    select(date, Q_predicted = mean_run, Q_pred_int_2.5 = lower_bound_2.5,
           Q_pred_int_97.5 = upper_bound_97.5) %>%
           # Q_neon_continuous_filtered = discharge_auto_qc,
           # Q_neon_continuous_raw = discharge_auto) %>%
    filter(! is.na(Q_predicted)) %>%
    arrange(date) %>%
    write_csv('in/lstm_out/predictions/TECR.csv')


#field observations of Q by year for neon sites ####

# for(s in neon_sites){
#     print(s)
#     neon_q_manual = read_csv(glue('../imputation/data/neon_field_Q/{s}.csv')) %>%
#         mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
#         rename(discharge_manual = discharge) %>%
#         distinct(datetime, .keep_all = TRUE) %>%
#         arrange(datetime)
#
#     print(table(year(neon_q_manual$datetime)))
#     readLines(n=1)
# }

# count missing predictions for each neon site ####

fs = list.files(file.path(outdir, 'q_lm_outdata/predictions'), full.names = TRUE)

#by record, by year
outmissing = list()
for(i in seq_along(fs)){

    f = fs[i]
    neon_site = str_match(f, '([A-Z]+).csv$')[, 2]
    d = read_csv(f)

    if(neon_site != 'COMO'){
        outmissing[[i]] = d %>%
            group_by(year = year(datetime)) %>%
            summarize(first_prediction = first(datetime),
                      last_prediction = last(datetime),
                      n_timesteps = n(),
                      n_missing_predictions = sum(is.na(Q_predicted)),
                      pct_missing = n_missing_predictions / n_timesteps * 100) %>%
            mutate(site_code = neon_site) %>%
            relocate(site_code, .before = year)
    } else {
        outmissing[[i]] = d %>%
            group_by(year = year(date)) %>%
            summarize(first_prediction = as_datetime(first(date)),
                      last_prediction = as_datetime(last(date)),
                      n_timesteps = n(),
                      n_missing_predictions = sum(is.na(Q_predicted)),
                      pct_missing = n_missing_predictions / n_timesteps * 100) %>%
            mutate(site_code = neon_site) %>%
            relocate(site_code, .before = year)
    }
}

outmissing = map_dfr(outmissing, bind_rows) %>%
    arrange(site_code, year)

group_by(outmissing, site_code) %>%
    summarize(mean_pct_missing = mean(pct_missing)) %>%
    print(n = 100)

outmissing %>%
    filter(year > 2016) %>%
    group_by(site_code) %>%
    summarize(mean_pct_missing = mean(pct_missing)) %>%
    print(n = 100)

filter(outmissing, site_code == 'CUPE')


#by day, by year
outmissing_byday = list()
for(i in seq_along(fs)){

    f = fs[i]
    neon_site = str_match(f, '([A-Z]+).csv$')[, 2]
    d = read_csv(f)

    if(neon_site != 'COMO'){
        outmissing_byday[[i]] = d %>%
            group_by(date = as_date(datetime)) %>%
            summarize(Q_predicted = mean(Q_predicted, na.rm = TRUE)) %>%
            ungroup() %>%
            group_by(year = year(date)) %>%
            summarize(first_prediction = first(date),
                      last_prediction = last(date),
                      n_timesteps = n(),
                      n_missing_predictions = sum(is.na(Q_predicted)),
                      pct_missing = n_missing_predictions / n_timesteps * 100) %>%
            mutate(site_code = neon_site) %>%
            relocate(site_code, .before = year)
    } else {
        outmissing_byday[[i]] = d %>%
            group_by(year = year(date)) %>%
            summarize(first_prediction = as_datetime(first(date)),
                      last_prediction = as_datetime(last(date)),
                      n_timesteps = n(),
                      n_missing_predictions = sum(is.na(Q_predicted)),
                      pct_missing = n_missing_predictions / n_timesteps * 100) %>%
            mutate(site_code = neon_site) %>%
            relocate(site_code, .before = year)
    }
}


outmissing_byday = map_dfr(outmissing_byday, bind_rows) %>%
    arrange(site_code, year)

group_by(outmissing_byday, site_code) %>%
    summarize(mean_pct_missing = mean(pct_missing)) %>%
    print(n = 100)

outmissing_byday %>%
    filter(year > 2016) %>%
    group_by(site_code) %>%
    summarize(mean_pct_missing = mean(pct_missing)) %>%
    print(n = 100)

filter(outmissing_byday, site_code == 'CUPE')

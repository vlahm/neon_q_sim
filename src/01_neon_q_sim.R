# Mike Vlah
# vlahm13@gmail.com
# last data retrieval: 2023-01-31
# last edit: 2023-01-31

# library(segmented)
# library(feather)
# library(nhdplusTools)
library(dataRetrieval)
library(data.table)
# library(imputeTS)
library(lubridate)
library(glue)
library(dygraphs)
library(xts)
library(hydroGOF)
# library(errors)
# library(ncdf4)
# library(foreach)
# library(doParallel)
# library(sf)
library(ggplot2)
library(gridExtra)
library(htmlwidgets)
library(DAAG)
# library(forecast)
library(tidyverse)
library(macrosheds)

# reticulate::use_condaenv('nh2')
# xr <- reticulate::import("xarray")
# pd <- reticulate::import("pandas")
# np <- reticulate::import("numpy")

#TODO: switch on addition of USGS Q record for TECR if we want to predict that site better in the past (~1980)
#download.file(neon sites) for watershed areas
#what about sites tkin and tkot?

options(readr.show_progress = FALSE,
        readr.show_col_types = FALSE,
        timeout = 3000)

# set working directory wherever you'd like to build the file structure for this analysis.
# no other paths need to be modified.
setwd('~/git/macrosheds/papers/q_sim')

source('src/00_helpers.R')

## 1. retrieve input data ####

# NEON site metadata as reproduced by Rhea et al. 2023
# (primary source here: https://www.neonscience.org/field-sites/explore-field-sites)

if(! file.exists('in/neon_site_info.csv')){
    download.file('https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/data/contents/neon_site_info.csv',
                  destfile = 'in/neon_site_info.csv')
}

neon_areas <- read_csv('in/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    mutate(ws_area_ha = WSAreaKm2 * 100) %>%
    select(site_code = SiteID, ws_area_ha)

neon_sites <- neon_areas$site_code

# NEON discharge data (field measurements and continuous)

if(! length(list.files('in/neon_continuous_Q/'))){
    get_neon_inst_discharge(neon_sites)
}

if(! length(list.files('in/neon_field_Q/'))){
    get_neon_field_discharge(neon_sites)
}

# NEON discharge evaluation results from Rhea et al. 2023

if(! file.exists('in/neon_q_eval.csv')){
    download.file('https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/data/contents/neon_q_eval.csv',
                  destfile = 'in/neon_q_eval.csv')
}

q_eval <- read_csv('in/neon_q_eval.csv') %>%
    filter(site %in% neon_sites)

# MacroSheds data from H. J. Andrews and Niwot domains: used for donor gauges in lieu of USGS data

if(! dir.exists('in/hjandrews') | ! dir.exists('in/niwot')){
    macrosheds::ms_download_core_data(macrosheds_root = './in',
                                      domains = c('hjandrews', 'niwot'))
}

ms_q <- macrosheds::ms_load_product(macrosheds_root = './in',
                                    prodname = 'discharge',
                                    domains = c('hjandrews', 'niwot')) %>%
    filter(ms_status == 0)

# relevant MacroSheds site metadata

ms_areas <- macrosheds::ms_load_sites() %>%
    filter(site_type == 'stream_gauge',
           ! is.na(ws_area_ha)) %>%
    select(site_code, ws_area_ha)


## 2. q eval? - this is in the helpers now. make sure it's legit ####

# ms_sites <- unique(ms_q$site_code)

ms_q <- group_split(ms_q, site_code)

## clean NEON Q according to Rhea at al. 2023

check1 <- is.na(q_eval$regression_status) | q_eval$regression_status %in% c('good', 'fair')
check2 <- is.na(q_eval$drift_status) | q_eval$drift_status %in% c('likely_no_drift', 'not_assessed')
check3 <- is.na(q_eval$rating_curve_status) | q_eval$rating_curve_status %in% c('Tier1', 'Tier2')

q_eval$keep <- check1 & check2 & check3
q_eval <- q_eval %>%
    group_by(site, year, month) %>%
    # summarize(n = n()) %>%
    # filter(n > 1) %>%
    # print(n=100)
    summarize(keep = any(keep), #not being strict here, since so few instances of and plenty of visual subsetting anyway
              .groups = 'drop')

original_neon_q <- list()
for(i in seq_along(neon_sites)){

    s <- neon_sites[i]
    ni <- which(ms_sites == s)

    qqq <- ms_q[[ni]] %>%
        mutate(year = year(date),
               month = month(date))

    original_neon_q[[i]] <- select(qqq, date, site_code, discharge)

    qqq <- qqq %>%
        left_join(select(q_eval, site_code = site, year, month, keep),
                  by = c('site_code', 'year', 'month'))

    if(! s %in% c('GUIL', 'WLOU', 'BLDE', 'TOMB')){
        qqq = mutate(qqq, discharge = ifelse(keep, discharge, NA))
    }

    # if(s == 'BIGC') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2019-10-01'), discharge, NA))
    if(s == 'CARI') qqq$discharge[qqq$date > as.Date('2020-02-23') & qqq$date < as.Date('2020-03-31')] = NA
    if(s == 'CUPE') qqq = mutate(qqq, discharge = ifelse(date >= as.Date('2020-01-01'), discharge, NA))
    if(s == 'KING') qqq$discharge[qqq$date < as.Date('2018-12-01') | qqq$date > as.Date('2019-05-06')] = NA
    if(s == 'LEWI') qqq$discharge[qqq$date > as.Date('2020-03-05') & qqq$date < as.Date('2020-07-25')] = NA
    if(s == 'COMO') qqq$discharge[qqq$date > as.Date('2021-01-16') & qqq$date < as.Date('2021-03-06')] = NA
    # if(s == 'MAYF') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2019-07-19'), discharge, NA))
    if(s == 'MCDI') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2020-07-13'), discharge, NA))
    if(s == 'MCRA') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2019-05-06'), discharge, NA))
    # if(s == 'POSE') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2019-08-05'), discharge, NA))
    # if(s == 'GUIL') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2019-10-01'), discharge, NA))
    if(s == 'BLDE') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2020-10-01'), discharge, NA))
    # if(s == 'HOPB') qqq = mutate(qqq, discharge = ifelse(date <= as.Date('2020-12-28'), discharge, NA)) #it's actually P that's messed up here


    print(s)
    print(dygraphs::dygraph(xts(select(qqq, original_neon_q[[i]]$discharge, discharge), order.by = qqq$date)))
    readLines(n = 1)

    ms_q[[ni]] <- select(qqq, site_code, date, discharge)
}
names(original_neon_q) = sapply(original_neon_q, function(x) x$site_code[1])

## manually clean a few sites (leading and trailing NAs will get trimmed later)

zzx <- which(ms_sites == 'GSCC01')
# zz = ms_q[[zzx]]
# dygraphs::dygraph(xts::xts(select(zz, discharge),
#                                  order.by = zz$date)) %>%
#     dygraphs::dyRangeSelector()
ms_q[[zzx]] <- filter(ms_q[[zzx]], date > as.Date('1990-01-01'))
zzx <- which(ms_sites == 'GSCC02')
ms_q[[zzx]] <- filter(ms_q[[zzx]], date > as.Date('1990-01-01'))
zzx <- which(ms_sites == 'GSCC03')
ms_q[[zzx]] <- filter(ms_q[[zzx]], date > as.Date('1990-01-01'))
zzx <- which(ms_sites == 'GSCC04')
ms_q[[zzx]] <- filter(ms_q[[zzx]], date > as.Date('1990-01-01'))
# zzx <- which(q_names == 'WS79')
# ms_q[[zzx]] <- filter(ms_q[[zzx]], date > as.Date('2000-01-01'))
zzx <- which(ms_sites == 'QS')
ms_q[[zzx]] <- filter(ms_q[[zzx]], date > as.Date('2001-01-01'))
zzx <- which(ms_sites == 'GFVN')
ms_q[[zzx]] <- filter(ms_q[[zzx]], date > as.Date('1990-01-01'))

## 2. reconstruct early NEON discharge via regression on donor gauges ####

formulaB <- 'discharge_log ~ `{paste(paste0(gagenums, "_log"), collapse = "`+`")}` + season'

results <- tibble(site_code = neon_sites, nse_logq = NA, nse_cv_logq = NA, bestmod_logq = NA,
                  adj_r_squared = NA)

# REDB ####
neon_site = 'REDB'; gagenums = '10172200'
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
lm_df = filter(lm_df, as.Date(datetime) != as.Date('2019-06-12')) #remove erroneous-looking outlier
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df, interactions = TRUE, through_origin = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

# HOPB ####
neon_site = 'HOPB'; gagenums = c('01174565')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df, interactions = TRUE, through_origin = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

# KING ####
neon_site = 'KING'; gagenums = c('06879650', '06879810', '06879100', '06878600');
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# lm_df$`06879650`[lm_df$`06879650` > 1200 & ! is.na(lm_df$discharge)] = NA
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

library(caret)

# BLUE ####
neon_site = 'BLUE'; gagenums = '07332390'
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

# CUPE ####
neon_site = 'CUPE'; gagenums = c('50136400', '50138000', '50144000') #'50128907' no area available
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

# GUIL ####
neon_site = 'GUIL'; gagenums = c('50028000', '50024950', '50126150', '50026025')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE,
    min_points_per_param = 15,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

# # TECR (using KREW donor gauges) ####
# neon_site = 'TECR'; gagenums = c('11216400') #gagenums = 'T003'
# lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# mods = generate_nested_formulae(
#     full_spec = as.formula(glue(formulaB)),
#     d = lm_df,
#     interactions = TRUE)
# best = eval_model_set(data = lm_df, model_list = mods)
# results = plots_and_results(neon_site, best, lm_df, results)

#SYCA ####

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
best = eval_model_set(data = lm_df, model_list = mods)
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
# # saveWidget(dg, glue('../imputation/out/lm_plots/pred/{neon_site}_log.html'))
#
# png(glue('../imputation/out/lm_plots/fit/{neon_site}_fit_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
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
# # write_csv(lm_df, glue('../imputation/out/lm_out/by_site/{neon_site}.csv'))
#
# class(m_piecewise) = 'lm'
# png(glue('../imputation/out/lm_plots/diag/{neon_site}_diag_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
# defpar = par(mfrow=c(2,2))
# plot(m_piecewise)
# par(defpar)
# dev.off()

# WALK NO UP-TO-DATE REFERENCE Q ####
neon_site = 'WALK'; gagenums = c('03535000', '03535400', '03495405') #gagenums = c('east_fork', 'west_fork');
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

#TOMB ####
# neon_site = 'TOMB'; gagenums = '02469525'
neon_site = 'TOMB'; gagenums = c('02469761', '02469525')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

#ARIK ####
neon_site = 'ARIK'; gagenums = c('06827000', '06823000')#, '06821500' produces rank deficiencies
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# include_sensor_data = TRUE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

#MCDI ####
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
best = eval_model_set(data = lm_df, model_list = mods)
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
# # saveWidget(dg, glue('../imputation/out/lm_plots/pred/{neon_site}_log.html'))
#
# png(glue('../imputation/out/lm_plots/fit/{neon_site}_fit_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
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
# # write_csv(lm_df, glue('../imputation/out/lm_out/by_site/{neon_site}.csv'))
#
# class(m_piecewise) = 'lm'
# png(glue('../imputation/out/lm_plots/diag/{neon_site}_diag_log.png'), 6, 6, 'in', type = 'cairo', res = 300)
# defpar = par(mfrow=c(2,2))
# plot(m_piecewise)
# par(defpar)
# dev.off()


#LECO ####
neon_site = 'LECO'; gagenums = '03497300'
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

#LEWI ####
neon_site = 'LEWI'; gagenums = c('01636316', '01616100', '01636464')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

#PRIN MAYBE SOME KIND OF INTERP THING IS POSSIBLE? ####
# neon_site = 'PRIN'; gagenums = c('08044000', '08042950', '08042800')
neon_site = 'PRIN'; gagenums = c('08044000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

#POSE ####
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
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

#BLDE ####
neon_site = 'BLDE'; gagenums = c('06190540', '06188000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# include_sensor_daterange = c('2019-01-01', '2020-09-01'))
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

#BLWA ####
neon_site = 'BLWA'; gagenums = c('02466030', '02465000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

#MCRA ####
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
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

# #BIGC - needs KREW update (using daily mode) ####
# neon_site = 'BIGC'; gagenums = '11237500' #gagenums = c('P300', 'P301', 'P304') #gagenums = c('11238250')
# lm_df = assemble_q_lm_df_daily(neon_site = neon_site, nearby_usgs_gages = gagenums)
# mods = generate_nested_formulae(
#     full_spec = as.formula(glue(formulaB)),
#     d = lm_df,
#     interactions = TRUE)
# best = eval_model_set(data = lm_df, model_list = mods)
# results = plots_and_results_daily(neon_site, best, lm_df, results) #STOP: this func needs update

#MAYF ####
neon_site = 'MAYF'; gagenums = c('02465493', '02465292', '02424000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

#OKSR ####
neon_site = 'OKSR'; gagenums = c('15905100', '15908000', '15564879', '15875000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

#CARI ####
neon_site = 'CARI'; gagenums = c('15514000', '15511000', '15493400')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

#COMO (composite model) ####
formulaC = "discharge_log ~ `{paste(paste0(gagenums, \"_log\"), collapse = \"`+`\")}`"

neon_site = 'COMO'; gagenums = c('ALBION', 'SADDLE')
ms_d = Filter(function(x) x$site_code[1] %in% gagenums, ms_q) %>%
    reduce(bind_rows) %>%
    pivot_wider(names_from = site_code, values_from = discharge) %>%
    arrange(date)
lm_df1 = assemble_q_lm_df_daily(neon_site = neon_site, ms_Q_data = ms_d)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaC)),
    d = lm_df1,
    interactions = TRUE)
best1 = eval_model_set(data = lm_df1, model_list = mods)

neon_site = 'COMO'; gagenums = c('ALBION', 'SADDLE', 'MARTINELLI')
ms_d = Filter(function(x) x$site_code[1] %in% gagenums, ms_q) %>%
    reduce(bind_rows) %>%
    pivot_wider(names_from = site_code, values_from = discharge) %>%
    arrange(date)
lm_df2 = assemble_q_lm_df_daily(neon_site = neon_site, ms_Q_data = ms_d)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaC)),
    d = lm_df2,
    interactions = TRUE)
best2 = eval_model_set(data = lm_df2, model_list = mods)

results = plots_and_results_daily_composite(neon_site, best1, best2, lm_df1, lm_df2, results)

#FLNT ####
neon_site = 'FLNT'; gagenums = c('02355662', '02353000', '02356000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

# MART ####
neon_site = 'MART'; gagenums = c('14138870', '14123500', '14120000')
lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)#, include_pct_highvals = 10)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results = plots_and_results(neon_site, best, lm_df, results)

# # WLOU ####
# neon_site = 'WLOU'; gagenums = c('09026500', '09025300', '09027100', '09034900')
# lm_df = assemble_q_lm_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# mods = generate_nested_formulae(
#     full_spec = as.formula(glue(formulaB)),
#     d = lm_df,
#     interactions = TRUE)
# best = eval_model_set(data = lm_df, model_list = mods)
# results = plots_and_results(neon_site, best, lm_df, results)

## write results ####

write_csv(results, '../imputation/out/lm_out/results.csv')

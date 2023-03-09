# Mike Vlah
# vlahm13@gmail.com
# last data retrieval: 2023-01-31 (all files except hjandrews_q.txt, retrieved 2022-04-14)
# last edit: 2023-02-16

library(dataRetrieval)
library(data.table)
library(lubridate)
library(glue)
library(dygraphs)
library(xts)
library(hydroGOF)
library(ggplot2)
library(gridExtra)
library(neonUtilities)
library(htmlwidgets)
library(DAAG)
library(tidyverse)
library(macrosheds)

#TODO: switch on addition of USGS Q record for TECR if we want to predict that site better in the past (~1980)
#download.file(neon sites) for watershed areas
#what about sites tkin and tkot?

options(readr.show_progress = FALSE,
        readr.show_col_types = FALSE,
        timeout = 3000)

# set working directory wherever you'd like to build the file structure for this analysis.
# no other paths need to be modified. (except working directory in the plotting scripts)
setwd('~/git/macrosheds/papers/q_sim')

source('src/00_helpers.R')

## 1. retrieve input data ####

# NEON site metadata as reproduced by Rhea et al. 2023
# (primary source here: https://www.neonscience.org/field-sites/explore-field-sites)

if(! file.exists('in/neon_site_info.csv')){
    download.file('https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/data/contents/neon_site_info.csv',
                  destfile = 'in/neon_site_info.csv')
}

if(! file.exists('in/neon_site_info2.csv')){
    #filename changes with every update, so might have to modify URL below
    download.file('https://www.neonscience.org/sites/default/files/NEON_Field_Site_Metadata_20220412.csv',
                  destfile = 'in/neon_site_info2.csv')
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

# MacroSheds Q data from Niwot domain: used for donor gauges in lieu of USGS data

if(! dir.exists('in/niwot')){
    macrosheds::ms_download_core_data(macrosheds_root = './in',
                                      domains = 'niwot')
}

ms_q <- macrosheds::ms_load_product(macrosheds_root = './in',
                                    prodname = 'discharge',
                                    domains = 'niwot') %>%
    filter(ms_status == 0)

# H. J. Andrews Experimental Forest Q data: used for donor gauges in lieu of USGS data

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

donor_gauges <- yaml::read_yaml('in/donor_gauges.yml')

## 2. run setup ####

build_dir_structure()

# formulaA <- 'discharge ~ {paste(gagenums, collapse = "+")} + season'
formulaB <- 'discharge_log ~ `{paste(paste0(gagenums, "_log"), collapse = "`+`")}` + season'
formulaC = "discharge_log ~ `{paste(paste0(gagenums, \"_log\"), collapse = \"`+`\")}`"

results_lm <- results_lm_noscale <- tibble(
    site_code = neon_sites, nse = NA, nse_cv = NA, kge = NA, kge_cv = NA,
    pbias = NA, pbias_cv = NA, bestmod = NA, adj_r_squared = NA
)

## 3. linear regression (lm) on specific discharge ####

# REDB ####
neon_site = 'REDB'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# remove erroneous outlier (both sites on same stream reach. one could not report >100L/s while the other reports <10L/s)
lm_df = filter(lm_df, as.Date(datetime) != as.Date('2019-06-12'))
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df, interactions = TRUE, through_origin = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# HOPB ####
neon_site = 'HOPB'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df, interactions = TRUE, through_origin = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# KING ####
neon_site = 'KING'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# lm_df$`06879650`[lm_df$`06879650` > 1200 & ! is.na(lm_df$discharge)] = NA
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# BLUE ####
neon_site = 'BLUE'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# CUPE ####
neon_site = 'CUPE'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# GUIL ####
neon_site = 'GUIL'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE,
    min_points_per_param = 15,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# TECR [not yet viable. requires updated KREW donor gauges] ####
# neon_site = 'TECR'; gagenums = c('11216400') #gagenums = 'T003'
# lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# mods = generate_nested_formulae(
#     full_spec = as.formula(glue(formulaB)),
#     d = lm_df,
#     interactions = TRUE)
# best = eval_model_set(data = lm_df, model_list = mods)
# results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# SYCA ####

# neon_site = 'SYCA'; gagenums = c('09510200', '09499000') 09510150
neon_site = 'SYCA'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# include_sensor_daterange = c('2019-06-01', '2020-07-01'))
# ggplot(lm_df, aes(x=`09510200_log`, y = `discharge_log`)) + geom_point()

# remove erroneous outlier (sites separated by < 10 km on same river with no major intervening tribs.
#   very unlikely for one to report >100L/s while the other reports ~15L/s). This measurement also followed a long hiatus.
lm_df = filter(lm_df, ! as.Date(datetime) == as.Date('2021-07-26'))

mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE,
    min_points_per_param = 15,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

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
# results_lm = plots_and_results(neon_site, best, lm_df, results_lm)
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
# # results_lm$bestmod_logq[results_lm$site_code == neon_site] = gsub('x1', '09510200_log', paste(as.character(formula(m_piecewise))[c(2, 1, 3)], collapse = ' '))
# # results_lm$nse_logq[results_lm$site_code == neon_site] = nse_syca
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

# WALK ####
neon_site = 'WALK'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# TOMB ####
# neon_site = 'TOMB'; gagenums = '02469525'
neon_site = 'TOMB'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# ARIK ####
neon_site = 'ARIK'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# include_sensor_data = TRUE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# MCDI ####
neon_site = 'MCDI'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# include_pct_highvals = 5)
# include_sensor_daterange = c('2018-11-01', '2019-10-01'))
# ggplot(lm_df, aes(x=`06888500_log`, y = `discharge_log`)) + geom_point()
# ggplot(lm_df, aes(x=`06879650_log`, y = `discharge_log`)) + geom_point()
# lm_df = filter(lm_df, ! datetime %in% ymd_hms(c('2018-08-27 14:19:00', '2017-07-05 13:19:00', '2018-09-26 15:30:00')))
lm_df = filter(lm_df, ! datetime %in% ymd_hms(c('2018-08-27 14:19:00', '2022-03-23 19:00:00')))

mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

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
# results_lm = plots_and_results(neon_site, best, lm_df, results_lm)
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
# # results_lm$bestmod_logq[results_lm$site_code == neon_site] = gsub('x1', '06888500_log', paste(as.character(formula(m_piecewise))[c(2, 1, 3)], collapse = ' '))
# # results_lm$nse_logq[results_lm$site_code == neon_site] = nse_mcdi
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


# LECO ####
neon_site = 'LECO'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# LEWI ####
neon_site = 'LEWI'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# PRIN [very poor. maybe interpolation would help] ####
# neon_site = 'PRIN'; gagenums = c('08044000', '08042950', '08042800')
neon_site = 'PRIN'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# POSE ####
neon_site = 'POSE'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# include_pct_highvals = 10)
# include_sensor_daterange = c('2017-01-25', '2018-08-01'))
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE,
    min_points_per_param = 15,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# BLDE ####
neon_site = 'BLDE'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# include_sensor_daterange = c('2019-01-01', '2020-09-01'))
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# BLWA ####
neon_site = 'BLWA'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# MCRA ####

neon_site = 'MCRA'; gagenums = donor_gauges[[neon_site]]
hja = read_csv('in/hjandrews_q.txt') %>%
    filter(SITECODE %in% gagenums) %>%
    select(datetime = DATE_TIME, SITECODE, INST_Q) %>%
    mutate(INST_Q = INST_Q * 28.317) %>% #cfs to L/s
    pivot_wider(names_from = SITECODE, values_from = INST_Q) %>%
    arrange(datetime)
lm_df = assemble_q_df(neon_site = neon_site, ms_Q_data = hja)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# BIGC [not yet viable. requires updated KREW donor gauges] ####
# neon_site = 'BIGC'; gagenums = '11237500' #gagenums = c('P300', 'P301', 'P304') #gagenums = c('11238250')
# lm_df = assemble_q_df_daily(neon_site = neon_site, nearby_usgs_gages = gagenums)
# mods = generate_nested_formulae(
#     full_spec = as.formula(glue(formulaB)),
#     d = lm_df,
#     interactions = TRUE)
# best = eval_model_set(data = lm_df, model_list = mods)
# results_lm = plots_and_results_daily(neon_site, best, lm_df, results_lm) #STOP: this func needs update

# MAYF ####
neon_site = 'MAYF'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# OKSR ####
neon_site = 'OKSR'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    interactions = TRUE,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# CARI ####
neon_site = 'CARI'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaC)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# COMO [composite prediction] ####

neon_site = 'COMO'; gagenums = donor_gauges[[neon_site]][1:2]
ms_d = ms_q %>%
    filter(site_code %in% gagenums) %>%
    mutate(date = as.Date(datetime)) %>%
    select(date, site_code, val) %>%
    pivot_wider(names_from = site_code, values_from = val) %>%
    arrange(date)
lm_df1 = assemble_q_df_daily(neon_site = neon_site, ms_Q_data = ms_d)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaC)),
    d = lm_df1,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best1 = eval_model_set(data = lm_df1, model_list = mods)

neon_site = 'COMO'; gagenums = donor_gauges[[neon_site]]
ms_d = ms_q %>%
    filter(site_code %in% gagenums) %>%
    mutate(date = as.Date(datetime)) %>%
    select(date, site_code, val) %>%
    pivot_wider(names_from = site_code, values_from = val) %>%
    arrange(date)
lm_df2 = assemble_q_df_daily(neon_site = neon_site, ms_Q_data = ms_d)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaC)),
    d = lm_df2,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best2 = eval_model_set(data = lm_df2, model_list = mods)

results_lm = plots_and_results_daily_composite(neon_site, best1, best2, lm_df1, lm_df2, results_lm)

# FLNT ####
neon_site = 'FLNT'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# MART ####
neon_site = 'MART'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)#, include_pct_highvals = 10)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods)
results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

# WLOU [not yet viable. poor prediction and donor gauges missing seasons] ####
# neon_site = 'WLOU'; gagenums = c('09026500', '09025300', '09027100', '09034900')
# lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums)
# mods = generate_nested_formulae(
#     full_spec = as.formula(glue(formulaB)),
#     d = lm_df,
#     interactions = TRUE)
# best = eval_model_set(data = lm_df, model_list = mods)
# results_lm = plots_and_results(neon_site, best, lm_df, results_lm)

## 4. linear regression (lm) on un-scaled discharge ####

rename_dir_structure()
build_dir_structure()
# results_lm_noscale <- read_csv('out/lm_out/results.csv')

# REDB
neon_site = 'REDB'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
# remove erroneous outlier (both sites on same stream reach. one could not report >100L/s while the other reports <10L/s)
lm_df = filter(lm_df, as.Date(datetime) != as.Date('2019-06-12'))
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df, interactions = TRUE, through_origin = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# HOPB
neon_site = 'HOPB'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df, interactions = TRUE, through_origin = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# KING
neon_site = 'KING'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
# lm_df$`06879650`[lm_df$`06879650` > 1200 & ! is.na(lm_df$discharge)] = NA
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# BLUE
neon_site = 'BLUE'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# CUPE
neon_site = 'CUPE'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# GUIL
neon_site = 'GUIL'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE,
    min_points_per_param = 15,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# TECR [not yet viable. requires updated KREW donor gauges]
# neon_site = 'TECR'; gagenums = c('11216400') #gagenums = 'T003'
# lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
# mods = generate_nested_formulae(
#     full_spec = as.formula(glue(formulaB)),
#     d = lm_df,
#     interactions = TRUE)
# best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
# results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# SYCA

# neon_site = 'SYCA'; gagenums = c('09510200', '09499000') 09510150
neon_site = 'SYCA'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
# include_sensor_daterange = c('2019-06-01', '2020-07-01'))
# ggplot(lm_df, aes(x=`09510200_log`, y = `discharge_log`)) + geom_point()

# remove erroneous outlier (sites separated by < 10 km on same river with no major intervening tribs.
#   very unlikely for one to report >100L/s while the other reports ~15L/s). This measurement also followed a long hiatus.
lm_df = filter(lm_df, ! as.Date(datetime) == as.Date('2021-07-26'))

mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE,
    min_points_per_param = 15,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# WALK
neon_site = 'WALK'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# TOMB
# neon_site = 'TOMB'; gagenums = '02469525'
neon_site = 'TOMB'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# ARIK
neon_site = 'ARIK'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
# include_sensor_data = TRUE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# MCDI
neon_site = 'MCDI'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
# include_pct_highvals = 5)
# include_sensor_daterange = c('2018-11-01', '2019-10-01'))
# ggplot(lm_df, aes(x=`06888500_log`, y = `discharge_log`)) + geom_point()
# ggplot(lm_df, aes(x=`06879650_log`, y = `discharge_log`)) + geom_point()
# lm_df = filter(lm_df, ! datetime %in% ymd_hms(c('2018-08-27 14:19:00', '2017-07-05 13:19:00', '2018-09-26 15:30:00')))
lm_df = filter(lm_df, ! datetime %in% ymd_hms(c('2018-08-27 14:19:00', '2022-03-23 19:00:00')))

mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# LECO
neon_site = 'LECO'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# LEWI
neon_site = 'LEWI'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# PRIN [very poor. maybe interpolation would help]
# neon_site = 'PRIN'; gagenums = c('08044000', '08042950', '08042800')
neon_site = 'PRIN'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# POSE
neon_site = 'POSE'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
# include_pct_highvals = 10)
# include_sensor_daterange = c('2017-01-25', '2018-08-01'))
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    interactions = TRUE,
    min_points_per_param = 15,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# BLDE
neon_site = 'BLDE'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
# include_sensor_daterange = c('2019-01-01', '2020-09-01'))
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# BLWA
neon_site = 'BLWA'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# MCRA

neon_site = 'MCRA'; gagenums = donor_gauges[[neon_site]]
hja = read_csv('in/hjandrews_q.txt') %>%
    filter(SITECODE %in% gagenums) %>%
    select(datetime = DATE_TIME, SITECODE, INST_Q) %>%
    mutate(INST_Q = INST_Q * 28.317) %>% #cfs to L/s
    pivot_wider(names_from = SITECODE, values_from = INST_Q) %>%
    arrange(datetime)
lm_df = assemble_q_df(neon_site = neon_site, ms_Q_data = hja, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# BIGC [not yet viable. requires updated KREW donor gauges]
# neon_site = 'BIGC'; gagenums = '11237500' #gagenums = c('P300', 'P301', 'P304') #gagenums = c('11238250')
# lm_df = assemble_q_df_daily(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
# mods = generate_nested_formulae(
#     full_spec = as.formula(glue(formulaB)),
#     d = lm_df,
#     interactions = TRUE)
# best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
# results_lm_noscale = plots_and_results_daily(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE) #STOP: this func needs update

# MAYF
neon_site = 'MAYF'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# OKSR
neon_site = 'OKSR'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    interactions = TRUE,
    max_interaction = 3)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# CARI
neon_site = 'CARI'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaC)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# COMO [composite prediction]

neon_site = 'COMO'; gagenums = donor_gauges[[neon_site]][1:2]
ms_d = ms_q %>%
    filter(site_code %in% gagenums) %>%
    mutate(date = as.Date(datetime)) %>%
    select(date, site_code, val) %>%
    pivot_wider(names_from = site_code, values_from = val) %>%
    arrange(date)
lm_df1 = assemble_q_df_daily(neon_site = neon_site, ms_Q_data = ms_d, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaC)),
    d = lm_df1,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best1 = eval_model_set(data = lm_df1, model_list = mods, unscale_q_by_area = FALSE)

neon_site = 'COMO'; gagenums = donor_gauges[[neon_site]]
ms_d = ms_q %>%
    filter(site_code %in% gagenums) %>%
    mutate(date = as.Date(datetime)) %>%
    select(date, site_code, val) %>%
    pivot_wider(names_from = site_code, values_from = val) %>%
    arrange(date)
lm_df2 = assemble_q_df_daily(neon_site = neon_site, ms_Q_data = ms_d, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaC)),
    d = lm_df2,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best2 = eval_model_set(data = lm_df2, model_list = mods, unscale_q_by_area = FALSE)

results_lm_noscale = plots_and_results_daily_composite(neon_site, best1, best2, lm_df1, lm_df2, results_lm_noscale, unscale_q_by_area = FALSE)

# FLNT
neon_site = 'FLNT'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# MART
neon_site = 'MART'; gagenums = donor_gauges[[neon_site]]
lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)#, include_pct_highvals = 10)
mods = generate_nested_formulae(
    full_spec = as.formula(glue(formulaB)),
    d = lm_df,
    min_points_per_param = 15,
    max_interaction = 3,
    interactions = TRUE)
best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

# WLOU [not yet viable. poor prediction and donor gauges missing seasons]
# neon_site = 'WLOU'; gagenums = c('09026500', '09025300', '09027100', '09034900')
# lm_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
# mods = generate_nested_formulae(
#     full_spec = as.formula(glue(formulaB)),
#     d = lm_df,
#     interactions = TRUE)
# best = eval_model_set(data = lm_df, model_list = mods, unscale_q_by_area = FALSE)
# results_lm_noscale = plots_and_results(neon_site, best, lm_df, results_lm_noscale, unscale_q_by_area = FALSE)

## 5. random forest regression (abandoned; poor predictions, esp. out of Q sample range) ####

# library(caret)
# library(ranger)

# REDB
neon_site = 'REDB'; gagenums = '10172200'
rf_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
rf_df = select(rf_df, -ends_with('_log')) %>% rename_with(~paste0('x', .), matches('^[0-9]+$'))
gagenums = paste0('x', gagenums)
best = eval_model_set_rf(data = rf_df)
plots_and_results_rf(neon_site, best, rf_df, results_rf)

# KING
neon_site = 'KING'; gagenums = c('06879650', '06879810', '06879100', '06878600');
rf_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
rf_df = select(rf_df, -ends_with('_log')) %>% rename_with(~paste0('x', .), matches('^[0-9]+$'))
gagenums = paste0('x', gagenums)
best = eval_model_set_rf(data = rf_df)
plots_and_results_rf(neon_site, best, rf_df, results_rf)

# GUIL
neon_site = 'GUIL'; gagenums = c('50028000', '50024950', '50126150', '50026025')
rf_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
rf_df = select(rf_df, -ends_with('_log')) %>% rename_with(~paste0('x', .), matches('^[0-9]+$'))
gagenums = paste0('x', gagenums)
best = eval_model_set_rf(data = rf_df)
plots_and_results_rf(neon_site, best, rf_df, results_rf)

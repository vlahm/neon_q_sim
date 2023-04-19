# Mike Vlah
# vlahm13@gmail.com
# last data retrieval dates given in comments below:
# last edit: 2023-04-17

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
library(glmnet)
library(segmented)
library(tidyverse)
library(macrosheds)
library(parallel)
library(doParallel)

options(readr.show_progress = FALSE,
        readr.show_col_types = FALSE,
        timeout = 3000)

#pre-bundled in/out data available at: [**]
if(! exists('ts_plot')) source('src/00_helpers.R')
if(! exists('ms_areas')) source('src/01_data_retrieval.R')

## 1. run setup ####

build_dir_structure()

results <- tibble(
    site_code = neon_sites, method = NA, nse = NA, nse_cv = NA, kge = NA,
    kge_cv = NA, pbias = NA, pbias_cv = NA, bestmod = NA, adj_r_squared = NA
)

## 2. OLS (lm), piecewise (segmented, or ridge (glmnet) regression on specific discharge ####

## "standard" scenarios

regress(neon_site = 'REDB', framework = 'lm')
regress(neon_site = 'HOPB', framework = 'lm')
regress(neon_site = 'BLUE', framework = 'lm')
regress(neon_site = 'KING', framework = 'glmnet')
regress(neon_site = 'CUPE', framework = 'glmnet')
regress(neon_site = 'GUIL', framework = 'glmnet')
regress(neon_site = 'WALK', framework = 'glmnet')
regress(neon_site = 'TOMB', framework = 'glmnet')
regress(neon_site = 'ARIK', framework = 'glmnet')
regress(neon_site = 'MCDI', framework = 'glmnet')
regress(neon_site = 'LEWI', framework = 'glmnet')
regress(neon_site = 'PRIN', framework = 'glmnet', bootstrap_ci = FALSE)
regress(neon_site = 'POSE', framework = 'glmnet')
regress(neon_site = 'BLDE', framework = 'glmnet')
regress(neon_site = 'BLWA', framework = 'glmnet')
regress(neon_site = 'MAYF', framework = 'glmnet')
regress(neon_site = 'OKSR', framework = 'glmnet', bootstrap_ci = FALSE)
regress(neon_site = 'CARI', framework = 'glmnet')
regress(neon_site = 'FLNT', framework = 'glmnet')
regress(neon_site = 'MART', framework = 'glmnet')
# TECR not yet viable. requires updated KREW donor gauges
# BIGC not yet viable. requires updated KREW donor gauges
# WLOU not yet viable. poor prediction and donor gauges missing seasons

## other scenarios

#SYCA: erroneous outlier (sites separated by < 10 km on same river with no major intervening
#   tribs. very unlikely for one to report >100L/s while the other reports ~15L/s).
#   The questionable measurement also followed a long sampling hiatus.
syca_d <- assemble_q_df(neon_site = 'SYCA', nearby_usgs_gages = donor_gauges[['SYCA']]) %>%
    filter(! as.Date(datetime) == as.Date('2021-07-26'))
regress(neon_site = 'SYCA', framework = 'segmented', precomputed_df = syca_d)

#LECO: a wildfire in 2016 changed the intercept!
leco_d <- assemble_q_df(neon_site = 'LECO', nearby_usgs_gages = donor_gauges[['LECO']]) %>%
    mutate(dummy = ifelse(datetime > as.POSIXct('2016-11-27'), 1, 0),
           dummy = as.factor(dummy))
# plot(leco_d$`03497300_log`, leco_d$discharge_log, col = leco_d$dummy)
# plot(leco_d$datetime, leco_d$discharge_log, type = 'b', col = leco_d$dummy, pch = 20)
regress(neon_site = 'LECO', framework = 'lm', precomputed_df = leco_d,
        custom_formula = 'discharge_log ~ `{paste0(donor_gauges[["LECO"]], "_log")}` * season + dummy',
        dummy_break = as.Date('2016-11-27'))

#MCRA: uses Andrews Experimental Forest data
mcra_d_ <- read_csv('in/hjandrews_q.txt') %>%
    filter(SITECODE %in% donor_gauges[['MCRA']]) %>%
    select(datetime = DATE_TIME, SITECODE, INST_Q) %>%
    mutate(INST_Q = INST_Q * 28.317) %>% #cfs to L/s
    pivot_wider(names_from = SITECODE, values_from = INST_Q) %>%
    arrange(datetime)
mcra_d <- assemble_q_df(neon_site = 'MCRA', ms_Q_data = mcra_d_)
regress(neon_site = 'MCRA', framework = 'glmnet', precomputed_df = mcra_d)

#COMO: uses MacroSheds data and composite glmnet model
neon_site <- 'COMO'

gauge_ids <- donor_gauges[[neon_site]][1:2]
como_d_ <- ms_q %>%
    filter(site_code %in% gauge_ids) %>%
    mutate(date = as.Date(datetime)) %>%
    select(date, site_code, val) %>%
    pivot_wider(names_from = site_code, values_from = val) %>%
    arrange(date)
como_d1 <- assemble_q_df_daily(neon_site = neon_site, ms_Q_data = como_d_) %>%
    select(-contains('MARTINELLI'))
como_best1 <- regress(
    neon_site = neon_site, framework = 'glmnet', precomputed_df = como_d1,
    custom_formula = "discharge_log ~ `{paste(paste0(gauge_ids[1:2], \"_log\"), collapse = \"`*`\")}` * season",
    custom_gauges = gauge_ids, no_write = TRUE
)

gauge_ids <- donor_gauges[[neon_site]]
como_d_ <- ms_q %>%
    filter(site_code %in% gauge_ids) %>%
    mutate(date = as.Date(datetime)) %>%
    select(date, site_code, val) %>%
    pivot_wider(names_from = site_code, values_from = val) %>%
    arrange(date)
como_d2 <- assemble_q_df_daily(neon_site = neon_site, ms_Q_data = como_d_)
como_best2 <- regress(
    neon_site = neon_site, framework = 'glmnet',
    precomputed_df = como_d2, no_write = TRUE
)

results <- plots_and_results_daily_composite(
    neon_site = neon_site,
    best1 = como_best1,
    best2 = como_best2, #best2 must be the case with more gauges included
    results = results
)

results$method[results$site_code == 'COMO'] <- 'glmnet composite'


## 3. the same, but on absolute discharge ####

rename_dir_structure()
build_dir_structure()

results_specificq <- results
results <- tibble(
    site_code = neon_sites, method = NA, nse = NA, nse_cv = NA, kge = NA,
    kge_cv = NA, pbias = NA, pbias_cv = NA, bestmod = NA, adj_r_squared = NA
)

## "standard" scenarios

regress(neon_site = 'REDB', framework = 'lm', scale_q_by_area = FALSE)
regress(neon_site = 'HOPB', framework = 'lm', scale_q_by_area = FALSE)
regress(neon_site = 'BLUE', framework = 'lm', scale_q_by_area = FALSE)
regress(neon_site = 'KING', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'CUPE', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'GUIL', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'WALK', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'TOMB', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'ARIK', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'MCDI', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'LEWI', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'PRIN', framework = 'glmnet', scale_q_by_area = FALSE, bootstrap_ci = FALSE)
regress(neon_site = 'POSE', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'BLDE', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'BLWA', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'MAYF', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'OKSR', framework = 'glmnet', scale_q_by_area = FALSE, bootstrap_ci = FALSE)
regress(neon_site = 'CARI', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'FLNT', framework = 'glmnet', scale_q_by_area = FALSE)
regress(neon_site = 'MART', framework = 'glmnet', scale_q_by_area = FALSE)
# TECR not yet viable. requires updated KREW donor gauges
# BIGC not yet viable. requires updated KREW donor gauges
# WLOU not yet viable. poor prediction and donor gauges missing seasons

## other scenarios

#SYCA: erroneous outlier (sites separated by < 10 km on same river with no major intervening
#   tribs. very unlikely for one to report >100L/s while the other reports ~15L/s).
#   The questionable measurement also followed a long sampling hiatus.
syca_d <- assemble_q_df(neon_site = 'SYCA', scale_q_by_area = FALSE,
                        nearby_usgs_gages = donor_gauges[['SYCA']]) %>%
    filter(! as.Date(datetime) == as.Date('2021-07-26'))
regress(neon_site = 'SYCA', framework = 'segmented', precomputed_df = syca_d,
        scale_q_by_area = FALSE)

#LECO: a wildfire in 2016 changed the intercept!
leco_d <- assemble_q_df(neon_site = 'LECO', scale_q_by_area = FALSE,
                        nearby_usgs_gages = donor_gauges[['LECO']]) %>%
    mutate(dummy = ifelse(datetime > as.POSIXct('2016-11-27'), 1, 0),
           dummy = as.factor(dummy))
# plot(leco_d$`03497300_log`, leco_d$discharge_log, col = leco_d$dummy)
# plot(leco_d$datetime, leco_d$discharge_log, type = 'b', col = leco_d$dummy, pch = 20)
regress(neon_site = 'LECO', framework = 'lm', precomputed_df = leco_d,
        custom_formula = 'discharge_log ~ `{paste0(donor_gauges[["LECO"]], "_log")}` * season + dummy',
        dummy_break = as.Date('2016-11-27'), scale_q_by_area = FALSE)

#MCRA: uses Andrews Experimental Forest data
mcra_d_ <- read_csv('in/hjandrews_q.txt') %>%
    filter(SITECODE %in% donor_gauges[['MCRA']]) %>%
    select(datetime = DATE_TIME, SITECODE, INST_Q) %>%
    mutate(INST_Q = INST_Q * 28.317) %>% #cfs to L/s
    pivot_wider(names_from = SITECODE, values_from = INST_Q) %>%
    arrange(datetime)
mcra_d <- assemble_q_df(neon_site = 'MCRA', ms_Q_data = mcra_d_, scale_q_by_area = FALSE)
regress(neon_site = 'MCRA', framework = 'glmnet', precomputed_df = mcra_d,
        scale_q_by_area = FALSE)

#COMO: uses MacroSheds data and composite glmnet model
neon_site <- 'COMO'

gauge_ids <- donor_gauges[[neon_site]][1:2]
como_d_ <- ms_q %>%
    filter(site_code %in% gauge_ids) %>%
    mutate(date = as.Date(datetime)) %>%
    select(date, site_code, val) %>%
    pivot_wider(names_from = site_code, values_from = val) %>%
    arrange(date)
como_d1 <- assemble_q_df_daily(neon_site = neon_site, ms_Q_data = como_d_,
                               scale_q_by_area = FALSE) %>%
    select(-contains('MARTINELLI'))
como_best1 <- regress(
    neon_site = neon_site, framework = 'glmnet', precomputed_df = como_d1,
    custom_formula = "discharge_log ~ `{paste(paste0(gauge_ids[1:2], \"_log\"), collapse = \"`*`\")}` * season",
    custom_gauges = gauge_ids, no_write = TRUE, scale_q_by_area = FALSE
)

gauge_ids <- donor_gauges[[neon_site]]
como_d_ <- ms_q %>%
    filter(site_code %in% gauge_ids) %>%
    mutate(date = as.Date(datetime)) %>%
    select(date, site_code, val) %>%
    pivot_wider(names_from = site_code, values_from = val) %>%
    arrange(date)
como_d2 <- assemble_q_df_daily(neon_site = neon_site, ms_Q_data = como_d_,
                               scale_q_by_area = FALSE)
como_best2 <- regress(
    neon_site = neon_site, framework = 'glmnet',
    precomputed_df = como_d2, no_write = TRUE, scale_q_by_area = FALSE
)

results <- plots_and_results_daily_composite(
    neon_site = neon_site,
    best1 = como_best1,
    best2 = como_best2, #best2 must be the case with more gauges included
    results = results,
    unscale_q_by_area = FALSE
)

results$method[results$site_code == 'COMO'] <- 'glmnet composite'

## 4. random forest regression (abandoned; poor predictions, esp. out of Q sample range) ####

# library(caret)
# library(ranger)

# # REDB
# neon_site = 'REDB'; gagenums = '10172200'
# rf_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
# rf_df = select(rf_df, -ends_with('_log')) %>% rename_with(~paste0('x', .), matches('^[0-9]+$'))
# gagenums = paste0('x', gagenums)
# best = eval_model_set_rf(data = rf_df)
# plots_and_results_rf(neon_site, best, rf_df, results_rf)
#
# # KING
# neon_site = 'KING'; gagenums = c('06879650', '06879810', '06879100', '06878600');
# rf_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
# rf_df = select(rf_df, -ends_with('_log')) %>% rename_with(~paste0('x', .), matches('^[0-9]+$'))
# gagenums = paste0('x', gagenums)
# best = eval_model_set_rf(data = rf_df)
# plots_and_results_rf(neon_site, best, rf_df, results_rf)
#
# # GUIL
# neon_site = 'GUIL'; gagenums = c('50028000', '50024950', '50126150', '50026025')
# rf_df = assemble_q_df(neon_site = neon_site, nearby_usgs_gages = gagenums, scale_q_by_area = FALSE)
# rf_df = select(rf_df, -ends_with('_log')) %>% rename_with(~paste0('x', .), matches('^[0-9]+$'))
# gagenums = paste0('x', gagenums)
# best = eval_model_set_rf(data = rf_df)
# plots_and_results_rf(neon_site, best, rf_df, results_rf)

# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-02-26

library(tidyverse)
library(macrosheds)
library(glue)
library(reticulate)

reticulate::use_condaenv('nh2')
xr <- reticulate::import("xarray")
pd <- reticulate::import("pandas")
np <- reticulate::import("numpy")

options(readr.show_progress = FALSE,
        readr.show_col_types = FALSE,
        timeout = 3000)

# set working directory to same location as in 01_neon_q_sim.R
setwd('~/git/macrosheds/papers/q_sim')
# specify location where NeuralHydrology runs are stored
nh_dir <- '../../qa_experimentation/imputation/src/nh_methods/runs'
# specify run IDs for all tested generalist models
generalist_runids <- 1468:1520
# specify run IDs for replicates of best model for each specialist type
specialist_runids <- c(1548:1627, 1748:1937)
pgdl_runids <- 2028:2117

source('src/00_helpers.R')

## 1. setup ####

neon_sites <- read_csv('in/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    pull(SiteID)

#lm result tables
results_specq <- read_csv('out/lm_out/results_specificq.csv')
results_q <- read_csv('out/lm_out/results.csv')

#matrices for holding individual eval metrics
lstm_gen_nse <- lstm_gen_kge <- lstm_spec_nse <- lstm_spec_kge <- lstm_pgdl_nse <- lstm_pdgl_kge <-
    matrix(NA_real_,
           nrow = length(neon_sites),
           ncol = length(generalist_runids),
           dimnames = list(neon_sites, NULL))

#plottable results data.frame
plotd <- matrix(
    NA_real_, nrow = 27, ncol = 11,
    dimnames = list(
        NULL,
        c('site', 'nse_lm', 'nse_lm_scaled', 'nse_gen', 'nse_spec', 'nse_pgdl',
          'kge_lm', 'kge_lm_scaled', 'kge_gen', 'kge_spec', 'kge_pgdl')
    )
) %>% as_tibble()

plotd$site <- neon_sites

## 2. compile generalist results ####

for(i in seq_along(generalist_runids)){

    td <- try(locate_test_results(nh_dir, generalist_runids[i]), silent = TRUE)
    if(inherits(td, 'try-error') || is.null(td)) next

    xx = reticulate::py_load_object(td)

    for(s in neon_sites){
        try({
            lstm_gen_nse[rownames(lstm_gen_nse) == s, i] <- xx[[paste0(s, '_MANUALQ')]]$`1D`$discharge_NSE
            lstm_gen_kge[rownames(lstm_gen_nse) == s, i] <- xx[[paste0(s, '_MANUALQ')]]$`1D`$discharge_KGE
        }, silent = TRUE)
    }
}

## 2. assemble table of compiled results ####

#lm

#UPDATE THIS: read results.csv and results_noscale.csv. separate bars for each.
#line plot might be more effective? nah.
#order by mean kge across all 5 models

fs <- list.files('out/lm_out/fit', full.names = TRUE)
for(s in plotd$site){
    i <- which(plotd[, 'site'] == s)
    # predobs <- try(read_csv(glue('out/lm_out/fit/{s}.csv')), silent = TRUE)
    plotd$nse_lm_scaled[i] <- filter(results_specq, site_code == !!s) %>% pull(nse)
    plotd$nse_lm[i] <- filter(results_q, site_code == !!s) %>% pull(nse)
    plotd$kge_lm_scaled[i] <- filter(results_specq, site_code == !!s) %>% pull(kge)
    plotd$kge_lm[i] <- filter(results_q, site_code == !!s) %>% pull(kge)

    # plotd$nse_gen <-
    xx = reticulate::py_load_object(file.path(nh_dir, 'run1361_1907_032422/test/model_epoch030/test_results.p'))

    pred = xx[[paste0(s, '_GAPPED')]]$`1D`$xr$discharge_sim$to_pandas()

}



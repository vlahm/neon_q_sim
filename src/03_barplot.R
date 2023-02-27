# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-02-26

library(tidyverse)
library(macrosheds)
library(glue)

options(readr.show_progress = FALSE,
        readr.show_col_types = FALSE,
        timeout = 3000)

# set working directory to same location as in 01_neon_q_sim.R
setwd('~/git/macrosheds/papers/q_sim')

source('src/00_helpers.R')

## 1. setup ####

results_specq <- read_csv('out/lm_out/results_specificq.csv')
results_q <- read_csv('out/lm_out/results.csv')

plotd <- matrix(
    NA_real_, nrow = 27, ncol = 11,
    dimnames = list(
        NULL,
        c('site', 'nse_lm', 'nse_lm_scaled', 'nse_gen', 'nse_spec', 'nse_pgdl',
          'kge_lm', 'kge_lm_scaled', 'kge_gen', 'kge_spec', 'kge_pgdl')
    )
) %>% as_tibble()

plotd$site <- read_csv('in/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    pull(SiteID)

## 2. assemble model results ####

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
}



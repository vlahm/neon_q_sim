# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-02-28

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

pal <- viridis::viridis(5, end = 1)

neon_sites <- read_csv('in/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    pull(SiteID)

#initialize results data.frame
plotd <- matrix(
    NA_real_, nrow = 27, ncol = 5,
    dimnames = list(
        NULL,
        c('site', 'nse_lm', 'nse_lm_scaled', 'kge_lm', 'kge_lm_scaled')
        # c('site', 'nse_lm', 'nse_lm_scaled', 'nse_gen', 'nse_spec', 'nse_pgdl',
        #   'kge_lm', 'kge_lm_scaled', 'kge_gen', 'kge_spec', 'kge_pgdl')
    )
) %>% as_tibble()

plotd$site <- neon_sites

# 2. compile lm results ####

results_specq <- read_csv('out/lm_out_specQ/results_specificq.csv')
results_q <- read_csv('out/lm_out/results.csv')

for(s in plotd$site){

    i <- which(plotd[, 'site'] == s)

    plotd$nse_lm_scaled[i] <- filter(results_specq, site_code == !!s) %>% pull(nse)
    plotd$nse_lm[i] <- filter(results_q, site_code == !!s) %>% pull(nse)
    plotd$kge_lm_scaled[i] <- filter(results_specq, site_code == !!s) %>% pull(kge)
    plotd$kge_lm[i] <- filter(results_q, site_code == !!s) %>% pull(kge)
}

## 3. compile LSTM results ####

gen_res <- retrieve_test_results(generalist_runids)
spec_res <- retrieve_test_results(specialist_runids)
pgdl_res <- retrieve_test_results(pgdl_runids)

c('site', 'nse_lm', 'nse_lm_scaled', 'nse_gen', 'nse_spec', 'nse_pgdl',
  'kge_lm', 'kge_lm_scaled', 'kge_gen', 'kge_spec', 'kge_pgdl')

plotd <- left_join(plotd, reduce_results(gen_res$nse, 'nse_gen'))
plotd <- left_join(plotd, reduce_results(gen_res$kge, 'kge_gen'))
plotd <- left_join(plotd, reduce_results(spec_res$nse, 'nse_spec'))
plotd <- left_join(plotd, reduce_results(spec_res$kge, 'kge_spec'))
plotd <- left_join(plotd, reduce_results(pgdl_res$nse, 'nse_pgdl'))
plotd <- left_join(plotd, reduce_results(pgdl_res$kge, 'kge_pgdl'))

# 4. plot ####
plotd_nse <- select(plotd, -contains('kge'))
plotd_nse$rowmax <- apply(plotd_nse[, -1], 1, max, na.rm = TRUE)
plotd_nse <- arrange(plotd_nse, desc(rowmax))
plotd_m <- as.matrix(select(plotd_nse, -site, -rowmax))
plotd_m[! is.na(plotd_m) & plotd_m < -0.05] <- -0.05
plotd_m <- t(plotd_m)
rownames(plotd_m) <- c('linreg', 'linreg scaled', 'LSTM generalist', 'LSTM specialist', 'LSTM process-guided')

png(width = 8, height = 4, units = 'in', type = 'cairo', res = 300,
    filename = 'figs/fig2.png')
barplot(plotd_m, beside = TRUE, ylim = c(0, 1), names.arg = plotd$site,
        col = pal, las = 2, ylab = 'Nash-Sutcliffe Efficiency',
        legend.text = TRUE, border = FALSE,
        args.legend = list(x = 163, y=1.2, bty = 'n', cex = 0.9, border = FALSE,
                           xpd = NA, ncol = 3))
dev.off()

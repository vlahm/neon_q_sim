# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-02-28

library(tidyverse)
library(macrosheds)
library(glue)
library(reticulate)
library(data.table)

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
# specialist_runids <- c(1548:1627, 1748:1937) #search runs
specialist_runids <- c(2118:2247, 2293:2397)
# pgdl_runids <- 2028:2117 #search runs
pgdl_runids <- 2248:2292

source('src/00_helpers.R')

## 1. setup ####

pal <- c('black', rev(viridis::viridis(5, begin = 0.1, end = 1)))

if(! file.exists('in/neon_site_info.csv')){
    download.file('https://www.hydroshare.org/resource/03c52d47d66e40f4854da8397c7d9668/data/contents/neon_site_info.csv',
                  destfile = 'in/neon_site_info.csv')
}

neon_sites <- read_csv('in/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    pull(SiteID)

if(! length(list.files('in/neon_continuous_Q/'))){
    get_neon_inst_discharge(neon_sites)
}

if(! length(list.files('in/neon_field_Q/'))){
    get_neon_field_discharge(neon_sites)
}

#initialize results data.frame
plotd <- matrix(
    NA_real_, nrow = 27, ncol = 7,
    dimnames = list(
        NULL,
        c('site', 'nse_neon', 'kge_neon', 'nse_lm', 'nse_lm_scaled', 'kge_lm', 'kge_lm_scaled')
        # c('site', 'nse_lm', 'nse_lm_scaled', 'nse_gen', 'nse_spec', 'nse_pgdl',
        #   'kge_lm', 'kge_lm_scaled', 'kge_gen', 'kge_spec', 'kge_pgdl')
    )
) %>% as_tibble()

plotd$site <- neon_sites

## 2. compile lm results ####

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

# c('site', 'nse_lm', 'nse_lm_scaled', 'nse_gen', 'nse_spec', 'nse_pgdl',
#   'kge_lm', 'kge_lm_scaled', 'kge_gen', 'kge_spec', 'kge_pgdl')

plotd <- left_join(plotd, reduce_results(gen_res$nse, 'nse_gen', max, na.rm = TRUE))
plotd <- left_join(plotd, reduce_results(gen_res$kge, 'kge_gen', max, na.rm = TRUE))
plotd <- left_join(plotd, reduce_results(spec_res$nse, 'nse_spec', mean, na.rm = TRUE))
plotd <- left_join(plotd, reduce_results(spec_res$kge, 'kge_spec', mean, na.rm = TRUE))
plotd <- left_join(plotd, reduce_results(pgdl_res$nse, 'nse_pgdl', mean, na.rm = TRUE))
plotd <- left_join(plotd, reduce_results(pgdl_res$kge, 'kge_pgdl', mean, na.rm = TRUE))

## 4. compile NEON results ####

for(s in plotd$site){

    i <- which(plotd[, 'site'] == s)

    neon_q_manual <- read_csv(glue('in/neon_field_Q/{s}.csv')) %>%
        filter(! is.na(discharge)) %>%
        mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
        select(-site_code) %>%
        group_by(datetime) %>%
        summarize(discharge = mean(discharge, na.rm = TRUE)) %>%
        ungroup()

    neon_q_auto <- read_csv(glue('in/neon_continuous_Q/{s}.csv')) %>%
        filter(! is.na(discharge)) %>%
        select(-site_code)

    if(s == 'TOMB'){

        cmp <- approxjoin_datetime(
                x = mutate(neon_q_manual, site_code = s, var = 'discharge') %>%
                    rename(val = discharge),
                y = mutate(neon_q_auto, site_code = s, var = 'discharge') %>%
                    rename(val = discharge),
                rollmax = '12:00:00'
            ) %>%
            as_tibble() %>%
            select(datetime, discharge.man = val_x, discharge.aut = val_y)

    } else {

        cmp <- left_join(neon_q_manual, neon_q_auto, by = 'datetime', suffix = c('.man', '.aut'))
    }

    plotd[i, 'nse_neon'] <- hydroGOF::NSE(cmp$discharge.aut, cmp$discharge.man)
    plotd[i, 'kge_neon'] <- hydroGOF::KGE(cmp$discharge.aut, cmp$discharge.man)
}

## 5. figure 2 ####

# plotd <- select(plotd, site, nse_neon, kge_neon, nse_gen, kge_gen, nse_spec, kge_spec, nse_pgdl,
#                 kge_pgdl, nse_lm, kge_lm, nse_lm_scaled, kge_lm_scaled)
plotd <- select(plotd, site, nse_neon, kge_neon, nse_lm, kge_lm, nse_lm_scaled,
                kge_lm_scaled, nse_gen, kge_gen, nse_spec, kge_spec, nse_pgdl,
                kge_pgdl)

plotd_nse <- select(plotd, -contains('kge'))
plotd_nse$rowmax <- apply(select(plotd_nse, -site), 1, max, na.rm = TRUE)
plotd_nse <- arrange(plotd_nse, desc(nse_neon))
# plotd_nse <- arrange(plotd_nse, desc(rowmax))
plotd_m <- as.matrix(select(plotd_nse, -site, -rowmax))
plotd_m[! is.na(plotd_m) & plotd_m < -0.05] <- -0.05
plotd_m <- t(plotd_m)
rownames(plotd_m) <- c('Published', 'Linreg', 'Linreg scaled', 'LSTM generalist', 'LSTM specialist', 'LSTM process-guided')

png(width = 8, height = 4, units = 'in', type = 'cairo', res = 300,
    filename = 'figs/fig2_withneon.png')
plot(1:230, rep(0.5, 230), ylim = c(0, 1), ann = FALSE, axes = FALSE, col = 'transparent')
gray_bar_seq <- seq(0.2, 238, 17)
for(i in 1:14){
    ix <- gray_bar_seq[i] + i / 11
    polygon(c(ix, ix + 8.3, ix + 8.3, ix), c(-0.1, -0.1, 1.04, 1.04), col = 'gray85', border = FALSE, xpd = NA)
}
par(new = TRUE)
barplot(plotd_m, beside = TRUE, ylim = c(0, 1), names.arg = plotd_nse$site,
        col = pal, las = 2, ylab = 'Nash-Sutcliffe Efficiency',
        legend.text = TRUE, border = 'transparent',
        args.legend = list(x = 163, y=1.2, bty = 'n', cex = 0.9, border = FALSE,
                           xpd = NA, ncol = 3))
dev.off()

plotd_kge <- select(plotd, -contains('nse'))
plotd_kge$rowmax <- apply(select(plotd_kge, -site), 1, max, na.rm = TRUE)
plotd_kge <- arrange(plotd_kge, desc(rowmax))
plotd_m <- as.matrix(select(plotd_kge, -site, -rowmax))
plotd_m[! is.na(plotd_m) & plotd_m < -0.05] <- -0.05
plotd_m <- t(plotd_m)
rownames(plotd_m) <- c('Published', 'Linreg', 'Linreg scaled', 'LSTM generalist', 'LSTM specialist', 'LSTM process-guided')

png(width = 8, height = 4, units = 'in', type = 'cairo', res = 300,
    filename = 'figs/fig2_kge_withneon.png')
plot(1:230, rep(0.5, 230), ylim = c(0, 1), ann = FALSE, axes = FALSE, col = 'transparent')
gray_bar_seq <- seq(0.2, 238, 17)
for(i in 1:14){
    ix <- gray_bar_seq[i] + i / 11
    polygon(c(ix, ix + 8.3, ix + 8.3, ix), c(-0.36, -0.36, 1.04, 1.04), col = 'gray85', border = FALSE, xpd = NA)
}
par(new = TRUE)
barplot(plotd_m, beside = TRUE, ylim = c(0, 1), names.arg = plotd_kge$site,
        col = pal, las = 2, ylab = 'Kling-Gupta Efficiency',
        legend.text = TRUE, border = 'transparent',
        args.legend = list(x = 163, y=1.2, bty = 'n', cex = 0.9, border = FALSE,
                           xpd = NA, ncol = 3))
dev.off()

## 6. table 4 ####

plotd %>%
    summarize(across(starts_with(c('nse', 'kge')),
                     list(median = ~median(., na.rm = T),
                          min = ~min(., na.rm = T),
                          max = ~max(., na.rm = T),
                          n = ~sum(! is.na(.))),
                     .names = '{.col}_{.fn}')) %>%
    pivot_longer(everything(), names_to = c('score', 'model'), names_pattern = '^(nse|kge)_(.+)') %>%
    tidyr::extract('model', c('model', 'stat'), '(.*?)_([^_]+)$') %>%
    pivot_wider(names_from = c('stat', 'score'), values_from = 'value') %>%
    mutate(across(where(is.numeric), ~round(., 3))) %>%
    write_csv('out/score_table.csv')

## 7. stats ####


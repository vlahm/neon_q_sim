# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-02-28

library(tidyverse)
library(macrosheds)
library(glue)
library(reticulate)
library(data.table)
library(lubridate)
library(hydroGOF)

#python interaction required for loading lstm results.
#see step 2 in src/lstm_dungeon/README.txt
reticulate::use_condaenv('nh2')
xr <- reticulate::import("xarray")
pd <- reticulate::import("pandas")
np <- reticulate::import("numpy")

options(readr.show_progress = FALSE,
        readr.show_col_types = FALSE,
        timeout = 3000)

#pre-bundled in/out data available at: [**]
if(! exists('ts_plot')) source('src/00_helpers.R')
if(! exists('ms_areas')) source('src/01_data_retrieval.R')
if(! dir.exists('out/lm_out')) source('src/02_regression.R', local = new.env())
if(! dir.exists('in/lstm_data')) source('src/03_organize_camels_macrosheds_nhm.R', local = new.env())
if(! dir.exists('out/lstm_runs')) stop("you need to run src/04_run_lstms.R. It will take many days unless run on a cluster. Or use our bundled results.")

## 1. setup ####

nh_dir <- 'out/lstm_runs'
dir.create('out/neon_wateryear_assess', showWarnings = FALSE)

# # specify run IDs for all tested generalist models
# generalist_runids <- 1468:1520
# tecr_runids <- 1718:1747
# bigc_runids <-
# # specify run IDs for replicates of best model for each specialist type
# specialist_runids <- c(2293:2422)
# pgdl_runids <- 2248:2292

pal <- c('black', rev(viridis::viridis(5, begin = 0.2, end = 1)))

neon_sites <- read_csv('in/NEON/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    pull(SiteID)

#initialize results data.frame
plotd <- matrix(
    NA_real_, nrow = 27, ncol = 11,
    dimnames = list(
        NULL,
        c('site', 'nse_neon_min', 'nse_neon_overall', 'nse_neon_max',
          'kge_neon_min', 'kge_neon_overall', 'kge_neon_max',
          'nse_lm', 'nse_lm_scaled', 'kge_lm', 'kge_lm_scaled')
          # 'nse_gen', 'nse_spec', 'nse_gen_pgdl', 'nse_spec_pgdl',
          # 'kge_gen', 'kge_spec', 'kge_gen_pgdl', 'kge_spec_pgdl')
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

param_search_results <- read_csv('out/lstm_out/param_search_skill.csv') %>%
    mutate(strategy = case_when(strategy == 'gen' ~ 'generalist',
                                strategy == 'spec' ~ 'specialist',
                                strategy == 'pgdl_gen' ~ 'generalist_pgdl',
                                strategy == 'pgdl_spec' ~ 'specialist_pgdl')) %>%
    select(site = site_code, KGE = kge, NSE = nse, strategy) %>%
    pivot_wider(names_from = strategy, values_from = c(KGE, NSE)) %>%
    bind_rows(tibble(site = setdiff(plotd$site, .$site))) %>%
    rename_with(tolower) %>%
    arrange(site)

ensemble_results <- read_csv('out/lstm_out/results.csv') %>%
    select(-pbias, site = site_code) %>%
    pivot_wider(names_from = strategy, values_from = c(KGE, NSE)) %>%
    mutate(KGE_specialist_pgdl = NA_real_,
           NSE_specialist_pgdl = NA_real_) %>%
    rename_with(tolower)

plotd <- plotd %>%
    left_join(ensemble_results, by = 'site') %>%
    arrange(site)

for(col in colnames(ensemble_results)){
    col_ <- plotd[, col]
    plotd[is.na(col_), col] <- param_search_results[is.na(col_), col]
}

## 4. compile NEON results ####

for(s in plotd$site){

    i <- which(plotd[, 'site'] == s)

    neon_q_manual <- read_csv(glue('in/NEON/neon_field_Q/{s}.csv')) %>%
        filter(! is.na(discharge)) %>%
        mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
        select(-site_code) %>%
        group_by(datetime) %>%
        summarize(discharge = mean(discharge, na.rm = TRUE)) %>%
        ungroup()

    neon_q_auto <- read_csv(glue('in/NEON/neon_continuous_Q/{s}.csv')) %>%
        filter(! is.na(discharge)) %>%
        select(-site_code)

    # print(mode_interval_dt(neon_q_auto$datetime))

    if(s == 'TOMB'){

        cmp <- approxjoin_datetime(
                x = rename(neon_q_manual, val = discharge),
                y = rename(neon_q_auto, val = discharge),
                rollmax = '00:30:00'
            ) %>%
            as_tibble() %>%
            select(datetime, discharge.man = val_x, discharge.aut = val_y)

    } else {

        cmp <- approxjoin_datetime(
                x = rename(neon_q_manual, val = discharge),
                y = rename(neon_q_auto, val = discharge),
                rollmax = '00:02:30'
            ) %>%
            as_tibble() %>%
            select(datetime, discharge.man = val_x, discharge.aut = val_y)
    }

    plotd[i, 'nse_neon_overall'] <- NSE(cmp$discharge.aut, cmp$discharge.man)
    plotd[i, 'kge_neon_overall'] <- KGE(cmp$discharge.aut, cmp$discharge.man)

    cmp_splt <- cmp %>%
        mutate(wateryr = year(datetime),
               wateryr = ifelse(month(datetime) >= 10, wateryr + 1, wateryr)) %>%
        group_by(wateryr) %>%
        group_split()

    #needed to build composite Q
    suppressWarnings({
        purrr::reduce(cmp_splt, function(x, y){
            y_ <- tibble(wateryr = y$wateryr[1],
                         nse = NSE(y$discharge.aut, y$discharge.man),
                         kge = KGE(y$discharge.aut, y$discharge.man))
            # if(nrow(y) < 5) { y_$nse <- NA; y_$kge <- NA }
            bind_rows(x, y_)
        }, .init = tibble()) %>%
            write_csv(glue('out/neon_wateryear_assess/{s}.csv'))
    })

    if(s == 'OKSR'){
        cmp_splt <- purrr::keep(cmp_splt, ~length(na.omit(.$discharge.aut)) >= 2)
    } else {
        cmp_splt <- purrr::keep(cmp_splt, ~length(na.omit(.$discharge.aut)) >= 5)
    }

    nse_wy <- purrr::map(cmp_splt, ~NSE(.$discharge.aut, .$discharge.man))
    kge_wy <- purrr::map(cmp_splt, ~KGE(.$discharge.aut, .$discharge.man))

    plotd[i, 'nse_neon_min'] <- purrr::reduce(nse_wy, min)
    plotd[i, 'nse_neon_max'] <- purrr::reduce(nse_wy, max)
    plotd[i, 'kge_neon_min'] <- purrr::reduce(kge_wy, min)
    plotd[i, 'kge_neon_max'] <- purrr::reduce(kge_wy, max)
}

## 5. figure 2 ####

statd <- plotd

#fig Sx (same as fig 2, but showing NSE)

plotd <- plotd %>%
    mutate(nse_lm = pmax(nse_lm, nse_lm_scaled),
           kge_lm = pmax(kge_lm, kge_lm_scaled)) %>%
    select(site, nse_neon_min, nse_neon_overall, nse_neon_max,
           kge_neon_min, kge_neon_overall, kge_neon_max,
           nse_lm, kge_lm,
           nse_generalist, kge_generalist, nse_specialist, kge_specialist,
           nse_generalist_pgdl, kge_generalist_pgdl,
           nse_specialist_pgdl, kge_specialist_pgdl)

plotd_nse <- select(plotd, -contains('kge'))
plotd_nse$rowmax <- apply(select(plotd_nse, -site, -ends_with('max')), 1, max, na.rm = TRUE)
plotd_nse <- arrange(plotd_nse, desc(nse_neon_overall))
# plotd_nse <- arrange(plotd_nse, desc(rowmax))
plotd_m <- as.matrix(select(plotd_nse, -site, -ends_with(c('max', 'min'))))
plotd_m[! is.na(plotd_m) & plotd_m < -0.05] <- -0.05
plotd_m <- t(plotd_m)
plotd_nse <- mutate(plotd_nse, nse_neon_min = ifelse(nse_neon_min < -0.05, -0.05, nse_neon_min))
rownames(plotd_m) <- c('Published', 'Linreg', 'LSTM generalist',
                       'LSTM specialist', 'LSTM process-guided generalist', 'LSTM process-guided specialist')

png(width = 8, height = 4, units = 'in', type = 'cairo', res = 300,
    filename = 'figs/fig2_nse.png')
par(mar = c(5, 2, 2.5, 0), oma = c(0, 0, 0, 0))
plot(1:230, rep(0.5, 230), ylim = c(0, 1), ann = FALSE, axes = FALSE, col = 'transparent')
gray_bar_seq <- seq(0.2, 238, 17)
for(i in 1:14){
    ix <- gray_bar_seq[i] + i / 11
    polygon(c(ix, ix + 8.3, ix + 8.3, ix), c(-0.36, -0.36, 1.04, 1.04), col = 'gray85', border = FALSE, xpd = NA)
}
par(new = TRUE)
barplot(plotd_m, beside = TRUE, ylim = c(0, 1), names.arg = plotd_nse$site,
        col = pal, las = 2, ylab = '',
        legend.text = TRUE, border = 'transparent', yaxt = 'n',
        args.legend = list(x = 163, y=1.2, bty = 'n', cex = 0.9, border = FALSE,
                           xpd = NA, ncol = 3, text.width = c(30, 40, 50)))
segments(0.5, -0.045, 190, -0.045, xpd = NA, lwd = 0.1, lty = 2)
minmax_seq <- seq(1.55, 245, by = 7)
for(i in 1:27){
    points(minmax_seq[i], plotd_nse$nse_neon_max[i], col = 'black', bg = 'white',
           xpd = NA, pch = 24, cex = 0.4, lwd = 0.25)
    points(minmax_seq[i], plotd_nse$nse_neon_min[i], col = 'black', bg = 'white',
           xpd = NA, pch = 25, cex = 0.4, lwd = 0.25)
}
mtext('NEON stream/river site', 1, 3.8, font = 2)
mtext('Nash-Sutcliffe Efficiency', 2, 1, font = 2)
axis(2, seq(0, 1, 0.1), line = -0.9, tcl = -0.3, padj = 1)
dev.off()

#actual fig 2

plotd_kge <- select(plotd, -contains('nse'))
plotd_kge$rowmax <- apply(select(plotd_kge, -site, -ends_with('max')), 1, max, na.rm = TRUE)
plotd_kge <- arrange(plotd_kge, desc(kge_neon_overall))
plotd_m <- as.matrix(select(plotd_kge, -site, -ends_with(c('max', 'min'))))
plotd_m[! is.na(plotd_m) & plotd_m < -0.05] <- -0.05
plotd_m <- t(plotd_m)
plotd_kge <- mutate(plotd_kge, kge_neon_min = ifelse(kge_neon_min < -0.05, -0.05, kge_neon_min))
rownames(plotd_m) <- c('Published', 'Linreg', 'LSTM generalist',
                       'LSTM specialist', 'LSTM process-guided generalist', 'LSTM process-guided specialist')

png(width = 8, height = 4, units = 'in', type = 'cairo', res = 300,
    filename = 'figs/fig2_kge.png')
par(mar = c(5, 2, 2.5, 0), oma = c(0, 0, 0, 0))
plot(1:230, rep(0.5, 230), ylim = c(0, 1), ann = FALSE, axes = FALSE, col = 'transparent')
gray_bar_seq <- seq(0.2, 238, 17)
for(i in 1:14){
    ix <- gray_bar_seq[i] + i / 11
    polygon(c(ix, ix + 8.3, ix + 8.3, ix), c(-0.36, -0.36, 1.04, 1.04), col = 'gray85', border = FALSE, xpd = NA)
}
par(new = TRUE)
barplot(plotd_m, beside = TRUE, ylim = c(0, 1), names.arg = plotd_kge$site,
        col = pal, las = 2, ylab = '',
        legend.text = TRUE, border = 'transparent', yaxt = 'n',
        args.legend = list(x = 163, y=1.2, bty = 'n', cex = 0.9, border = FALSE,
                           xpd = NA, ncol = 3, text.width = c(30, 40, 50)))
segments(0.5, -0.045, 190, -0.045, xpd = NA, lwd = 0.1, lty = 2)
minmax_seq <- seq(1.55, 245, by = 7)
for(i in 1:27){
    points(minmax_seq[i], plotd_kge$kge_neon_max[i], col = 'black', bg = 'white',
           xpd = NA, pch = 24, cex = 0.4, lwd = 0.25)
    points(minmax_seq[i], plotd_kge$kge_neon_min[i], col = 'black', bg = 'white',
           xpd = NA, pch = 25, cex = 0.4, lwd = 0.25)
}
mtext('NEON stream/river site', 1, 3.8, font = 2)
mtext('Kling-Gupta Efficiency', 2, 1, font = 2)
axis(2, seq(0, 1, 0.1), line = -0.9, tcl = -0.3, padj = 1)
dev.off()

## 6. table 4 ####

statd %>%
    select(-ends_with(c('min', 'max'))) %>%
    summarize(across(starts_with(c('nse', 'kge')),
                     list(median = ~median(., na.rm = T),
                          mean = ~mean(., na.rm = T),
                          min = ~min(., na.rm = T),
                          max = ~max(., na.rm = T),
                          n = ~sum(! is.na(.))),
                     .names = '{.col}_{.fn}')) %>%
    pivot_longer(everything(), names_to = c('score', 'model'), names_pattern = '^(nse|kge)_(.+)') %>%
    tidyr::extract('model', c('model', 'stat'), '(.*?)_([^_]+)$') %>%
    pivot_wider(names_from = c('stat', 'score'), values_from = 'value') %>%
    mutate(across(where(is.numeric), ~round(., 3))) %>%
    select(-n_nse) %>%
    rename(n_sites = n_kge) %>%
    write_csv('out/score_table.csv')

## 7. stats ####

#number of field measurements per site
n_field_meas <- sapply(
    list.files('in/NEON/neon_field_Q', full.names = TRUE),
    function(x) length(na.omit(read_csv(x)$discharge))
) %>% sort()

range(n_field_meas)
mean(n_field_meas)

#best scores across models
maxd_ <- select(statd, -ends_with(c('neon_min', 'neon_max')))
maxd <- maxd_ %>%
    rowwise() %>%
    mutate(max_nse = max(c_across(starts_with('nse') & -contains('neon')), na.rm = TRUE),
           max_kge = max(c_across(starts_with('kge') & -contains('neon')), na.rm = TRUE),
           bestmod_nse = which.max(c_across(starts_with('nse') & -contains('neon'))),
           bestmod_kge = which.max(c_across(starts_with('kge') & -contains('neon')))) %>%
    select(site, contains('neon'), starts_with(c('max', 'bestmod')))
maxd$bestmod_nse <- grep('^nse_', colnames(maxd_), value = TRUE)[-1][maxd$bestmod_nse]
maxd$bestmod_nse <- substr(maxd$bestmod_nse, 5, nchar(maxd$bestmod_nse))
maxd$bestmod_kge <- grep('^kge_', colnames(maxd_), value = TRUE)[-1][maxd$bestmod_kge]
maxd$bestmod_kge <- substr(maxd$bestmod_kge, 5, nchar(maxd$bestmod_kge))

print(maxd, n = 50)

#distribution of best model scores
plot(density(maxd$max_nse))
plot(density(maxd$max_kge))

#median across best model for each site
group_by(maxd) %>%
    summarize(median_nse = median(max_nse, na.rm = TRUE),
              median_kge = median(max_kge, na.rm = TRUE))

#and for neon data
group_by(maxd) %>%
    summarize(median_nse = median(nse_neon_overall, na.rm = TRUE),
              median_kge = median(kge_neon_overall, na.rm = TRUE))

#at how many sites did we achieve better efficiencies?
filter(maxd, max_kge > kge_neon_overall)
filter(maxd, max_nse > nse_neon_overall)

#how many neon sites' published records have efficiencies < 0.7?
filter(maxd, kge_neon_overall < 0.7)
filter(maxd, nse_neon_overall < 0.7)

#mean KGE of the sites we raised above the 0.7 mark
maxd %>%
    filter(kge_neon_overall < 0.7,
           max_kge > 0.7) %>%
    pull(max_kge) %>%
    mean()

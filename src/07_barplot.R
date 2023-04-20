# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-02-28

library(tidyverse)
library(macrosheds)
library(glue)
library(reticulate)
library(data.table)

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

# specify run IDs for all tested generalist models
generalist_runids <- 1468:1520
tecr_runids <- 1718:1747
bigc_runids <-
# specify run IDs for replicates of best model for each specialist type
specialist_runids <- c(2293:2422)
pgdl_runids <- 2248:2292

pal <- c('black', rev(viridis::viridis(5, begin = 0.2, end = 1)))

neon_sites <- read_csv('in/NEON/neon_site_info.csv') %>%
    filter(! SiteType == 'Lake') %>%
    pull(SiteID)

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

tecr_res <- retrieve_test_results(tecr_runids)

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

#fig Sx (same as fig 2, but showing NSE)

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
                           xpd = NA, ncol = 3))
mtext('NEON stream/river site', 1, 3.8, font = 2)
mtext('Nash-Sutcliffe Efficiency', 2, 1, font = 2)
axis(2, seq(0, 1, 0.1), line = -0.9, tcl = -0.3, padj = 1)
dev.off()

#actual fig 2

plotd_kge <- select(plotd, -contains('nse'))
plotd_kge$rowmax <- apply(select(plotd_kge, -site), 1, max, na.rm = TRUE)
plotd_kge <- arrange(plotd_kge, desc(kge_neon))
plotd_m <- as.matrix(select(plotd_kge, -site, -rowmax))
plotd_m[! is.na(plotd_m) & plotd_m < -0.05] <- -0.05
plotd_m <- t(plotd_m)
rownames(plotd_m) <- c('Published', 'Linreg', 'Linreg scaled', 'LSTM generalist', 'LSTM specialist', 'LSTM process-guided')

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
                           xpd = NA, ncol = 3))
mtext('NEON stream/river Site', 1, 3.8, font = 2)
mtext('Kling-Gupta Efficiency', 2, 1, font = 2)
axis(2, seq(0, 1, 0.1), line = -0.9, tcl = -0.3, padj = 1)
dev.off()

## 6. table 4 ####

plotd %>%
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
maxd <- plotd %>%
    rowwise() %>%
    mutate(max_nse = max(c_across(starts_with('nse') & -nse_neon), na.rm = TRUE),
           max_kge = max(c_across(starts_with('kge') & -kge_neon), na.rm = TRUE),
           bestmod_nse = which.max(c_across(starts_with('nse') & -nse_neon)),
           bestmod_kge = which.max(c_across(starts_with('kge') & -kge_neon))) %>%
    select(site, ends_with('neon'), starts_with(c('max', 'bestmod')))
maxd$bestmod_nse <- grep('^nse_', colnames(plotd), value = TRUE)[-1][maxd$bestmod_nse]
maxd$bestmod_nse <- substr(maxd$bestmod_nse, 5, nchar(maxd$bestmod_nse))
maxd$bestmod_kge <- grep('^kge_', colnames(plotd), value = TRUE)[-1][maxd$bestmod_kge]
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
    summarize(median_nse = median(nse_neon, na.rm = TRUE),
              median_kge = median(kge_neon, na.rm = TRUE))

#at how many sites did we achieve better efficiencies?
filter(maxd, max_kge > kge_neon)
# filter(maxd, max_nse > nse_neon)

#how many neon sites' published records have KGE < 0.7?
filter(maxd, kge_neon < 0.7)
# filter(maxd, nse_neon < 0.7)

#mean KGE of the 7 sites we raised above the 0.7 mark *approx
maxd %>%
    filter(kge_neon < 0.7,
           max_kge > 0.69) %>%  #*
    pull(max_kge) %>%
    mean()

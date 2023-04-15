# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-04

library(tidyverse)
library(reticulate)
library(glue)
#four more packages loaded in section 6

#see step 2 in src/lstm_dungeon/README.txt for installing conda environment
use_condaenv('nh2')
xr <- import("xarray")
pd <- import("pandas")
np <- import("numpy")

#pre-bundled in/out data available at: [**]
if(! exists('ts_plot')) source('src/00_helpers.R')
if(! exists('ms_areas')) source('src/01_data_retrieval.R')
if(! dir.exists('in/lstm_data')) source('src/03_organize_camels_macrosheds_nhm.R', local = new.env())

# 1. setup (inputs) ####

#establish run IDs
param_search <- list(generalist = list(1468:1507, #batch 1
                                       1508:1520), #batch 2
                     specialist = list(1548:1627, #batch 1
                                       1748:1937), #batch 2
                     # specialist = list(2293:2307, #batch 1
                     #                   2308:2422), #batch 2
                     pgdl = list(2028:2117)) #but also 1628:1657! just didn't find anything good there i guess
#anyway the solution must be to include SOME or NO search runs?

#just for plotting. any sets that perform well in this stage get full ensembles
mini_ensembles <- list(specialist = list(), #list(2293:2307, 2308:2397),
                       pgdl = list(2248:2292))

ensembles <- list(
    TECR = list(2423:2452),
    BIGC = list(2453:2482),
    MART = list(2453:2482), #same ensemble as BIGC
    LECO = list(2453:2482), #same ensemble as BIGC
    WALK = list(2483:2512),
    MCRA = list(2513:2542), #
    COMO = list(2543:2572),
    HOPB = list(2573:2602),
    # FLNT = list(2603:2632), #wasn't actually using ms NHM for continue, but has legit validation continue tb
    FLNT = list(2633:2662), #PGDL
    # BLDE = list()  #SPECIALIST
)

# 2. setup (rigamarole) ####

#establish i/o directory locations
confdir <- file.path(getwd(), 'in/lstm_configs')
datadir <- file.path(getwd(), 'in/lstm_data')
rundir <- file.path(getwd(), 'out/lstm_runs')

#insert i/o directories into model/config text files
cfgs <- list.files('.', pattern = '.*\\.yml$', recursive = TRUE, full.names = TRUE)
cfgs <- c(cfgs, list.files('.', pattern = 'pretrained_model_loc',
                           recursive = TRUE, full.names = TRUE))
for(cfg_ in cfgs){

    read_file(cfg_) %>%
        str_replace_all('PLACEHOLDER3', datadir) %>%
        str_replace_all('PLACEHOLDER2', rundir) %>%
        str_replace_all('PLACEHOLDER', confdir) %>%
        write_file(cfg_)
}

#build output dir structure
dir.create('out/lstm_out', showWarnings = FALSE)
dir.create('out/lstm_out/predictions', showWarnings = FALSE)
dir.create('out/lstm_out/fit', showWarnings = FALSE)
dir.create('figs/lstm_plots', showWarnings = FALSE)
dir.create('figs/lstm_plots/pred', showWarnings = FALSE)
dir.create('figs/lstm_plots/val', showWarnings = FALSE)

#field Q measured after this date is in the holdout set
holdout_cutuff <- as.Date('2019-12-31')

# 3. run LSTMs parameter searches (each may take several days!) ####

run_lstm('generalist', param_search$generalist[[1]]) #batch 1
run_lstm('generalist', param_search$generalist[[2]]) #batch 2
run_lstm('specialist', param_search$specialist[[1]]) #batch 1
run_lstm('specialist', param_search$specialist[[2]]) #batch 2
run_lstm('pgdl', param_search$pgdl[[1]])

# 4. run ensembles for potentially skilled models (also several days apiece!) ####

run_lstm('generalist', ensembles$TECR[[1]])
run_lstm('generalist', ensembles$BIGC[[1]])
run_lstm('generalist', ensembles$WALK[[1]])
run_lstm('generalist', ensembles$MCRA[[1]])
run_lstm('generalist', ensembles$COMO[[1]])
run_lstm('generalist', ensembles$HOPB[[1]])
run_lstm('pgdl', ensembles$FLNT[[1]])
# run_lstm('specialist', ensembles$BLDE[[1]])

# 5. gather results of ensembles ####

metrics <- matrix(
    NA, nrow = length(ensembles), ncol = 6,
    dimnames = list(
        names(ensembles),
        c('NSE_holdout', 'KGE_holdout', 'pbias_holdout', 'NSE', 'KGE', 'pbias')
    )
)
pred_q <- list()
for(i in seq_along(ensembles)){

    neon_site <- names(ensembles)[i]
    runlist <- unlist(ensembles[[i]])

    neon_q_manual <- read_csv(glue('in/NEON/neon_field_Q/{neon_site}.csv')) %>%
        mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
        rename(discharge_manual = discharge) %>%
        distinct(datetime, .keep_all = TRUE) %>%
        mutate(date = as.Date(datetime))

    rundirs <- list.files('out/lstm_runs/', pattern = '^run') %>%
        str_match(glue('run(?:{rl})_[0-9_]+', rl = paste(runlist, collapse = '|'))) %>%
        {.[, 1]} %>%
        na.omit() %>%
        as.vector()

    ws_area_ha <- na.omit(neon_areas$ws_area_ha[neon_areas$site_code == neon_site])

    ensemble_results <- list()
    for(j in seq_along(rundirs)){

        epoch_dir <- list.files(glue('out/lstm_runs/{r}/test', r = rundirs[j]),
                                full.names = TRUE)
        if(length(epoch_dir) > 1) stop('handle this')

        lstm_out <- reticulate::py_load_object(file.path(epoch_dir, 'test_results.p'))

        pred <- lstm_out[[paste0(neon_site, '_MANUALQ')]]$`1D`$xr$discharge_sim$to_pandas()
        pred <- tibble(date = as.Date(rownames(pred)), Q = pred$`0`) %>%
            rename(discharge_sim = Q) %>%
            #      L/s            mm/d             L/m^3  m^2/ha mm/m   s/d     ha
            mutate(discharge_sim = discharge_sim * 1000 * 1e4 / 1000 / 86400 * ws_area_ha)

        ensemble_results[[j]] <- pred
    }

    ensemble_m <- map(ensemble_results, ~.$discharge_sim) %>%
        reduce(cbind) %>%
        as.matrix() %>%
        unname()
    npreds <- nrow(ensemble_m)
    ensemble_out <- tibble(
        date = ensemble_results[[1]]$date,
        lower_bound_2.5 = rep(NA, npreds),
        mean_run = rep(NA, npreds),
        upper_bound_97.5 = rep(NA, npreds))

    for(j in seq_len(npreds)){

        if(all(is.na(ensemble_m[j, ]))) next
        ensemble_out$mean_run[j] <- mean(ensemble_m[j, ])
        ensemble_out$lower_bound_2.5[j] <- quantile(ensemble_m[j, ], probs = 0.025)
        ensemble_out$upper_bound_97.5[j] <- quantile(ensemble_m[j, ], probs = 0.975)
    }

    pred_q[[neon_site]] <- full_join(neon_q_manual, ensemble_out, by = 'date') %>%
        arrange(date)
    pred_q_filt <- pred_q[[neon_site]] %>%
        filter(! is.na(datetime))

    metrics[i, 1] <- hydroGOF::NSE(sim = pred_q_filt$mean_run[pred_q_filt$date > holdout_cutuff],
                                   obs = pred_q_filt$discharge_manual[pred_q_filt$date > holdout_cutuff])
    metrics[i, 2] <- hydroGOF::KGE(sim = pred_q_filt$mean_run[pred_q_filt$date > holdout_cutuff],
                                   obs = pred_q_filt$discharge_manual[pred_q_filt$date > holdout_cutuff])
    metrics[i, 3] <- hydroGOF::pbias(sim = pred_q_filt$mean_run[pred_q_filt$date > holdout_cutuff],
                                     obs = pred_q_filt$discharge_manual[pred_q_filt$date > holdout_cutuff])
    metrics[i, 4] <- hydroGOF::NSE(sim = pred_q_filt$mean_run,
                                   obs = pred_q_filt$discharge_manual)
    metrics[i, 5] <- hydroGOF::KGE(sim = pred_q_filt$mean_run,
                                   obs = pred_q_filt$discharge_manual)
    metrics[i, 6] <- hydroGOF::pbias(sim = pred_q_filt$mean_run,
                                     obs = pred_q_filt$discharge_manual)

    pred_q_filt %>%
        select(site_code, datetime, Q_neon_field = discharge_manual,
               Q_predicted = mean_run) %>%
        write_csv(glue('out/lstm_out/fit/{neon_site}.csv'))
}

as.data.frame(metrics) %>%
    rownames_to_column('site_code') %>%
    select(site_code, NSE, NSE_holdout, KGE, KGE_holdout, pbias, pbias_holdout) %>%
    write_csv('out/lstm_out/results.csv')

# 6. generate output plots and datasets ####

library(lubridate)
library(xts)
library(dygraphs)
library(htmlwidgets)

for(i in seq_along(pred_q)){

    neon_site <- names(pred_q)[i]
    pred_q_site <- pred_q[[i]]

    neon_q_auto <- read_csv(glue('in/NEON/neon_continuous_Q/{neon_site}.csv')) %>%
        filter(! is.na(discharge)) %>%
        rename(discharge_auto = discharge)

    nadts <- is.na(pred_q_site$datetime)
    pred_q_site$datetime[nadts] <- as_datetime(pred_q_site$date[nadts])
    pred_q_site <- pred_q_site %>%
        full_join(neon_q_auto, by = 'datetime') %>%
        filter(! is.na(discharge_manual) | minute(datetime) %% 5 == 0)

    pred_q_site %>%
        mutate(date = as.Date(datetime)) %>%
        select(date,
               Q_predicted = mean_run,
               Q_pred_int_2.5 = lower_bound_2.5,
               Q_pred_int_97.5 = upper_bound_97.5,
               Q_neon_continuous = discharge_auto) %>%
        filter(! is.na(Q_predicted)) %>%
        arrange(date) %>%
        write_csv(glue('out/lstm_out/predictions/{neon_site}.csv'))

    ## plot predictions vs continuous NEON Q

    plotdata <- pred_q_site %>%
        select(datetime, Q_predicted = mean_run,
               Q_neon_field = discharge_manual,
               Q_neon_continuous = discharge_auto) %>%
        arrange(datetime)

    dygraphs::dygraph(xts(x = select(plotdata, -datetime) %>% tail(5e5),
                                order.by = tail(plotdata$datetime, 5e5))) %>%
        dyRangeSelector() %>%
        saveWidget(glue('figs/lstm_plots/pred/{neon_site}.html'))

    ## plot predictions versus field measurements

    plotdata <- plotdata %>%
        filter(! is.na(Q_neon_field)) %>%
        group_by(datetime) %>%
        summarize(across(everything(), ~mean(., na.rm = TRUE))) %>%
        ungroup()

    axlim <- c(0, max(c(plotdata$Q_predicted, plotdata$Q_neon_field), na.rm = TRUE))

    png(glue('figs/lstm_plots/val/{neon_site}_obs_v_pred.png'),
        6, 6, 'in', type = 'cairo', res = 300)
    plot(plotdata$Q_neon_field, plotdata$Q_predicted, xlab = 'NEON Field Discharge (L/s)',
         ylab = 'Predicted Discharge (L/s)',
         main = glue('Site: {neon_site}; NSE overall: {nse1}; NSE holdout: {nse2}',
                     nse1 = round(metrics[rownames(metrics) == neon_site, 'NSE'], 2),
                     nse2 = round(metrics[rownames(metrics) == neon_site, 'NSE_holdout'], 2)),
         xlim = axlim, ylim = axlim, xaxs = 'i', yaxs = 'i')
    abline(a = 0, b = 1, col = 'blue')
    legend('topleft', legend = '1:1', lty = 1, col = 'blue', bty = 'n')
    dev.off()
}

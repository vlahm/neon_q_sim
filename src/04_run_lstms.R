# Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-04

library(tidyverse)
library(reticulate)
library(glue)
library(hydroGOF)
#four more packages loaded in section 6

#NOTE: IF YOU ARE NOT USING OUR BUNDLED INPUT DATA, THIS FILE REQUIRES EDITS,
#   AS LSTM RUN IDs AND DIRECTORY NAMES WILL CHANGE. YOU'LL ALSO HAVE TO
#   MODIFY THE CONFIGURATION FILES TO POINT TO THE CORRECT PRETRAIN LOCATIONS,
#IF YOU INTEND TO RUN LSTMS ON AN HPC CLUSTER, USE
#   src/lstm_dungeon/run_lstms_hpc.py instead of run_lstms_local.py

#see step 2 in src/lstm_dungeon/README.txt for installing conda environment
use_condaenv('nh2')
xr <- import("xarray")
pd <- import("pandas")
np <- import("numpy")

#pre-bundled in/out data available at: [**]
if(! exists('ts_plot')) source('src/00_helpers.R')
if(! exists('ms_areas')) source('src/01_data_retrieval.R')
if(! dir.exists('in/lstm_data')) source('src/03_organize_camels_macrosheds_nhm.R', local = new.env())

# 1. setup ####

options(readr.show_progress = FALSE,
        readr.show_col_types = FALSE,
        timeout = 3000)

#establish i/o directory locations
confdir <- file.path(getwd(), 'in/lstm_configs')
datadir <- file.path(getwd(), 'in/lstm_data')
rundir <- file.path(getwd(), 'out/lstm_runs')

#insert i/o directories into model/config text files
cfgs <- list.files('in/lstm_configs', pattern = '.*\\.yml$',
                   recursive = TRUE, full.names = TRUE)
cfgs <- c(cfgs, list.files('out/lstm_runs', pattern = 'config.yml',
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


# 2. run parameter searches (each may take several days! use HPC if possible) ####

#establish run IDs
param_search <- list(generalist = 3003:3032,
                     specialist = 3033:3812,
                     pgdl_generalist = 3813:3842,
                     pgdl_specialist = 3843:4112)

run_lstm('generalist', param_search$generalist)
run_lstm('specialist', param_search$specialist)
run_lstm('generalist', param_search$pgdl_generalist)
run_lstm('specialist', param_search$pgdl_specialist)

# 3. evaluate all models. determine where ensembles are warranted ####

all_generalists <- eval_lstms('generalist', param_search$generalist) %>%
    identify_best_models(kge_thresh = -Inf)
skilled_generalists <- filter(all_generalists, kge >= 0.6)

all_specialists <- eval_lstms('specialist', param_search$specialist) %>%
    identify_best_models(kge_thresh = -Inf)
skilled_specialists <- filter(all_specialists, kge >= 0.6)

all_pgdl_generalists <- eval_lstms('generalist', param_search$pgdl_generalist) %>%
    identify_best_models(kge_thresh = -Inf)
skilled_pgdl_generalists <- filter(all_pgdl_generalists, kge >= 0.6)

all_pgdl_specialists <- eval_lstms('specialist', param_search$pgdl_specialist) %>%
    identify_best_models(kge_thresh = -Inf)
skilled_pgdl_specialists <- filter(all_pgdl_specialists, kge >= 0.6)

best_by_strategy <- bind_rows(
        mutate(all_generalists, strategy = 'gen'),
        mutate(all_specialists, strategy = 'spec'),
        mutate(all_pgdl_generalists, strategy = 'pgdl_gen'),
        mutate(all_pgdl_specialists, strategy = 'pgdl_spec')
    ) %>%
    group_by(site_code, strategy) %>%
    filter(kge == max(kge)) %>%
    ungroup() %>%
    arrange(strategy, desc(kge))

write_csv(best_by_strategy, 'out/lstm_out/param_search_skill.csv')

skilled <- bind_rows(
        mutate(skilled_generalists, strategy = 'gen'),
        mutate(skilled_specialists, strategy = 'spec'),
        mutate(skilled_pgdl_generalists, strategy = 'pgdl_gen'),
        mutate(skilled_pgdl_specialists, strategy = 'pgdl_spec')
    ) %>%
    group_by(site_code) %>%
    filter(kge == max(kge)) %>%
    ungroup() %>%
    arrange(desc(kge)) %>%
    print(n = 100)

# skilled = read_csv('~/Desktop/dcc_runs/skilled_searches.csv')
# distinct(skilled, run, .keep_all = T)

# 4. build ensemble configs ####

unique_runs <- unique(skilled$run)
for(run in unique_runs){

    #some generalists perform well on multiple sites
    sites <- filter(skilled, run == !!run)$site_code

    build_ensemble_config(sites = sites, runid = run, param_search = param_search)
}

# 5. run ensembles for potentially skilled models (also several days apiece!) ####

ensembles_gen <- list(
    TECR = 4113:4142, #same ensemble as FLNT
    # FLNT = 4113:4142, #same ensemble as TECR
    HOPB = 4173:4202,
    MART = 4293:4322,
    MCRA = 4323:4352, #same ensemble as BLDE and MCRA
    # BLDE = 4323:4352, #
    # MAYF = 4323:4352, #
    LECO = 4473:4502,
    WALK = 4623:4652
)

ensembles_spec <- list(
    BIGC = 4143:4172,
    LEWI = 4203:4232,
    BLUE = 4383:4412,
    BLWA = 4413:4442,
    CUPE = 4443:4472,
    CARI = 4563:4592
)

ensembles_gen_pgdl <- list(
    WLOU = 4233:4262,
    KING = 4263:4292,
    MCDI = 4353:4382,
    COMO = 4503:4532,
    REDB = 4533:4562,
    POSE = 4593:4622
)

## probably want to run one at a time
# for(ensemb_site in ensembles_gen){
#     run_lstm('generalist', ensemb_site, ensemble = TRUE)
# }
# for(ensemb_site in ensembles_spec){
#     run_lstm('specialist', ensemb_site, ensemble = TRUE)
# }
# for(ensemb_site in ensembles_gen_pgdl){
#     run_lstm('generalist', ensemb_site, ensemble = TRUE)
# }

run_lstm('generalist', ensembles_gen$TECR, ensemble = TRUE)
# run_lstm('generalist', ensembles_gen$FLNT, ensemble = TRUE)
run_lstm('generalist', ensembles_gen$HOPB, ensemble = TRUE)
run_lstm('generalist', ensembles_gen$MART, ensemble = TRUE)
run_lstm('generalist', ensembles_gen$MCRA, ensemble = TRUE)
# run_lstm('generalist', ensembles_gen$BLDE, ensemble = TRUE)
# run_lstm('generalist', ensembles_gen$MAYF, ensemble = TRUE)
run_lstm('generalist', ensembles_gen$LECO, ensemble = TRUE)
run_lstm('generalist', ensembles_gen$WALK, ensemble = TRUE)

run_lstm('specialist', ensembles_gen$BIGC, ensemble = TRUE)
run_lstm('specialist', ensembles_gen$LEWI, ensemble = TRUE)
run_lstm('specialist', ensembles_gen$BLUE, ensemble = TRUE)
run_lstm('specialist', ensembles_gen$BLWA, ensemble = TRUE)
run_lstm('specialist', ensembles_gen$CUPE, ensemble = TRUE)
run_lstm('specialist', ensembles_gen$CARI, ensemble = TRUE)

run_lstm('generalist', ensembles_gen$WLOU, ensemble = TRUE)
run_lstm('generalist', ensembles_gen$KING, ensemble = TRUE)
run_lstm('generalist', ensembles_gen$MCDI, ensemble = TRUE)
run_lstm('generalist', ensembles_gen$COMO, ensemble = TRUE)
run_lstm('generalist', ensembles_gen$REDB, ensemble = TRUE)
run_lstm('generalist', ensembles_gen$POSE, ensemble = TRUE)

## you might need to sort out runs that failed due to race conditions
## if you ran models in parallel

# status_gen <- tibble()
# for(ensemb_site in ensembles_gen){
#     status_gen <- bind_rows(status_gen, identify_incomplete_lstms('generalist', ensemb_site))
# }
# status_spec <- tibble()
# for(ensemb_site in ensembles_spec){
#     status_spec <- bind_rows(status_spec, identify_incomplete_lstms('specialist', ensemb_site))
# }
# status_gen_pgdl <- tibble()
# for(ensemb_site in ensembles_gen_pgdl){
#     status_gen_pgdl <- bind_rows(status_gen_pgdl, identify_incomplete_lstms('generalist', ensemb_site))
# }


# 6. evaluate ensembles; gather results ####

## if you ran ensembles on a cluster, you'll need to adjust paths for your local machine
## using this commented chunk

# system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_data|PLACEHOLDER3|g' -i"))
# system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_runs|PLACEHOLDER2|g' -i"))
# system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/runs|PLACEHOLDER2|g' -i"))
# system(paste0("find out/lstm_runs -name 'config.yml' | xargs sed -e 's|/hpc/home/mjv22/q_sim/lstm_configs|PLACEHOLDER|g' -i"))
# cfgs <- list.files('out/lstm_runs', pattern = 'config.yml', recursive = TRUE, full.names = TRUE)
# for(cfg_ in cfgs){
#     read_file(cfg_) %>%
#         str_replace_all('PLACEHOLDER3', datadir) %>%
#         str_replace_all('PLACEHOLDER2', rundir) %>%
#         str_replace_all('PLACEHOLDER', confdir) %>%
#         write_file(cfg_)
# }

for(ensemb_site in ensembles_gen){
    eval_lstms('generalist', ensemb_site)
}
for(ensemb_site in ensembles_spec){
    eval_lstms('specialist', ensemb_site)
}
for(ensemb_site in ensembles_gen_pgdl){
    eval_lstms('generalist', ensemb_site)
}

metrics <- build_metrics_skeleton(ensembles_gen, ensembles_spec, ensembles_gen_pgdl)
pred_q <- list()
for(i in seq_len(nrow(metrics))){
# cntr <- 1
# for(i in seq_along(ensembles)){

    neon_site <- metrics$site[i]
    phase <- metrics$phase[i]
    runlist <- metrics$runlist[[i]]

    neon_q_manual <- read_csv(glue('in/NEON/neon_field_Q/{neon_site}.csv')) %>%
        mutate(discharge = ifelse(discharge < 0, 0, discharge)) %>%
        rename(discharge_manual = discharge) %>%
        distinct(datetime, .keep_all = TRUE) %>%
        mutate(date = as.Date(datetime))

    rundirs <- list.files('out/lstm_runs/', pattern = paste0('^', phase)) %>%
        str_match(glue('{phase}(?:{rl})_[0-9_]+', rl = paste(runlist, collapse = '|'))) %>%
        {.[, 1]} %>%
        na.omit() %>%
        as.vector()

    ensemble_results <- list()
    for(j in seq_along(rundirs)){

        epoch_dir <- list.files(glue('out/lstm_runs/{r}/test', r = rundirs[j]),
                                full.names = TRUE)
        if(length(epoch_dir) > 1) stop('handle this')

        lstm_out <- reticulate::py_load_object(file.path(epoch_dir, 'test_results.p'))

        ws_area_ha <- na.omit(neon_areas$ws_area_ha[neon_areas$site_code == neon_site])

        pred <- lstm_out[[paste0(neon_site, '_MANUALQ')]]$`1D`$xr$discharge_sim$to_pandas()
        pred <- tibble(date = as.Date(rownames(pred)), Q = pred$`0`) %>%
            rename(discharge_sim = Q) %>%
            #      L/s            mm/d             L/m^3  m^2/ha mm/m   s/d     ha
            mutate(discharge_sim = discharge_sim * 1000 * 1e4 / 1000 / 86400 * ws_area_ha,
                   site_code = neon_site)

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

    # plot(pred_q_filt$datetime, pred_q_filt$mean_run, type = 'l')
    # points(pred_q_filt$datetime, pred_q_filt$discharge_manual)

    metrics$KGE[i] <- KGE(sim = pred_q_filt$mean_run, obs = pred_q_filt$discharge_manual)
    metrics$NSE[i] <- NSE(sim = pred_q_filt$mean_run, obs = pred_q_filt$discharge_manual)
    metrics$pbias[i] <- pbias(sim = pred_q_filt$mean_run, obs = pred_q_filt$discharge_manual)

    pred_q_filt %>%
        select(site_code, datetime, Q_neon_field = discharge_manual,
               Q_predicted = mean_run) %>%
        write_csv(glue('out/lstm_out/fit/{neon_site}.csv'))
}

metrics %>%
    select(site_code = site, strategy, KGE, NSE, pbias) %>%
    arrange(desc(KGE)) %>%
    write_csv('out/lstm_out/results.csv')

# 7. generate output plots and datasets ####

library(lubridate)
library(xts)
library(dygraphs)
library(htmlwidgets)

#don't save predictions for sites with weak ensembles.
legit_ensemb <- filter(metrics, KGE >= 0.5, NSE >= 0.4) %>% pull(site)
pred_q <- Filter(function(x) Find(function(y) ! is.na(y), x$site_code) %in% legit_ensemb,
                 pred_q)

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
               Q_neon_release2023_qcpass = discharge_auto) %>%
        filter(! is.na(Q_predicted)) %>%
        arrange(date) %>%
        distinct(date, .keep_all = TRUE) %>%
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
         main = glue('Site: {neon_site}; KGE: {kge}; NSE: {nse}; pbias: {pb}',
                     nse = round(metrics$NSE[metrics$site == neon_site], 2),
                     kge = round(metrics$KGE[metrics$site == neon_site], 2),
                     pb = round(metrics$pbias[metrics$site == neon_site], 2)),
         xlim = axlim, ylim = axlim, xaxs = 'i', yaxs = 'i')
    abline(a = 0, b = 1, col = 'blue')
    legend('topleft', legend = '1:1', lty = 1, col = 'blue', bty = 'n')
    dev.off()
}

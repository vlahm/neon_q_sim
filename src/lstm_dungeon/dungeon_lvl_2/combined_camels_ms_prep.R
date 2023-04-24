library(feather)
# library(plyr)
library(tidyverse)
library(lubridate)
library(foreach)
library(doParallel)
library(glue)
library(yaml)
library(reticulate)
library(fuzzyjoin)
library(data.table)
library(errors)
library(zoo)
library(ncdf4)

options(readr.show_progress = FALSE)
options(readr.show_col_types = FALSE)

# rundeets_ = read_lines(file.path(grep('run1468', list.files('~/git/macrosheds/papers/q_sim/out/lstm_runs', full.names = T), value = T), 'config.yml'))
# grep('batch|hidden_|output_d|epochs', rundeets_, value = T)
# xid = grep('finetune_modules', rundeets_)
# rundeets_[xid:(xid+2)]
# xid = grep('learning_rate', rundeets_)
# rundeets_[xid:(xid+3)]

#FOR GENERALIST ENSEMBLES (AND MAYBE ONE PGDL?)
# xsite='TECR'
# hidden_size_ = 30L; output_dropout_ = 0.4; batch_size_ = 32L; epochs_ = 20L
# finetune_modules_ = list(list('head', 'lstm')))
# learning_rate_ = list(list(a = '1e-4', b = '1e-5', c = '1e-6'))
# xsite='BIGC'
# hidden_size_ = 30L; output_dropout_ = 0.2; batch_size_ = 128L; epochs_ = 30L
# finetune_modules_ = list(list('head', 'lstm')))
# learning_rate_ = list(list(a = '1e-4', b = '1e-5', c = '1e-6'))
# xsite='WALK'
# batch_size_ = 256L; epochs_ = 30L; hidden_size_ = 30L; output_dropout_ = 0.5;
# finetune_modules_ = list(list('head', 'lstm'))
# learning_rate_ = list(list(a = '1e-4', b = '1e-5', c = '1e-6'))
# xsite='MCRA'
# batch_size_ = 256L; epochs_ = 30L; hidden_size_ = 30L; output_dropout_ = 0.2;
# finetune_modules_ = list(list('head', 'lstm'))
# learning_rate_ = list(list(a = '1e-3', b = '1e-4', c = '1e-5'))
# xsite='COMO'
# batch_size_ = 64L; epochs_ = 40L; hidden_size_ = 30L; output_dropout_ = 0.5;
# finetune_modules_ = list(list('head', 'lstm'))
# learning_rate_ = list(list(a = '1e-3', b = '1e-4', c = '1e-5'))
# xsite='HOPB'
# batch_size_ = 32L; epochs_ = 20L; hidden_size_ = 30L; output_dropout_ = 0.5;
# finetune_modules_ = list(list('head', 'lstm'))
# learning_rate_ = list(list(a = '1e-3', b = '1e-4', c = '1e-5'))
# xsite='FLNT'
# batch_size_ = 512L; epochs_ = 30L; hidden_size_ = 30L; output_dropout_ = 0.5;
# finetune_modules_ = list(list('head', 'lstm'))
# learning_rate_ = list(list(a = '1e-2', b = '1e-3', c = '1e-4'))

#SPECIALIST MINI-ENSEMBLE? these params get overwritten in postprocessing
# batch_size_ = 128L; epochs_ = 30L; hidden_size_ = 30L; output_dropout_ = 0.2;
# finetune_modules_ = list(list('head', 'lstm'))
# learning_rate_ = list(list(a = '1e-4', b = '1e-5', c = '1e-6'))

#FOR SPECIALIST ENSEMBES, PIGGYBACKING ON GENERALIST ENSEMBLES
# xsite='HOPB'
# run_range = 2573:2602 #generalist ensemble range to piggyback from
# batch_size_ = 512L; epochs_ = 20L; hidden_size_ = 30L; output_dropout_ = 0.3;
# finetune_modules_ = list(list('head', 'lstm'))
# learning_rate_ = list(list(a = '1e-4', b = '1e-5', c = '1e-6'))
# xsite='WALK'
# run_range = 2483:2512
# batch_size_ = 128L; epochs_ = 30L; hidden_size_ = 30L; output_dropout_ = 0.2;
# finetune_modules_ = list(list('head', 'lstm'))
# learning_rate_ = list(list(a = '1e-5', b = '1e-6', c = '1e-7'))
#xsite='MCRA'
#run_range = 2453:2482 #same ensemble as BIGC
#batch_size_ = 1024L; epochs_ = 30L; hidden_size_ = 30L; output_dropout_ = 0.4;
#finetune_modules_ = list(list('head', 'lstm'))
#learning_rate_ = list(list(a = '1e-4', b = '1e-5', c = '1e-6'))
# xsite='COMO'
# run_range = 2543:2572
# batch_size_ = 128L; epochs_ = 40L; hidden_size_ = 30L; output_dropout_ = 0.3;
# finetune_modules_ = list(list('head', 'lstm'))
# learning_rate_ = list(list(a = '1e-4', b = '1e-5', c = '1e-6'))
# xsite='BIGC'
# run_range = 2453:2482
# batch_size_ = 256L; epochs_ = 30L; hidden_size_ = 30L; output_dropout_ = 0.3;
# finetune_modules_ = list(list('head', 'lstm'))
# learning_rate_ = list(list(a = '1e-3', b = '1e-4', c = '1e-5'))
 # xsite='MART'
 # run_range = 2513:2542
 # batch_size_ = 1024L; epochs_ = 30L; hidden_size_ = 30L; output_dropout_ = 0.4;
 # finetune_modules_ = list(list('head', 'lstm'))
 # learning_rate_ = list(list(a = '1e-4', b = '1e-5', c = '1e-6'))

#SPECIALIST: borrow from generalist search to construct 30*26 specialist search
# run_range = 3003:3032 #i don't think this ended up being used?

#NOTE that for a custom SPECIALIST, like NHC+UNHC, you may need to follow up with these:
#f=test.txt; tee run*/$f < ../runs_1044-1044/run1044/$f (for copying custom finetune train/val/test ranges and test.txt)
#find . -name 'run*' | xargs -n1 -I {} touch {}/finetune.txt (for creating teeable destinations if files don't exist)
#tee run*/finetune.txt < ../runs_1045-1054/finetune.txt (for copying finetune.txt)

#for a custom GENERALIST, prepare it just like custom specialist, then use these:
#find . -name 'finetune*' -delete
#f=test.txt; tee run*/$f < ../runs_1044-1044/run1044/$f
#tee run*/test_ranges.csv < ../runs_1105-1134/finetune_test_ranges.csv

#
#TODO: search clip_gradient_norm and seq_length

set_seeds = FALSE
if(set_seeds) set.seed(10)

### 0a setup (REQUIRES INPUT, possible modification) ####

ms_vsn = 1

ms_domains_avail <- list.files('~/git/macrosheds/qa_experimentation/data/ms_in_camels_format')
ms_domains_avail <- ms_domains_avail[! grepl('camels|NC', ms_domains_avail)]

# read_csv('/home/mike/git/macrosheds/portal/data/general/site_data.csv') %>%
#     filter(domain == 'baltimore',
#            site_type == 'stream_gauge') %>%
#     pull(site_code)

neon_neighbors = read_csv('/home/mike/git/macrosheds/qa_experimentation/data/informative_stuff/nearest_neon_neighbors.csv')
neon_neighbors = distinct(neon_neighbors, nearest_neighbor, .keep_all = TRUE)
nghbs = neon_neighbors$nearest_neighbor

sites_with_wonky_q <- read_csv('/home/mike/git/macrosheds/qa_experimentation/data/informative_stuff/ms_q_status.csv') %>%
    filter(! q_status %in% c('okay', 'short') |
           site_code == 'Rustlers') #too short to run
mbwn <- read_csv('/home/mike/git/macrosheds/qa_experimentation/data/informative_stuff/ms_basins_with_netcdfs.csv') %>%
    mutate(domain = ifelse(is.na(domain), 'neon', domain))
    # filter(! domain %in% c('santa_barbara', 'neon', 'bear', 'NC'))
ms_basins_with_netcdfs <- as.list(mbwn$site_code)
names(ms_basins_with_netcdfs) <- mbwn$domain


#NOT REMOVING SITES WITH WONKY Q                                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# ms_basins_with_netcdfs = ms_basins_with_netcdfs[! ms_basins_with_netcdfs %in% sites_with_wonky_q$site_code]
    # c('RedondoMeadow', 'UpperRedondo', 'Avery', 'Bigelow', 'FLUME_MCZOB', 'Roche_Moutonnee_Creek_Main', 'Gothic', 'Trevor_Creek_Main')]

#set paths (don't use "~")
camels_dir <- '/home/mike/git/macrosheds/qa_experimentation/data/CAMELS/basin_dataset_public_v1p2'
# combined_camels_ms_dir <- '/home/mike/git/macrosheds/qa_experimentation/data/CAMELS_macrosheds_combined'
combined_camels_ms_dir <- '/home/mike/git/macrosheds/papers/q_sim/in/lstm_data'
ms_dir <- glue('/home/mike/git/macrosheds/data_acquisition/macrosheds_dataset_v{ms_vsn}')
# write_dir <- '/home/mike/git/macrosheds/qa_experimentation/imputation/data/nh_methods/run_data'
config_dir <- '/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs'
kratzert_basin_list_loc <- '/home/mike/git/macrosheds/qa_experimentation/imputation/data/nh_methods/kratzert_basin_list.txt'
custom_basin_list_loc <- '/home/mike/git/macrosheds/qa_experimentation/imputation/data/nh_methods/custom_basin_list.txt'
python_executable_path <- '/home/mike/anaconda3/envs/nh/bin/python'
ms_sitelist_loc <- '/home/mike/git/macrosheds/qa_experimentation/data/informative_stuff/site_data2.csv'

kratzert_basin_list <- read_lines(kratzert_basin_list_loc)
# kratzert_basin_list <- union(kratzert_basin_list, addtl_config_list$test_basins)

ms_site_data <- read_csv(ms_sitelist_loc)

#runs are automatically named "runX" where X is an incrementing integer

neon_sites <- c('FLNT', 'WALK', 'WLOU', 'LEWI', 'HOPB', 'MART',
               'BIGC', 'REDB', 'COMO', 'BLDE', 'BLWA', 'PRIN', 'POSE',
               'MCDI', 'MCRA', 'CARI', 'ARIK', 'GUIL', 'KING',
               'MAYF', 'LECO', 'BLUE', 'SYCA', 'TECR', 'CUPE',
               'OKSR')#, 'TOMB')

nearest_neighbors_gapped = read_csv('~/git/macrosheds/qa_experimentation/data/informative_stuff/nearest_neon_neighbors.csv') %>%
    distinct(nearest_neighbor) %>% pull() %>%
    paste0(., '_GAPPED')

nearest_neighbors_gapped_ms = grep('^[0-9]+_GAPPED$', nearest_neighbors_gapped, value = TRUE, invert = TRUE) %>%
    str_match(., '(.*?)_GAPPED$') %>%
    {.[, 2]}
nngm_sitedom = ms_site_data %>%
    filter(site_code %in% nearest_neighbors_gapped_ms) %>%
    select(domain, site_code) %>%
    distinct()
gap_target_sites = as.list(paste0(nngm_sitedom$site_code, '_GAPPED'))
names(gap_target_sites) = nngm_sitedom$domain

# strategy = 'generalist'
strategy = 'specialist'
# strategy = ''
ugly_fix_12000 = TRUE #for specialist runs, to ensure train/val/test ranges are correct in section 4aQ

#use this to run the SPECIALIST for each of many sites (nhm+camels, nhm+ms-target, nhm finetune, test) OR
#the GENERALIST for each of many sites (nhm+camels, nhm+ms-target, test).
#this is not needed for BASIC runs like: nhm+camels, nhm+ms-xyz, test)
#SEE NOTE NEAR THE TOP FOR CUSTOM FILE MANIPULATIONS
# run_prediction_routine_across_sites <- TRUE #SPECIALIST, or configure_GENERALIST
run_prediction_routine_across_sites = if(strategy %in% c('specialist', 'generalist')) T else F
new_generalist <- if(strategy == 'generalist') T else F
# sites_for_routine <- 'NHC'
# sites_for_routine <- unlist(ms_basins_with_netcdfs, use.names = FALSE)
# sites_for_routine = nearest_neighbors_gapped
# sites_for_routine <- paste0(neon_sites, '_extrapolate') #GENERALIST
# sites_for_routine <- paste0(c('ARIK', 'BIGC', 'BLUE', 'BLWA', 'FLNT', 'HOPB', 'KING', 'PRIN'),#, 'SYCA', 'TOMB'),
#                             '_extrapolate') #SPECIALIST
manualqsites = list.files(file.path(combined_camels_ms_dir, 'time_series'), pattern = 'MANUALQ.nc$');
# sites_for_routine <- sub('.nc', '', grep('(?:SYCA|REDB)', manualqsites, value = TRUE))
# sites_for_routine <- 'REDB_MANUALQ'
# sites_for_routine <- grep('TOMB', paste0(neon_sites, '_MANUALQ'), invert = T, value = T)[1] #useful for generalist (still runs a mod per site)
# sites_for_routine <- grep('CUPE|PRIN|OKSR|CARI|ARIK|BIGC|WLOU|BLDE', paste0(neon_sites, '_MANUALQ'), value = T)
# sites_for_routine <- grep('FLNT|TOMB|BLWA|BLUE|HOPB|BIGC|SYCA|PRIN|ARIK|KING', paste0(neon_sites, '_MANUALQ'), value = T) #all NHM runs
# sites_for_routine <- grep('CUPE|TOMB|FLNT|SYCA|HOPB|TECR|BLUE|REDB|KING|WALK|BLWA|MCDI|MCRA|MAYF|LECO|COMO|MART|LEWI|POSE|GUIL', paste0(neon_sites, '_MANUALQ'), value = T) #remainder of specialist runs
sites_for_routine <- grep('WLOU|BIGC|BLDE|PRIN|CARI|ARIK|OKSR|CUPE|TOMB|FLNT|SYCA|HOPB|TECR|BLUE|REDB|KING|WALK|BLWA|MCDI|MCRA|MAYF|LECO|COMO|MART|LEWI|POSE|GUIL', paste0(neon_sites, '_MANUALQ'), value = T) #all sites
# sites_for_routine <- grep('WLOU|BIGC|BLDE|PRIN|CARI|ARIK|OKSR|TOMB|FLNT|SYCA|HOPB|BLUE|KING|WALK|BLWA|MCDI|MCRA|MAYF|LECO|COMO|MART|LEWI|POSE|GUIL', paste0(neon_sites, '_MANUALQ'), value = T)
# sites_for_routine <- grep(xsite, paste0(neon_sites, '_MANUALQ'), value = T)
# sites_for_routine <- grep('REDB|TECR|CUPE', paste0(neon_sites, '_MANUALQ'), value = T)
# sites_for_routine <- rep(grep('TECR', paste0(neon_sites, '_MANQ_EXTRAP'), value = T), 1)
# 'SYCA', # 'OKSR', # 'TOMB',
if(! run_prediction_routine_across_sites) sites_for_routine <- 'placeholder'

# arctic_sites_to_drop <- c('Kuparuk_River_-0.1', 'Kuparuk_River_-0.177', 'Kuparuk_River_-0.3','Kuparuk_River_-0.47',
#                           'Kuparuk_River_-0.7', 'Kuparuk_River_0.3', 'Kuparuk_River_0.74',
#                           'Kuparuk_River_1.39', 'Kuparuk_River_1.5', 'Kuparuk_River_1.8', 'Kuparuk_River_0', 'Kuparuk_River_2.5',
#                           'Kuparuk_River_2', 'Kuparuk_River_3', 'Kuparuk_River_4.1', 'Kuparuk_River_4', 'Oksrukuyik_Creek_-0.1',
#                           'Oksrukuyik_Creek_-0.3', 'Oksrukuyik_Creek_-0.7', 'Oksrukuyik_Creek_-1.2', 'Oksrukuyik_Creek_0.23',
#                           'Oksrukuyik_Creek_0.48', 'Oksrukuyik_Creek_1.06', 'Oksrukuyik_Creek_1.37', 'Oksrukuyik_Creek_1.7',
#                           'Oksrukuyik_Creek_2.5', 'Oksrukuyik_Creek_2.7',
#                           'Kuparuk_River_-0', 'Oksrukuyik_Creek_-0', 'Oksrukuyik_Creek_-1',
#                           'Oksrukuyik_Creek_0', 'Oksrukuyik_Creek_2')
arctic_sites_to_drop = c()

# ix = 1
siterange <- if(strategy == 'generalist') 1 else seq_along(sites_for_routine)
for(ix in siterange){

    routine_site <- sites_for_routine[ix]

    #the site to use for final finetune (if the target site is different)
    specialist_finetune_site_elect = sub('_MANUALQ', '', routine_site) #----------------------------!!!!!!!!!!! SPECIALIST and NHM...
    specialist_finetune_site_elect = paste0('NHM_', specialist_finetune_site_elect)#--------------!!!!!!!!!!! NHM
    # specialist_finetune_site_elect = NA #---------------------------------------------------------!!!!!!!!!!! GENERALIST

#"additional" configuration parameters are set here. the standard neuralhydrology
#   parameters are set in section 3.
addtl_config_list <- list(

    #add a descriptive note for this run
    run_note = 'DCC PGDL generalist search',
    # run_note = paste("specialist ensembles forreals:", xsite),
    # run_note = "specialist subset (all but tecr, cupe, redb) with fixed params 5x",
    # run_note = "all pgdl with fixed params 5x",
    # run_note = "wave2 2 WALK specialist with doe pretrain (manual config edit)",
    # run_note = "NHC+UNHC 10 (custom generalist)",

    #if doing NHM run, this will be ignored.
    # pretrained_model = '/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/runs/runs_1628-1657_pretrain_2309_212115/', #nhm with pet (NOT TRUE! no nhm sites in pretrain)
    # pretrained_model = '/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/runs/runs_1468-1507_pretrain_2209_135501/', #specialist or generalist
    # pretrained_model = '/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/lstm_runs/generalist_specialist_pretrain_2104_203937/', #new generalist and spec (doesn't actually exist at this location, but on dcc)
    pretrained_model = '/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/lstm_runs/generalist_specialist_pretrain_NHM_2304_170417/', #new NHM generalist and spec (doesn't actually exist at this location, but on dcc)
    # pretrained_model = NA,

    #choose CAMELS subset to work from. must be one of:
    #   'kratzert' (531 CAMELS sites)
    #   'full' (674 CAMELS sites)
    #   'custom' (CAMELS set designed to omit basins with questionable areas, and
    #       include a nice distribution of basin areas even if n_train_val_basins
    #       is small)
    basin_set = 'kratzert',
    #choose which CAMELS basins to test on (these will be combined with
    #   ms_test_basins if include_ms_data is TRUE).
    # camels_test_basins = c('06879650', '14158790', '09510200', '01187300'),
    # camels_test_basins = c("02198100", "11176400", "13235000", "06344600", "05507600", "04127918", "06447500", "05592575", "14305500",
    #                        "06906800", "01440000", "02077200", "14182500", "08178880", "06889200", "11480390", "09505350", "02108000",
    #                        "09066200", "09352900", "01181000", "11143000", "01510000", "08164600", "07145700", "01606500", "02143000",
    #                        "01466500", "04015330", "01350000", "08158810", "03504000", "04224775", "08165300", "01451800", "12025700",
    #                        "06339100", "08200000", "05595730", "02450250", "01139800", "02202600", "01516500", "09065500", "01195100",
    #                        "03028000", "03340800", "04197170", "08104900", "02074500", "11468500", "06879650", "05399500"),
    # camels_test_basins = grep('^[0-9]+_GAPPED$', nearest_neighbors_gapped, value = TRUE),
    camels_test_basins = c(),
    camels_only_train = TRUE, #SPECIALIST: T, GENERALIST: T, BASIC: T?

    #choose whether to include MacroSheds sites in this run (TRUE/FALSE)
    include_ms_data = TRUE,
    #domain must be omitted here if you want any of its sites to be omitted from train/val
    ms_train_domains = c('hbef',
                         'hjandrews','konza', 'santa_barbara', #'mcmurdo', #!santa
                         'baltimore', 'arctic',
                         'niwot','bonanza', 'plum', #!plum
                         'neon',
                         'boulder','catalina_jemez','calhoun','luquillo',
                         'fernow','east_river','usgs','krew','suef',
                         'santee','krycklan','walker_branch','shale_hills'),
    #in addition to the above (if you want to specify individual sites)
    # ms_train_basins = c('w1', 'w3', 'w5', 'w7', 'w9'),
    # ms_train_basins = c('GFCP', 'GFGL', 'POBR', 'BARN', 'MAWI'),
    ms_train_basins = c(), #GENERALIST and SPECIALIST and BASIC?
    # ms_train_basins = unlist(ms_basins_with_netcdfs, use.names = F), #rly not SPECIALIST?
    #choose basins to test on (ignored if include_ms_data is FALSE)
    # ms_test_basins = ms_basins_with_netcdfs,
    # ms_test_basins = NULL,
    # ms_test_basins = list(shale_hills = routine_site),
    # ms_test_basins = list(shale_hills = 'SH_weir'),
    # ms_test_basins = list(neon = 'WALK'), #****
    # ms_test_basins = list(baltimore = 'GFGB', #****
    #                       hjandrews = 'GSLOOK',
    #                       # neon = 'TOMB', #no Q?
    #                       # neon = 'BLWA', #no Q?
    #                       neon = 'BLUE',
    #                       neon = 'HOPB',
    #                       neon = 'BIGC',
    #                       neon = 'SYCA',
    #                       neon = 'PRIN',
    #                       neon = 'ARIK',
    #                       neon = 'KING',
    #                       konza = 'N01B',
    #                       # plum = 'upper_ipswich', #no Q
    #                       boulder = 'BC_SW_4',
    #                       calhoun = 'weir_4',
    #                       east_river = 'Gothic',
    #                       santee = 'WS78',
    #                       fernow = 'WS-2',
    #                       # east_river = 'EBC', #no Q
    #                       shale_hills = 'SCAL'),
    # ms_test_basins = list(baltimore = 'GFGB', #****
    #                       bonanza = 'C4',
    #                       arctic = 'Kuparuk_River_-0.1',
    #                       suef = 'GSCC02',
    #                       east_river = 'Gothic',
    #                       boulder = 'BC_SW_4',
    #                       catalina_jemez = 'UpperJaramillo',
    #                       konza = 'N04D',
    #                       walker_branch = 'east_fork',
    #                       santee = 'WS78',
    #                       baltimore = 'GFCP',
    #                       fernow = 'WS-2',
    #                       neon = 'WALK',
    #                       calhoun = 'weir_4'),

    # ms_test_basins = list(NC = 'NHC'), #NHC modified SPECIALIST, not gen
    # ms_test_basins = list(NC = 'NHC',
    #                       NC = 'UNHC'),
    ms_test_basins = list(neon = routine_site), #GENERALIST AND SPECIALIST
    # ms_test_basins = list(neon = 'FLNT_MANUALQ',
    #                       neon = 'WALK_MANUALQ',
    #                       neon = 'WLOU_MANUALQ',
    #                       neon = 'LEWI_MANUALQ',
    #                       neon = 'TECR_MANUALQ',
    #                       neon = 'HOPB_MANUALQ',
    #                       neon = 'MART_MANUALQ',
    #                       neon = 'BIGC_MANUALQ',
    #                       neon = 'REDB_MANUALQ',
    #                       neon = 'SYCA_MANUALQ',
    #                       neon = 'COMO_MANUALQ',
    #                       neon = 'BLDE_MANUALQ',
    #                       neon = 'BLWA_MANUALQ',
    #                       neon = 'PRIN_MANUALQ',
    #                       neon = 'POSE_MANUALQ',
    #                       neon = 'MCDI_MANUALQ',
    #                       neon = 'MCRA_MANUALQ',
    #                       neon = 'CARI_MANUALQ',
    #                       neon = 'OKSR_MANUALQ',
    #                       neon = 'TOMB_MANUALQ',
    #                       neon = 'ARIK_MANUALQ',
    #                       neon = 'GUIL_MANUALQ',
    #                       neon = 'CUPE_MANUALQ',
    #                       neon = 'KING_MANUALQ',
    #                       neon = 'MAYF_MANUALQ',
    #                       neon = 'LECO_MANUALQ',
    #                       neon = 'BLUE_MANUALQ'),
    # ms_test_basins = list(neon = 'FLNT',
    #                       neon = 'WALK',
    #                       neon = 'WLOU',
    #                       neon = 'LEWI',
    #                       neon = 'TECR',
    #                       neon = 'HOPB',
    #                       neon = 'MART',
    #                       neon = 'BIGC',
    #                       neon = 'REDB',
    #                       # neon = 'SYCA',
    #                       neon = 'COMO',
    #                       neon = 'BLDE',
    #                       neon = 'BLWA',
    #                       neon = 'PRIN',
    #                       neon = 'POSE',
    #                       neon = 'MCDI',
    #                       neon = 'MCRA',
    #                       neon = 'CARI',
    #                       # neon = 'OKSR',
    #                       # neon = 'TOMB',
    #                       neon = 'ARIK',
    #                       neon = 'GUIL',
    #                       neon = 'CUPE',
    #                       neon = 'KING',
    #                       neon = 'MAYF',
    #                       neon = 'LECO',
    #                       neon = 'BLUE'),
    # ms_test_basins = gap_target_sites,

    #should test basins be included in the training/validation set?
    learn_from_test_basins = FALSE, #might be useful for SPECIALIST? but so far F
    #should transfer learning be used to specialize the model toward the test basins?
    #   see section 2b for fine-tuning config parameters
    #   SEE ALSO include_test_basins_in_finetune (which must be TRUE?)
    finetune_to_test_basins = ifelse(strategy == 'specialist', T, F), #writes two extra config files #SPECIALIST
    #pretrain on one set, finetune on another, evaluate on a third. this way, you can use all
    #available data for each pretrain and finetune site in train/val, and save all the testing
    #for the test set of sites. probably best to manually specify sites?
    #NOTE: make sure learn_from_test_basins is set to FALSE. that param governs
    #some of the same logic in section 3a
    separate_pretrain_finetune_test = TRUE, #SPECIALIST?  and GENERALIST? and BASIC
    #alternative to finetune_to_test_basins. specify basins to finetune to.
    #set to NA to turn off (must turn off both to avoid finetuning)
    # finetune_basins = NA,
    finetune_basins = unlist(ms_basins_with_netcdfs, use.names = FALSE), #**** #is this for finetune 1?
    # finetune_basins = c('w1', 'w2', 'w3', 'w4', 'w5'), #, ****
    # finetune_basins = c( #126 total #****
    # 'GFCP', 'GFGB', 'GFGL', 'GFVN', 'POBR', 'DRKR', 'BARN', 'MCDN', 'MAWI', #baltimore
    # 'C2', 'C3', 'C4', #bonanza
    # 'BC_SW_4', 'BC_SW_20', #boulder
    # 'weir_4', #calhoun
    # 'RedondoMeadow', 'HistoryGrove', 'UpperJaramillo', 'LowerJaramillo', 'UpperLaJara', 'LowerLaJara', 'LowerRedondo', 'MarshallGulch', 'Bigelow', 'OracleRidge', 'UpperRedondo', 'FLUME_MCZOB', #catalina_jemez
    # 'PH', 'Copper', 'Avery', 'Marmot', 'Gothic', 'Rock', 'Bradley', 'Rustlers', 'EAQ', 'Quigley', #east_river
    # 'WS-5', 'WS-4', 'WS-2', 'WS-3', 'WS-10', 'WS-1', 'WS-13', 'WS-7', 'WS-6', #fernow
    # 'w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'w7', 'w8', 'w9', #hbef
    # 'GSLOOK', 'GSWS01', 'GSWS02', 'GSWS03', 'GSWS06', 'GSWS07', 'GSWS08', 'GSWS09', 'GSWS10', 'GSMACK', #hjandrews
    # 'N04D', 'N20B', 'N01B', 'N02B', #konza
    # 'P300', 'P301', 'P303', 'P304', 'D102', 'B200', 'B201', 'B203', 'B204', 'T003', #krew
    # 'Q1', 'Q2', 'Q3', 'QG', 'QS', 'RI', 'MPR', #luquillo
    # 'ALBION', 'GREEN4', 'GREEN5', 'MARTINELLI', 'NAVAJO', 'SADDLE', #niwot
    # # 'cart_creek', 'saw_mill_brook', 'bear_meadow', #plum
    # # 'AB00', 'AT07', 'BC02', 'DV01', 'GV01', 'HO00', 'MC00', 'ON02', 'RG01', 'RS02', 'FK00', 'TE03', 'SM04', 'CP00', 'SM01', 'RN01', #santa_barbara
    # 'WS77', 'WS78', 'WS79', 'WS80', #santee
    # 'GRO', 'SCAL', 'SCO', 'SH_weir', #suef
    # 'GSCC01', 'GSCC02', 'GSCC03', 'GSCC04', #krycklan
    # 'black_earth_creek', #usgs
    # 'east_fork', 'west_fork' #walker_branch
    # ),
    include_finetune_basins_in_pretrain = FALSE,
    #if FALSE, this ensures that test basins will not be included in finetuning,
    #even if they're specified in finetune_basins above
    include_test_basins_in_finetune = FALSE, #might be useful for SPECIALIST
    #what percent of the training/validation set should be comprised of large
    #   (> 2000 km^2) basins? Ignored unless basin_set == 'full'
    pct_large_basins = 10,

    # predictors <- 'precipitation'
    # static_attrs <- c('te_slope_mean', 'te_elev_mean', 'te_aspect_mean', 'hd_bfi_mean',
    #                   'lg_nlcd_ice_snow')

    #choose percent of dates to use for training
    #   the rest will be for validation if testing is done on separate basins.
    #   otherwise, the rest will be split between validation and test
    train_split_pct = 80, #IGNORED NOW
    #how many of the 671 basins to use for training and validation? (kratzert/newman pubs use just 531 that are < 2000km^2)
    #   this does not apply to fine-tuning.
    #ms without arctic = 185, ms + camels = 716
    #If using NHM data, this limit doesn't count toward NHM counterpart sites
    n_train_val_basins = 531,#531,#750,#716,#531,
    #should the test and validation periods be switcheroo'd for finetuning?
    reverse_finetune_val_test = FALSE,

    #if an integer (not NULL), model paramater specifications will be ignored (and chosen
    #   automatically to sample parameter options)
    # param_search_iters = NULL, #****
    param_search_iters = 1, #20
    #if ensemble is an integer, that number of models will be fit and their
    #   model objects combined (not implemented)
    # ensemble = NULL,
    ensemble = 10,

    #for finetuning to fill gaps, might as well learn from whatever data we have (this def overfits as implemented)
    finetune_predict_on_trainseq = FALSE,

    #NHM stuff
    #each site's netcdf has a column for nhm vs. true. let model see that?
    include_nhm_column = FALSE,
    #1 to train on [NHM-CAMELS, NHM-MS], then [CAMELS, MS], then finetune on NHM-target, and predict on target
    #2 to include nhm with the pretrain, continue, and ft phases of a SPECIALIST run (CAMELS+NHM pretrain, MS+NHM finetune, target NHM finetune, target eval)
    #NA to exclude nhm
    nhm_lineup = 2,
    # nhm_lineup = NA,
    #for running SPECIALISTS or GENERALISTS without NHM. note, you may have to adjust finetune and/or test
    #ranges manually in the config files #****?
    configure_specialist = ifelse(strategy == 'specialist', T, F),
    configure_generalist = ifelse(strategy == 'generalist', T, F), #start from the same camels pretrain each time

    #discharge and tmin as targets
    multitarget = TRUE,

    #something got messed up and separate_pretrain_finetune_test doesn't result in full dateranges for the test set, so this patch:
    # guarantee_full_test_daterange = TRUE #not built; see section 5b

    #specify the finetuning modules for the first finetune of a specialist run
    # continue_modules = list('head')
    continue_modules = NULL
)

if(addtl_config_list$include_finetune_basins_in_pretrain && ! addtl_config_list$finetune_predict_on_trainseq){
    stop(paste('only use include_finetune_basins_in_pretrain if also using finetune_predict_on_trainseq.',
               'this script is not yet set up to prevent overfitting otherwise. maybe that is not necessary?'))
}

if(is.na(addtl_config_list$finetune_basins[1]) && addtl_config_list$finetune_to_test_basins){
    addtl_config_list$finetune_basins <- c(addtl_config_list$camels_test_basins,
                                            unlist(addtl_config_list$ms_test_basins,
                                                   use.names = FALSE))
}

clst_type <- ifelse(.Platform$OS.type == 'windows', 'PSOCK', 'FORK')

## 0b helper functions ####

source(file.path(ms_dir, '..', 'src', 'output_dataset_convenience_functions',
                 'load_product.R'))

# extract_var_prefix <- function(x){
#
#     prefix <- substr(x, 1, 2)
#
#     return(prefix)
# }
#
# drop_var_prefix <- function(x){
#
#     unprefixed <- substr(x, 4, nchar(x))
#
#     return(unprefixed)
# }
#
# source(file.path(ms_dir, '..', 'src', 'dev', 'dev_helpers.R'))

# list_all_product_dirs <- function(prodname, location){
#
#     prodname_dirs <- list.dirs(path = location,
#                                full.names = TRUE,
#                                recursive = TRUE)
#
#     prodname_dirs <- grep(pattern = paste0('derived/', prodname, '__'),
#                           x = prodname_dirs,
#                           value = TRUE)
#
#     return(prodname_dirs)
# }

supplement_with_nhm <- function(x, nhm_list){

    if(is.null(x) || length(x) == 0) return(x)

    sups <- nhm_list[substr(nhm_list, 5, nchar(nhm_list)) %in% x]
    supped <- c(x, sups)
    return(supped)
}

## 1a load (or compile) start and end dates for each target and feature file (CAMELS) ####

q_files <- list.files(file.path(camels_dir, 'usgs_streamflow'),
                      recursive = TRUE,
                      full.names = TRUE)

forcing_files <- list.files(file.path(camels_dir, 'basin_mean_forcing/daymet'),
                            recursive = TRUE,
                            full.names = TRUE)

stored_dateranges_cmls <- file.path(config_dir,
                                    '..', '..', '..',
                                    'data', 'nh_methods', 'dateranges_cmls.csv')

if(file.exists(stored_dateranges_cmls)){

    message('using stored camels dateranges')
    dateranges_cmls <- read_csv(stored_dateranges_cmls)

} else {

    clst <- parallel::makeCluster(spec = parallel::detectCores(),
                                  type = clst_type)
    doParallel::registerDoParallel(clst)

    dateranges_cmls <- foreach(i = seq_along(q_files),
                          .combine = bind_rows) %dopar% {

        #might want to look for missing values that belie apparent temporal ranges

        q_file <- q_files[i]
        basin_id <- str_match(q_file,
                              'usgs_streamflow/[0-9]+/([0-9]+)_streamflow_qc.txt')[, 2]
        forcing_file <- grep(basin_id, forcing_files, value = TRUE)

        dateranges_q <- read.fwf(
            q_file,
            widths = c(9, 4, 3, 3, 9, 2),
            header = FALSE,
            strip.white = TRUE,
            colClasses = 'character',
            col.names = c('basin_id', 'year', 'month', 'day', 'Q', 'flag')
        ) %>%
            as_tibble() %>%
            mutate(date = ymd(paste(year, month, day))) %>%
            summarize(basin_id = first(basin_id),
                      mindate_q = min(date),
                      maxdate_q = max(date))

        dateranges_f  <- read.fwf(
            forcing_file,
            skip = 4,
            widths = c(5, 3, 3, 100),
            header = FALSE,
            strip.white = TRUE,
            colClasses = 'character',
            col.names = c('year', 'month', 'day', 'hour', 'a', 'b', 'c', 'd', 'e', 'f', 'g')
        ) %>%
            as_tibble() %>%
            mutate(date = ymd(paste(year, month, day))) %>%
            summarize(mindate_f = min(date),
                      maxdate_f = max(date))

        dateranges_cmls <- bind_cols(dateranges_q, dateranges_f)

        return(dateranges_cmls)
    }

    write_csv(dateranges_cmls, stored_dateranges_cmls)

    try(parallel::stopCluster(clst), silent = TRUE)
}

## 1b load (or compile) start and end dates for each target and feature file (MacroSheds) ####

stored_dateranges_ms <- file.path(config_dir,
                                  '..', '..', '..',
                                  'data', 'nh_methods', 'dateranges_ms.csv')

if(file.exists(stored_dateranges_ms)){

    message('using stored ms dateranges')
    dateranges_ms <- read_csv(stored_dateranges_ms) %>%
        filter(! basin_id %in% arctic_sites_to_drop)

} else {

    # q_data <- load_product(macrosheds_root = ms_dir,
    #                        prodname = 'discharge') %>%
    #     filter(substr(var, 1, 1) == 'I') #only sensor/weir Q measurements for now
    #
    # forcing_data <- tibble()
    # for(d in ms_domains_avail){
    #
    #     forcing_data <- read_feather(file.path(camels_dir,
    #                                            '..', '..',
    #                                            'ms_in_camels_format',
    #                                            d,
    #                                            'daymet_full_climate.feather')) %>%
    #         bind_rows(forcing_data)
    # }
    #
    # forcing_data <- forcing_data %>%
    #     filter(! site_code %in% setdiff(unique(forcing_data$site_code),
    #                                     unique(q_data$site_code)))
    #
    #ms_sites <- unique(forcing_data$site_code)

    ms_files <- list.files(file.path(combined_camels_ms_dir, 'time_series')) %>%
        str_subset('^NHM_', negate = TRUE) %>%
        str_subset('^[0-9]+\\.nc', negate = TRUE)
    ms_sites <- str_extract(ms_files, '(.+?)\\.nc$', 1)

    clst <- parallel::makeCluster(spec = parallel::detectCores(),
                                  type = clst_type)
    doParallel::registerDoParallel(clst)

    dateranges_ms <- foreach(i = seq_along(ms_sites),
                             .combine = bind_rows) %do% {

         s = ms_sites[i]
         xx <- ncdf4::nc_open(glue('~/git/macrosheds/papers/q_sim/in/lstm_data/time_series/{ss}.nc',
                                   ss = s))
         dd <- bind_cols(
             date = as.Date(ncdf4::ncvar_get(xx, 'date')),
             q = ncdf4::ncvar_get(xx, 'discharge'),
             f = ncdf4::ncvar_get(xx, 'dayl'))

         ncdf4::nc_close(xx)

         first_q <- Position(function(x) ! is.na(x), dd$q)
         last_q <- Position(function(x) ! is.na(x), dd$q, right = TRUE)
         first_f <- Position(function(x) ! is.na(x), dd$f)
         last_f <- Position(function(x) ! is.na(x), dd$f, right = TRUE)

         dde <- tibble(basin_id = s,
                       mindate_q = dd$date[first_q],
                       maxdate_q = dd$date[last_q],
                       mindate_f = dd$date[first_f],
                       maxdate_f = dd$date[last_f])

        return(dde)
        # dateranges_q <- q_data %>%
        #     filter(site_code == !!ms_sites[i]) %>%
        #     summarize(basin_id = first(site_code),
        #               mindate_q = as.Date(min(datetime)),
        #               maxdate_q = as.Date(max(datetime)))
        #
        # dateranges_f <- forcing_data %>%
        #     filter(site_code == !!ms_sites[i]) %>%
        #     summarize(mindate_f = as.Date(min(date)),
        #               maxdate_f = as.Date(max(date)))
        #
        # dateranges_comb <- bind_cols(dateranges_q, dateranges_f)
        #
        # return(dateranges_comb)
    }

    try(parallel::stopCluster(clst), silent = TRUE)

    # for(s in neon_sites){
    #
    #     xx <- ncdf4::nc_open(glue('~/git/macrosheds/papers/q_sim/in/lstm_data/time_series/{ss}_MANUALQ.nc',
    #                              ss = s))
    #     dd <- bind_cols(
    #         date = as.Date(ncdf4::ncvar_get(xx, 'date')),
    #         q = ncdf4::ncvar_get(xx, 'discharge'),
    #         f = ncdf4::ncvar_get(xx, 'dayl'))
    #
    #     ncdf4::nc_close(xx)
    #
    #     first_q <- Position(function(x) ! is.na(x), dd$q)
    #     last_q <- Position(function(x) ! is.na(x), dd$q, right = TRUE)
    #     first_f <- Position(function(x) ! is.na(x), dd$f)
    #     last_f <- Position(function(x) ! is.na(x), dd$f, right = TRUE)
    #
    #     dde <- tibble(basin_id = paste0(s, '_MANUALQ'),
    #                   mindate_q = dd$date[first_q],
    #                   maxdate_q = dd$date[last_q],
    #                   mindate_f = dd$date[first_f],
    #                   maxdate_f = dd$date[last_f])
    #
    #     dateranges_ms <- bind_rows(dateranges_ms, dde)
    # }

    # dateranges_ms = dateranges_ms %>%
    #     filter(basin_id %in% nghbs) %>%
    #     mutate(basin_id = paste0(basin_id, '_GAPPED')) %>%
    #     bind_rows(dateranges_ms) %>%
    #     distinct()
        # write_csv('../imputation/data/nh_methods/dateranges_ms.csv')

    write_csv(dateranges_ms, stored_dateranges_ms)
}


## 1c load (or compile) start and end dates for netcdf (NHM) ####

nhm_files <- list.files(file.path(combined_camels_ms_dir, 'time_series'),
                        recursive = TRUE,
                        full.names = TRUE,
                        pattern = '^NHM_')

stored_dateranges_nhm <- file.path(config_dir,
                                   '..', '..', '..',
                                   'data', 'nh_methods', 'dateranges_nhm.csv')

if(file.exists(stored_dateranges_nhm)){

    dateranges_nhm <- read_csv(stored_dateranges_nhm)

} else {

    clst <- parallel::makeCluster(spec = parallel::detectCores(),
                                  type = clst_type)
    doParallel::registerDoParallel(clst)

    dateranges_nhm <- foreach(
        i = seq_along(nhm_files),
        .combine = bind_rows) %dopar% {

            nhm_file <- nhm_files[i]
            basin_id <- str_match(nhm_file, '([^/]+)\\.nc')[, 2]

            nhm_con <- ncdf4::nc_open(nhm_file)
            nhm_d <- tibble(date = as.Date(ncdf4::ncvar_get(nhm_con, 'date')),
                            discharge = ncdf4::ncvar_get(nhm_con, 'discharge'),
                            prcp = ncdf4::ncvar_get(nhm_con, 'prcp'))
            ncdf4::nc_close(nhm_con)

            # if(any(is.na(nhm_d$prcp)) || nhm_d$prcp == -999) stop()
            # if(any(is.na(nhm_d$discharge)) || nhm_d$discharge == -999) stop()
            #no missing values, so no need to separately compute Q and forcing ranges

            dateranges_nhm <- nhm_d %>%
                summarize(mindate_q = min(date),
                          maxdate_q = max(date)) %>%
                mutate(mindate_f = mindate_q,
                       maxdate_f = maxdate_q,
                       basin_id = !!basin_id) %>%
                relocate(basin_id, .before = 'mindate_q')

            return(dateranges_nhm)
        }

    write_csv(dateranges_nhm, stored_dateranges_nhm)

    try(parallel::stopCluster(clst), silent = TRUE)
}


## 1d combine macrosheds, CAMELS, and NHM dateranges ####

attr_data <- read_csv(file.path(combined_camels_ms_dir,
                                'attributes',
                                'ms_attributes.csv'))

if(addtl_config_list$learn_from_test_basins) stop('test basins are being filtered just below. do they get added in again later?')
ms_trainval_basins <- ms_site_data %>%
    filter(site_type == 'stream_gauge') %>%
    select(domain, site_code) %>%
    filter(domain %in% addtl_config_list$ms_train_domains,
           site_code %in% dateranges_ms$basin_id,
           site_code %in% attr_data$site_code,
           ! site_code %in% unname(sapply(addtl_config_list$ms_test_basins,
                                          function(x) x[[1]]))) %>%
    pull(site_code)

ms_trainval_basins <- c(ms_trainval_basins, addtl_config_list$ms_train_basins)

extrap_basins <- grep('_extrapolate$', dateranges_ms$basin_id, value = TRUE)

dateranges_ms <- dateranges_ms %>%
    filter(basin_id %in% c(ms_trainval_basins,
                           unlist(addtl_config_list$ms_test_basins,
                                  use.names = FALSE),
                           addtl_config_list$finetune_basins,
                           extrap_basins))

dateranges_nhm_ms <- dateranges_nhm %>%
    filter(substr(basin_id, 5, nchar(basin_id)) %in% c(dateranges_ms$basin_id,
                                                       neon_sites))
dateranges_nhm_cmls <- dateranges_nhm %>%
    filter(substr(basin_id, 1, 4) == 'NHM_',
           substr(basin_id, 5, nchar(basin_id)) %in% dateranges_cmls$basin_id)

if(! is.na(addtl_config_list$nhm_lineup)){
    dateranges_ms <- bind_rows(dateranges_ms, dateranges_nhm_ms)
    dateranges_cmls <- bind_rows(dateranges_cmls, dateranges_nhm_cmls)
    dateranges <- bind_rows(dateranges_cmls, dateranges_ms)
} else {
    dateranges <- bind_rows(dateranges_cmls, dateranges_ms)
}

dateranges_copy <- dateranges

## 2a adjust a few things for NHM if desired ####
nhm_cmls_ids <- dateranges_nhm_cmls$basin_id
nhm_ms_ids <- dateranges_nhm_ms$basin_id
nhm_cmls_ids <- nhm_cmls_ids[substr(nhm_cmls_ids, 5, nchar(nhm_cmls_ids)) %in% kratzert_basin_list]

if(is.na(addtl_config_list$nhm_lineup)){
    nhm_cmls_ids <- NULL
    nhm_ms_ids <- NULL
}


## 2b investigate and clean up results ####

#identify typical temporal extents of Q data for each basin, and deviations (made more sense when only NEON sites from MS were included)
(late_starts_q <- filter(dateranges, mindate_q != min(dateranges$mindate_q)))
(early_ends_q <- filter(dateranges, maxdate_q != max(dateranges$maxdate_q)))
(both_q <- intersect(late_starts_q, early_ends_q))
(neither_q <- dateranges %>%
    select(basin_id, mindate_q, maxdate_q) %>%
    anti_join(union(late_starts_q, early_ends_q)))

#choose which collection of basins to subset from
all_camels_ids <- str_match(q_files,
                            'usgs_streamflow/[0-9]+/([0-9]+)_streamflow_qc.txt')[, 2]
all_camels_ids <- c(all_camels_ids, nhm_cmls_ids)
kratzert_basin_list <- c(kratzert_basin_list, nhm_cmls_ids)
omitted_basins <- setdiff(all_camels_ids, kratzert_basin_list)
setdiff(kratzert_basin_list, all_camels_ids) #empty is good
length(kratzert_basin_list) + length(omitted_basins) #674 is good, or 1200 if nhm included

if(addtl_config_list$basin_set == 'kratzert'){

    #ms_data gets added back on below

    #drop basins not included by Newman (2017) and Kratzert (2019)
    #   https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019WR026065
    dateranges <- filter(dateranges,
                         basin_id %in% kratzert_basin_list | grepl('^[0-9]+_GAPPED$', basin_id))
                                         # paste0(kratzert_basin_list, '_GAPPED')))

} else if(addtl_config_list$basin_set == 'custom'){

    if(! is.na(addtl_config_list$nhm_lineup)) stop('not set up for nhm')

    #ms_data gets added back on below
    custom_basin_list <- read_lines(custom_basin_list_loc)
    dateranges <- filter(dateranges,
                         basin_id %in% custom_basin_list)

} else if(addtl_config_list$basin_set != 'full'){
    if(! is.na(addtl_config_list$nhm_lineup)) stop('not set up for nhm')
    stop('addtl_config_list$basin_set must be "custom", "kratzert" or "full"')
}

if(addtl_config_list$include_ms_data){

    dateranges <- bind_rows(dateranges, dateranges_ms)
}

#for rows with different start dates for Q and forcings, adjust earlier
#   start date forward to match the later one
startdate_overhangs <- abs(dateranges$mindate_q - dateranges$mindate_f) > 200
startdate_overhangs[grepl('_extrapolate$', dateranges$basin_id)] <- FALSE
mindate_cols <- c('mindate_q', 'mindate_f')
for(sdo in which(startdate_overhangs)){
    overhang_row <- dateranges[sdo, ]
    later_start <- names(which.max(overhang_row[, mindate_cols]))
    earlier_start <- mindate_cols[mindate_cols != later_start]
    dateranges[sdo, earlier_start] <- dateranges[sdo, later_start] - 200
}

#for rows with different end dates for Q and forcings, adjust later
#   end date back to match the earlier one
enddate_overhangs <- dateranges$maxdate_q != dateranges$maxdate_f
enddate_overhangs[grepl('_extrapolate$', dateranges$basin_id)] <- FALSE
maxdate_cols <- c('maxdate_q', 'maxdate_f')
for(edo in which(enddate_overhangs)){
    overhang_row <- dateranges[edo, ]
    earlier_end <- names(which.min(overhang_row[, maxdate_cols]))
    later_end <- maxdate_cols[maxdate_cols != earlier_end]
    dateranges[edo, later_end] <- dateranges[edo, earlier_end]
}

#sanity check (both should be FALSE)
any(dateranges$mindate_q > dateranges$maxdate_q)
any(dateranges$mindate_f > dateranges$maxdate_f)

## 3a set up next model run ####

existing_runs <- list.dirs(config_dir,
                           recursive = FALSE,
                           full.names = FALSE)

combined_run_bool <- grepl('^runs', existing_runs)
combined_run_dirs <- existing_runs[combined_run_bool]
combined_run_bool <- combined_run_bool | grepl('^[a-zA-Z_]+$', existing_runs)
normal_run_dirs <- existing_runs[! combined_run_bool]

combined_run_ids <- lapply(combined_run_dirs, function(x){
    run_range <- as.numeric(str_match(x, 'runs_([0-9]+)-([0-9]+)')[, 2:3])
    return(run_range[1]:run_range[2])
}) %>%
    unlist()

if(ix == 1){
    if(length(normal_run_dirs)){

        if(length(combined_run_ids)){

            last_combined_runid <- max(combined_run_ids)
            last_normal_runid <- max(as.numeric(str_match(normal_run_dirs, 'run([0-9]+)')[, 2]))
            last_runid <- max(last_combined_runid, last_normal_runid)

        } else {
            last_runid <- max(as.numeric(str_match(normal_run_dirs, 'run([0-9]+)')[, 2]))
        }

    } else {
        last_runid <- 0
    }

    runid <- last_runid + 1

    doing_param_search <- ! is.null(addtl_config_list$param_search_iters)
    if(! doing_param_search){

        if(run_prediction_routine_across_sites) stop('unaccounted for')

        dir.create(file.path(config_dir, glue('run{ri}', ri = runid)))

        n_iters <- 1
        run_seq <- seq(runid, runid + n_iters - 1)

    } else {

        if(run_prediction_routine_across_sites && ! new_generalist){

            #this gets overridden
            n_mega_iters <- length(sites_for_routine)
            n_iters <- addtl_config_list$param_search_iters * n_mega_iters
            search_params_pretrain <- list(
                validate_every = c(1L),
                validate_n_random_basins = c(1000L),
                save_weights_every = 1L,
                metrics = list('NSE'),#, 'KGE', 'pbias'),
                output_activation = c('linear'),
                hidden_size = c(30L),
                initial_forget_bias = c(5L),
                output_dropout = c(0.5), #0.2
                loss = c('NSE'),
                clip_gradient_norm = c(5L),
                predict_last_n = c(1L),
                learning_rate = list(
                    # list(a = '1e-4', b = '5e-5', c = '1e-5')
                    list(a = '5e-4', b = '1e-4', c = '5e-5')
                ),
                batch_size = 512L,#c(1024L), #256
                epochs = c(30L), #40
                log_tensorboard = 'True',
                log_n_figures = 'False',
                save_validation_results = 'True',
                statics_embedding = list(
                    list(type = 'fc', hiddens = list(20L), activation = 'tanh', dropout = 0.4)
                ),
                dynamics_embedding = list(
                    list(type = 'fc', hiddens = list(200L), dropout = 0.4)
                ))

        } else {
            n_mega_iters <- 1
            n_iters <- addtl_config_list$param_search_iters
        }
        run_seq <- seq(runid, runid + n_iters - 1)

        runset_parent_dir <- file.path(config_dir, glue('runs_{a}-{b}',
                                                        a = run_seq[1],
                                                        b = run_seq[length(run_seq)]))
        dir.create(runset_parent_dir)
        lapply(file.path(runset_parent_dir, paste0('run', run_seq)),
               dir.create)

        # #uncomment for focused param search or even duplicate runs (old)
        # search_params <- list(validate_every = c(1L),
        #                       validate_n_random_basins = c(1000L),
        #                       save_weights_every = 1L,
        #                       metrics = list('NSE'),
        #                       output_activation = c('linear'),
        #                       hidden_size = c(30L),
        #                       initial_forget_bias = c(5L),
        #                       output_dropout = c(0.4),
        #                       # output_dropout = c(0.2),
        #                       loss = c('NSE'),
        #                       clip_gradient_norm = c(5L),
        #                       predict_last_n = c(1L),
        #                       # seq_length
        #                       learning_rate = list(
        #                                            # list(a = '1e-3',  #.001
        #                                            #      b = '5e-4',  #.0005
        #                                            #      c = '1e-4'), #.0001
        #                                            # list(a = '5e-3',  #.005
        #                                            #      b = '1e-3',  #.001
        #                                            #      c = '5e-4'), #.0005
        #                                            # list(a = '5e-4',  #.0005
        #                                            #      b = '1e-4',  #.0001
        #                                            #      c = '5e-5'),  #.00005
        #                                            # list(a = '1e-4',  #.0001
        #                                            #      b = '5e-5',  #.00005
        #                                            #      c = '1e-5'), #.00001
        #                                            list(a = '5e-5',  #.00005
        #                                                 b = '1e-5',  #.00001
        #                                                 c = '5e-6') #.000005
        #                                            # list(a = '5e-6',  #.000005
        #                                            #      b = '1e-6',  #.000001
        #                                            #      c = '5e-7') #.0000005
        #                                            # list(a = '5e-7',  #.0000005
        #                                            #      b = '1e-7',  #.0000001
        #                                            #      c = '5e-8')  #.00000005
        #                       ),
        #                       batch_size = c(64L),
        #                       # batch_size = c(256L), #1024L
        #                       epochs = c(40L),
        #                       log_tensorboard = 'True',
        #                       log_n_figures = 'False',
        #                       save_validation_results = 'True',
        #                       # statics_embedding = list(
        #                       #     list(type = 'fc', hiddens = list(10L), dropout = 0.4),
        #                       #     list(type = 'fc', hiddens = list(50L), dropout = 0.4)
        #                       # ),
        #                       # dynamics_embedding = list(
        #                       #     list(type = 'fc', hiddens = list(50L, 50L, 50L), dropout = 0.4),
        #                       #     list(type = 'fc', hiddens = list(100L, 50L), dropout = 0.4),
        #                       #     list(type = 'fc', hiddens = list(200L), dropout = 0.4)
        #                       # ),
        #                       #comment the following 3 if not finetuning
        #                       # learning_rate2 = list(list(a = '5e-2',
        #                       #                            b = '1e-2')),
        #                       epochs2 = c(20L),
        #                       # epochs2 = c(10L, 20L, 30L),
        #                       # finetune_modules = list(list('lstm')))
        #                       finetune_modules = list(list('lstm')))

        # #uncomment for focused param search
        # search_params <- list(validate_every = c(1L),
        #                       validate_n_random_basins = c(1000L),
        #                       save_weights_every = 1L,
        #                       metrics = list('NSE'),
        #                       output_activation = c('linear'),
        #                       hidden_size = c(30L),
        #                       initial_forget_bias = c(5L),
        #                       output_dropout = 0.4,
        #                       loss = c('NSE'),
        #                       clip_gradient_norm = c(5L),
        #                       predict_last_n = c(1L),
        #                       learning_rate = list(
        #                                            list(a = '1e-4',  #.0001
        #                                                 b = '5e-5',  #.00005
        #                                                 c = '1e-5') #.00001
        #                       ),
        #                       batch_size = c(128L),
        #                       epochs = c(30L),
        #                       log_tensorboard = 'True',
        #                       log_n_figures = 'False',
        #                       save_validation_results = 'True',
        #                       # statics_embedding = list(
        #                       #     list(type = 'fc', hiddens = list(10L), dropout = 0.4),
        #                       #     list(type = 'fc', hiddens = list(50L), dropout = 0.4)
        #                       # ),
        #                       # dynamics_embedding = list(
        #                       #     list(type = 'fc', hiddens = list(50L, 50L, 50L), dropout = 0.4),
        #                       #     list(type = 'fc', hiddens = list(100L, 50L), dropout = 0.4),
        #                       #     list(type = 'fc', hiddens = list(200L), dropout = 0.4)
        #                       # ),
        #                       #comment the following 3 if not finetuning
        #                       # learning_rate2 = list(list(a = '5e-2',
        #                       #                            b = '1e-2')),
        #                       epochs2 = c(10L),
        #                       # finetune_modules = list(list('lstm')))
        #                       finetune_modules = list(list('head', 'lstm')))

        # #uncomment for targeted param search (i.e. not universally good params, but good for one run)
        # search_params <- list(validate_every = c(1L),
        #                       validate_n_random_basins = c(1000L),
        #                       save_weights_every = 1L,
        #                       metrics = list('NSE'),
        #                       output_activation = c('linear'),
        #                       hidden_size = c(30L),
        #                       initial_forget_bias = c(5L),
        #                       output_dropout = 0.4,
        #                       loss = c('NSE'),
        #                       clip_gradient_norm = c(5L),
        #                       predict_last_n = c(1L),
        #                       learning_rate = list(
        #                                            list(a = '1e-4',  #.0001
        #                                                 b = '5e-5',  #.00005
        #                                                 c = '1e-5') #.00001
        #                       ),
        #                       batch_size = c(64L),
        #                       epochs = c(30L),
        #                       log_tensorboard = 'True',
        #                       log_n_figures = 'False',
        #                       save_validation_results = 'True',
        #                       epochs2 = c(30L),
        #                       finetune_modules = list(list('lstm')))
        #                       # finetune_modules = list(list('head', 'lstm')))

        #uncomment for wide param search
        search_params <- list(validate_every = c(1L),
                              validate_n_random_basins = c(1000L),
                              save_weights_every = 1L,
                              metrics = list('NSE'),# 'KGE', 'pbias'),
                              output_activation = c('linear'),
                              # hidden_size = hidden_size_,
                              hidden_size = 30L,
                              # hidden_size = c(20L, 30L, 40L),
                              initial_forget_bias = c(5L),
                              # output_dropout = output_dropout_,
                              output_dropout = c(0.2, 0.3, 0.4, 0.5),
                              loss = c('NSE'),
                              clip_gradient_norm = c(5L),
                              predict_last_n = c(1L),
                              # learning_rate = learning_rate_,
                              learning_rate = list(
                                  list(a = '1e-2',
                                       b = '1e-3',
                                       c = '1e-4'),
                                  list(a = '1e-3',
                                       b = '1e-4',
                                       c = '1e-5'),
                                  list(a = '1e-4',
                                       b = '1e-5',
                                       c = '1e-6'),
                                  list(a = '1e-5',
                                       b = '1e-6',
                                       c = '1e-7')
                              ),
                              batch_size = c(32L, 64L, 128L, 256L, 512L),
                              # batch_size = batch_size_,
                              epochs = c(20L, 30L, 40L),
                              # epochs = epochs_,
                              log_tensorboard = 'True',
                              log_n_figures = 'False',
                              save_validation_results = 'True',
                              # statics_embedding = list(
                              #     list(type = 'fc', hiddens = list(10L), dropout = 0.4),
                              #     list(type = 'fc', hiddens = list(50L), dropout = 0.4)
                              # ),
                              # dynamics_embedding = list(
                              #     list(type = 'fc', hiddens = list(50L, 50L, 50L), dropout = 0.4),
                              #     list(type = 'fc', hiddens = list(100L, 50L), dropout = 0.4),
                              #     list(type = 'fc', hiddens = list(200L), dropout = 0.4)
                              # ),
                              #comment the following 3 if not finetuning
                              # learning_rate2 = list(list(a = '5e-2',
                              #                            b = '1e-2')),
                              epochs2 = c(10L, 20L, 30L), ###
                              # epochs2 = c(20L), ###
                              # finetune_modules = list(list('lstm'))
                              # finetune_modules = list(list('head', 'lstm'))
                              # finetune_modules = finetune_modules_
                              finetune_modules = list(list('head'), list('lstm'), list('head', 'lstm'))
        )

        # #uncomment for wide param search (old)
        # search_params <- list(validate_every = c(1L, 3L),
        #                       validate_n_random_basins = c(5L, 25L),
        #                       save_weights_every = 1L,
        #                       metrics = list('NSE'),
        #                       # activation = c('tanh', 'sigmoid', 'linear'), #hidden layer activation, needs to be a dict
        #                       output_activation = c('linear', 'softplus'),
        #                       hidden_size = c(20L, 30L, 40L, 50L),
        #                       initial_forget_bias = c(2L, 3L, 5L),
        #                       output_dropout = c(0.1, 0.2, 0.4, 0.6),
        #                       loss = c('NSE', 'MSE', 'RMSE'),
        #                       clip_gradient_norm = c(1L, 5L, 10L),
        #                       predict_last_n = c(1L),# 30L, 365L),
        #                       # seq_length
        #                       learning_rate = list(
        #                                            list(a = '1e-3',  #.001
        #                                                 b = '5e-4',  #.0005
        #                                                 c = '1e-4'), #.0001
        #                                            list(a = '5e-3',  #.005
        #                                                 b = '1e-3',  #.001
        #                                                 c = '5e-4'), #.0005
        #                                            list(a = '5e-4',  #.0005
        #                                                 b = '1e-4',  #.0001
        #                                                 c = '5e-5') #.00005
        #                                            # list(a = '1e-4',  #.0001
        #                                            #      b = '5e-5',  #.00005
        #                                            #      c = '1e-5'), #.00001
        #                                            # list(a = '5e-5',  #.00005
        #                                            #      b = '1e-5',  #.00001
        #                                            #      c = '5e-6'), #.000005
        #                                            # list(a = '5e-6',  #.000005
        #                                            #      b = '1e-6',  #.000001
        #                                            #      c = '5e-7'), #.0000005
        #                                            # list(a = '5e-7',  #.0000005
        #                                            #      b = '1e-7',  #.0000001
        #                                            #      c = '5e-8')#.00000005
        #                       ),
        #                       batch_size = c(128L, 512L),
        #                       # batch_size = c(512L, 1024L, 2048L),
        #                       epochs = c(20L, 30L, 40L),
        #                       log_tensorboard = 'True',
        #                       log_n_figures = 'False',
        #                       save_validation_results = 'False',
        #                       #comment the following if not finetuning
        #                       # learning_rate2 = list(list(a = '5e-4',  #.0005
        #                       #                            b = '1e-4'), #.0001
        #                       #                       list(a = '1e-4',  #.0001
        #                       #                            b = '5e-5'), #.00005
        #                       #                       list(a = '5e-3',  #.005
        #                       #                            b = '1e-3'), #.001
        #                       #                       list(a = '5e-6',  #.000005
        #                       #                            b = '1e-6'), #.000001
        #                       #                       list(a = '5e-2',  #.05
        #                       #                            b = '1e-2')),#.01
        #                       epochs2 = c(20L, 30L, 40L),
        #                       finetune_modules = list(list('lstm'),
        #                                               list('head', 'lstm'))
        #                       # #if using transformer, uncomment the following:
        #                       # ,
        #                       # model = 'transformer',
        #                       # transformer_nlayers = c(1L, 2L, 4L, 8L),
        #                       # transformer_positional_encoding_type = list('sum', 'concatenate'),
        #                       # transformer_dim_feedforward = c(1L, 2L, 4L, 8L),
        #                       # transformer_positional_dropout = c(0.1, 0.2, 0.4, 0.6),
        #                       # transformer_dropout = c(0.1, 0.2, 0.4, 0.6),
        #                       # transformer_nhead = c(1L, 2L, 3L)
        #                       )

        n_param_options_per <- vapply(search_params, length, 1)
        n_param_options <- sum(n_param_options_per)

        run_param_list <- lapply(search_params,
                                 function(x){
                                     if(length(x) == 1){
                                         rep(x, n_iters)
                                     } else {
                                         sample(x, size = n_iters, replace = TRUE)
                                     }
                                 })

        run_param_list$learning_rate2 <- lapply(run_param_list$learning_rate,
                                               function(x){
                                                   xx <- x[2:3]
                                                   names(xx) <- c('a', 'b')
                                                   return(xx)
                                                })

        run_param_list <- lapply(1:n_iters, function(x){
            lapply(run_param_list, function(y){
                y[x]
            })
        })
    }
}

## 3b create new config file with desired specs (REQUIRES INPUT unless parameter searching) ####

cfg <- list()

# cfg$experiment_name <- glue('run', runid)

if(! doing_param_search){
    tbf <- glue('{cdr}/run{ri}/run{ri}_train_validate.txt',
                cdr = config_dir,
                ri = runid)
    pbtpf <- glue('{cdr}/run{ri}/train_ranges.pkl',
                  cdr = config_dir,
                  ri = runid)
} else {
    tbf <- glue('{rpd}/train_validate.txt',
                rpd = runset_parent_dir)
    pbtpf <- glue('{rpd}/train_ranges.pkl',
                  rpd = runset_parent_dir)
}

cfg$train_basin_file <- tbf
cfg$validation_basin_file <- cfg$train_basin_file
cfg$test_basin_file <- gsub('train_validate', 'test', cfg$train_basin_file)
cfg$per_basin_train_periods_file <- pbtpf
cfg$per_basin_validation_periods_file <- gsub('train_ranges',
                                              'validation_ranges',
                                              cfg$per_basin_train_periods_file)
cfg$per_basin_test_periods_file <- gsub('train_ranges',
                                        'test_ranges',
                                        cfg$per_basin_train_periods_file)
if(set_seeds) cfg$seed <- 1L
cfg$device <- 'cuda:0'
cfg$validate_every <- 3L
cfg$validate_n_random_basins <- 1000L #1L
cfg$metrics <- list('MSE') #NSE
cfg$model <- 'cudalstm'
# cfg$model <- 'transformer'
# if(cfg$model == 'transformer'){...
cfg$head <- 'regression'
cfg$output_activation <- 'linear'
cfg$hidden_size <- 30L #40L#pub=256L; takes longer; probably better; uses a bit more RAM
cfg$initial_forget_bias <- 5L#3L
cfg$output_dropout <- 0.4
cfg$optimizer <- 'Adam'
cfg$loss <- 'MSE' #'NSE' in pub
# cfg$learning_rate <- list(`0` = '5e-2',
#                           `10` = '1e-3',
#                           `20` = '5e-3')
# cfg$learning_rate <- list('1e-3') #pub? see line 377 of main.py for contraditcion
cfg$learning_rate <- list(`0` = '5e-3',
                          `10` = '5e-4',
                          `20` = '5e-5')
# cfg$learning_rate <- list(`0` = '1e-2',
#                           `30` = '5e-3',
#                           `40` = '1e-3')
cfg$batch_size <- 256L #2000L in the paper; too large and this overloads GPU
cfg$epochs <- 30L
cfg$clip_gradient_norm <- 1L
cfg$predict_last_n <- 1L
cfg$seq_length <- 200L#365L #pub 270
cfg$num_workers <- as.integer(parallel::detectCores())
cfg$log_interval <- 1L
cfg$log_tensorboard <- 'True'
cfg$save_validation_results <- 'True'
cfg$log_n_figures <- 1L
cfg$save_weights_every <- 10L
cfg$dataset <- ifelse(addtl_config_list$include_ms_data, 'generic', 'camels_us')
cfg$data_dir <- ifelse(addtl_config_list$include_ms_data,
                       combined_camels_ms_dir,
                       camels_dir)
cfg$allow_subsequent_nan_losses <- 2L
cfg$statics_embedding <- list(type = 'fc',
                              hiddens = list(20L),
                              activation = 'tanh',
                              dropout = 0.4)
# cfg$dynamics_embedding <- list(type = 'fc',
#                                hiddens = list(100L, 100L, 50L),
#                                activation = 'tanh',
#                                dropout = 0.4)
cfg$dynamics_embedding = list(type = 'fc', hiddens = list(200L), activation = 'tanh', dropout = 0.4)
# cfg$mc_dropout = 'True'
# cfg$n_samples = 10L
# cfg$negative_sample_handling = 'clip'

cfg$forcings <- if(addtl_config_list$include_ms_data) NULL else list('nldas')
# cfg$forcings <- list('maurer', 'nldas', 'daymet')
if(length(cfg$forcings) == 1 && cfg$forcings[[1]] %in% c('maurer', 'nldas')){
    cfg$dynamic_inputs <- list('PRCP(mm/day)', 'SRAD(W/m2)', 'Tmax(C)', 'Tmin(C)', 'Vp(Pa)')
} else if(is.null(cfg$forcings) || length(cfg$forcings) == 1){
    # cfg$dynamic_inputs <- list('prcp')

    if(addtl_config_list$include_nhm_column){
        cfg$dynamic_inputs <- list('prcp', 'srad', 'tmax', 'tmin', 'vp', 'dayl', 'q_source', 'pet', 'swe')
    } else {
        if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2){
            cfg$dynamic_inputs <- list('prcp', 'srad', 'tmax', 'tmin', 'vp', 'dayl', 'pet', 'swe') #used to be missing pet. now good.
        } else {
            cfg$dynamic_inputs <- list('prcp', 'srad', 'tmax', 'tmin', 'vp', 'dayl', 'pet', 'swe')#,
        }
                                   # 'QmaExtraNA', 'QmaExtraNAind') #Qbin0, QalternatingNA, QrandomNA
        # variable_features = 'existent'
        # if(ix %% 10 == 1){
        #     cfg$dynamic_inputs <- list('prcp', 'srad', 'tmax', 'tmin', 'vp', 'dayl', 'pet', 'swe')
        # } else if(ix %% 10 == 2){
        #     cfg$dynamic_inputs <- list('prcp', 'srad', 'tmax', 'tmin', 'vp', 'dayl', 'pet', 'swe', 'QalternatingNA')
        # } else if(ix %% 10 == 3){
        #     cfg$dynamic_inputs <- list('prcp', 'srad', 'tmax', 'tmin', 'vp', 'dayl', 'pet', 'swe', 'Qlag0')
        # } else if(ix %% 10 == 4){
        #     cfg$dynamic_inputs <- list('prcp', 'srad', 'tmax', 'tmin', 'vp', 'dayl', 'pet', 'swe', 'Qlag0', 'Qbin0')
        # } else if(ix %% 10 == 5){
        #     cfg$dynamic_inputs <- list('prcp', 'srad', 'tmax', 'tmin', 'vp', 'dayl', 'pet', 'swe', 'Qlag1')
        # } else if(ix %% 10 == 6){
        #     cfg$dynamic_inputs <- list('prcp', 'srad', 'tmax', 'tmin', 'vp', 'dayl', 'pet', 'swe', 'Qlag1', 'Qbin1')
        # } else if(ix %% 10 == 7){
        #     cfg$dynamic_inputs <- list('prcp', 'srad', 'tmax', 'tmin', 'vp', 'dayl', 'pet', 'swe', 'Qma', 'QmaNAind')
        # } else if(ix %% 10 == 8){
        #     cfg$dynamic_inputs <- list('prcp', 'srad', 'tmax', 'tmin', 'vp', 'dayl', 'pet', 'swe', 'QmaExtraNA', 'QmaExtraNAind')
        # } else if(ix %% 10 == 9){
        #     cfg$dynamic_inputs <- list('prcp', 'srad', 'tmax', 'tmin', 'vp', 'dayl', 'pet', 'swe', 'QrandomNA')
        # } else if(ix %% 10 == 0){
        #     cfg$dynamic_inputs <- list('prcp', 'srad', 'tmax', 'tmin', 'vp', 'dayl', 'pet', 'swe', 'QrandomNA', 'QrandomNAind')
        # }
    }

    if(addtl_config_list$multitarget){
        cfg$dynamic_inputs[cfg$dynamic_inputs == 'tmin'] <- NULL
    }

} else if(length(cfg$forcings) == 3){
    cfg$dynamic_inputs <-
        list('PRCP(mm/day)_maurer', 'SRAD(W/m2)_maurer', 'Tmax(C)_maurer', 'Tmin(C)_maurer', 'Vp(Pa)_maurer',
             'PRCP(mm/day)_nldas', 'SRAD(W/m2)_nldas', 'Tmax(C)_nldas', 'Tmin(C)_nldas', 'Vp(Pa)_nldas',
             'prcp(mm/day)_daymet', 'srad(W/m2)_daymet', 'tmax(C)_daymet', 'tmin(C)_daymet', 'vp(Pa)_daymet')
} else {
    stop('address this case')
}
# cfg$dynamic_inputs <- list('PRCP(mm/day)_nldas', 'SRAD(W/m2)_nldas', 'Tmax(C)_nldas', 'Tmin(C)_nldas', 'Vp(Pa)_nldas',
#                            'PRCP(mm/day)_maurer', 'SRAD(W/m2)_maurer', 'Tmax(C)_maurer', 'Tmin(C)_maurer', 'Vp(Pa)_maurer',
#                            'prcp(mm/day)_daymet', 'srad(W/m2)_daymet', 'tmax(C)_daymet', 'tmin(C)_daymet', 'vp(Pa)_daymet')

if(addtl_config_list$multitarget){
    cfg$target_variables <- list('discharge', 'tmin')
} else {
    cfg$target_variables <- if(addtl_config_list$include_ms_data) list('discharge') else list('QObs(mm/d)')
}
cfg$clip_targets_to_zero <- cfg$target_variables[1]

# invalid_attr <- list( #list of attributes deemed unusable in pub
#     'gauge_name', 'area_geospa_fabric', 'geol_1st_class', 'glim_1st_class_frac', 'geol_2nd_class',
#     'glim_2nd_class_frac', 'dom_land_cover_frac', 'dom_land_cover', 'high_prec_timing',
#     'low_prec_timing', 'huc', 'q_mean', 'runoff_ratio', 'stream_elas', 'slope_fdc',
#     'baseflow_index', 'hfd_mean', 'q5', 'q95', 'high_q_freq', 'high_q_dur', 'low_q_freq',
#     'low_q_dur', 'zero_q_freq', 'geol_porostiy', 'root_depth_50', 'root_depth_99', 'organic_frac',
#     'water_frac', 'other_frac'
# )

# cfg$static_attributes <- list(#this is the list from the finetuning example
#     'elev_mean' 'slope_mean' 'area_gages2' 'frac_forest' 'lai_max' 'lai_diff'
#     'gvf_max' 'gvf_diff' 'soil_depth_pelletier' 'soil_depth_statsgo' 'soil_porosity'
#     'soil_conductivity' 'max_water_content' 'sand_frac' 'silt_frac' 'clay_frac'
#     'carbonate_rocks_frac' 'geol_permeability' 'p_mean' 'pet_mean' 'aridity'
#     'frac_snow' 'high_prec_freq' 'high_prec_dur' 'low_prec_freq' 'low_prec_dur'
# )

#list Spencer grabbed
# [1] "p_mean"              "pet_mean"            "aridity"             "p_seasonality"
# [5] "frac_snow_daily"     "high_prec_freq"      "high_prec_dur"       "high_prec_timing"
# [9] "low_prec_freq"       "low_prec_dur"        "low_prec_timing"     "geol_1st_class"
# [13] "glim_1st_class_frac" "geol_2nd_class"      "glim_2nd_class_frac" "carb_rocks_frac"
# [17] "geol_porosity"       "geol_permeability"   "sand_frac"           "silt_frac"
# [21] "clay_frac"           "organic_frac"        "gauge_lat"           "gauge_lon"
# [25] "area"                "elev_mean"           "slope_mean"          "forest_frac"
# [29] "dom_land_cover_frac" "dom_land_cover"      "root_depth_50"       "root_depth_99"

#lines with # at the end might contain attrs not used in kratzert 2019
# cfg$static_attributes <- NULL
cfg$static_attributes <- list(

    # THESE ARE CHILL
    # 'geol_permeability', #used to be chill
    'p_mean', 'pet_mean','p_seasonality','aridity',
    'high_prec_freq','high_prec_dur', 'low_prec_freq', 'low_prec_dur',
    'elev_mean',
    # 'sand_frac','silt_frac','clay_frac',
    'gauge_lat', #
    # 'gauge_lon', #
    'slope_mean', #

    ## THESE NAMES HAVE BEEN CORRECTED
    # 'carbonate_rocks_frac', #(carb_rocks_frac) #used to be chill
    'frac_snow', #(frac_snow_daily)
    #'frac_forest', #(forest_frac) #ONLY NaN FOR arctic, luquillo; questionable? see Spencer dm on 9/8/21

    ifelse(addtl_config_list$include_ms_data, 'area', 'area_gages2') #("area_gages2" has been changed to "area" in the CAMELS set)

    ## THESE HAVE ISSUES OR ARE MISSING FOR MS SITES
    #'soil_depth_pelletier','soil_depth_statsgo','max_water_content','soil_conductivity'(unavailable?)
    #'lai_max','gvf_max','lai_diff','gvf_diff',(hard to get. see Spencer dm on 9/8/21)#
    # 'soil_porosity', (?)

    ## TRY THESE AT SOME POINT (not included in 2019 pub, but available)
    # 'high_prec_timing',
    # 'low_prec_timing', #problematic with np.isnan. would need to be encoded; missing for SH_weir
    # 'geol_1st_class','glim_1st_class_frac', #same deal
    # 'geol_2nd_class', 75 MS sites have NA
    # 'glim_2nd_class_frac',

    # 'organic_frac' missing for all the same sites that are missing sand, silt, clay

    # 'dom_land_cover_frac', 'dom_land_cover' #same deal again, also missing for arctic and luq
    # 'root_depth_50','root_depth_99'
)

## 3c write config and basin list files ####

if(! doing_param_search){
    if(run_prediction_routine_across_sites) stop('unaccounted for')
    path_extra <- ''
    n_configs <- 1
} else {
    path_extra <- glue('runs_{a}-{b}/',
                       a = run_seq[1],
                       b = run_seq[length(run_seq)])

    n_configs <- addtl_config_list$param_search_iters
}

if(run_prediction_routine_across_sites & ! new_generalist){
    cfg_pretrain <- cfg
    cfg_pretrain$experiment_name <- gsub('/', '_pretrain', path_extra)
    for(j in seq_along(search_params_pretrain)){
        eps <- search_params_pretrain$epochs
        names(search_params_pretrain$learning_rate[[1]]) <- floor(seq(0, eps, length.out = 4))[1:3]

        p <- search_params_pretrain[j]
        paramj <- names(p)
        if(paramj == 'learning_rate'){
            cfg_pretrain[[paramj]] <- p[[1]][[1]]
        } else if(paramj == 'metrics'){
            cfg_pretrain[[paramj]] <- as.list(p[[1]][[1]])
        } else if(! paramj %in% c('learning_rate2', 'epochs2', 'finetune_modules')){
            cfg_pretrain[[paramj]] <- p[[1]]
        }

        # cfg_writepath <- file.path(config_dir,
        #                            glue('{ex}{ex2}_pretrain.yml',
        #                                 ex = path_extra,
        #                                 ex2 = gsub('/', '', path_extra)))
        cfg_writepath <- file.path(config_dir,
                                   glue('{ex}pretrain.yml',
                                        ex = path_extra))
        write_yaml(cfg_pretrain, cfg_writepath)
        cfg_str <- read_file(cfg_writepath) %>% str_replace_all("'", '')
        cfg_str <- gsub('- type', '  type', cfg_str)
        write_file(cfg_str, cfg_writepath)
    }
}


for(ii in seq_len(n_configs)){
    cfg0 <- cfg
    i <- ix * n_configs + ii - n_configs
    runid <- run_seq[i]

    cfg0$experiment_name <- paste0('run', run_seq[i])

    if(doing_param_search){

        run_params <- run_param_list[[i]]
        for(j in seq_along(run_params)){

            #finalize learning rate steps based on the number of training epochs
            eps <- run_params$epochs
            names(run_params$learning_rate[[1]]) <- floor(seq(0, eps, length.out = 4))[1:3]

            p <- run_params[j]
            paramj <- names(p)
            if(paramj == 'learning_rate'){
                cfg0[[paramj]] <- p[[1]][[1]]
            } else if(paramj == 'metrics'){
                cfg0[[paramj]] <- as.list(p[[1]][[1]])
            } else if(! paramj %in% c('learning_rate2', 'epochs2', 'finetune_modules')){
                cfg0[[paramj]] <- p[[1]]
            }
        }
    }

    cfg_writepath <- file.path(config_dir,
                               glue('{ex}run{cr}/{type}{cr}.yml',
                                    ex = path_extra,
                                    type = ifelse(run_prediction_routine_across_sites,
                                                  'continue', 'run'),
                                    cr = runid))
    if(run_prediction_routine_across_sites){
        cfg0$train_basin_file <- sub('/train_validate', glue('/{xx}/continue', xx = cfg0$experiment_name), cfg0$train_basin_file)
        cfg0$validation_basin_file <- cfg0$train_basin_file
        cfg0$test_basin_file <- sub('continue.txt', 'test.txt', cfg0$train_basin_file)
        cfg0$per_basin_train_periods_file <- sub('/train', glue('/{xx}/continue_train', xx = cfg0$experiment_name), cfg0$per_basin_train_periods_file)
        cfg0$per_basin_validation_periods_file <- sub('/validation', glue('/{xx}/continue_validation', xx = cfg0$experiment_name), cfg0$per_basin_validation_periods_file)
        cfg0$per_basin_test_periods_file <- sub('/test', glue('/{xx}/test', xx = cfg0$experiment_name), cfg0$per_basin_test_periods_file)
        if(! is.null(addtl_config_list$continue_modules)){
            cfg0$finetune_modules = addtl_config_list$continue_modules[[1]]
        } else {
            cfg0$finetune_modules = list('head', 'lstm')
        }
    }

    write_yaml(cfg0, cfg_writepath)

    if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 1){
        nhm_cfg <- cfg0
        nhm_cfg$train_basin_file <- sub('/train', '/nhm_train', nhm_cfg$train_basin_file)
        nhm_cfg$validation_basin_file <- sub('/train', '/nhm_train', nhm_cfg$train_basin_file)
        # nhm_cfg$test_basin_file <- cfg0$test_basin_file
        nhm_cfg$per_basin_train_periods_file <- sub('/train', '/nhm_train', nhm_cfg$per_basin_train_periods_file)
        nhm_cfg$per_basin_validation_periods_file <- sub('/validation', '/nhm_validation', nhm_cfg$per_basin_validation_periods_file)
        # nhm_cfg$per_basin_test_periods_file <- cfg0$per_basin_test_periods_file
        write_yaml(nhm_cfg, sub('\\.yml$', 'nhm.yml', cfg_writepath))
    }
}

#save notes file containing "additional config" from the setup section
if(doing_param_search){
    notes_writepath <- file.path(config_dir,
                                 glue('{ex}addtl_notes.yml',
                                      ex = path_extra,
                                      cr = runid))
} else {
    notes_writepath <- file.path(config_dir,
                                 glue('run{cr}/addtl_notes.yml',
                                      cr = runid))
}

write_yaml(addtl_config_list, notes_writepath)

#write basin list files
if(addtl_config_list$basin_set == 'kratzert'){

    # if(! is.null(addtl_config_list$camels_test_basins)) stop('not set up for nhm here')
    if(! is.null(addtl_config_list$camels_test_basins) && ! is.na(addtl_config_list$nhm_lineup)) stop('not set up for nhm here')

    train_validate_basins <- setdiff(kratzert_basin_list, addtl_config_list$camels_test_basins)

    if(addtl_config_list$include_ms_data){

        ms_trainval_basins <- ms_site_data %>%
            filter(site_type == 'stream_gauge') %>%
            select(domain, site_code) %>%
            filter(domain %in% addtl_config_list$ms_train_domains,
                   site_code %in% dateranges_ms$basin_id,
                   site_code %in% attr_data$site_code) %>%
            pull(site_code)


        if(! is.na(addtl_config_list$nhm_lineup)){ #either lineup
            ms_trainval_basins <- supplement_with_nhm(ms_trainval_basins, nhm_ms_ids)
            addtl_config_list$ms_train_basins <- supplement_with_nhm(addtl_config_list$ms_train_basins, nhm_ms_ids)
        }

        ms_trainval_basins <- c(ms_trainval_basins, addtl_config_list$ms_train_basins)
        if(! addtl_config_list$camels_only_train){
            train_validate_basins <- c(train_validate_basins, ms_trainval_basins)
        }
        trainval_basins_unmod <- train_validate_basins
    }

    if(addtl_config_list$n_train_val_basins <= length(train_validate_basins)){

        if(! is.na(addtl_config_list$nhm_lineup)){
            warning('multiplying addtl_config_list$n_train_val_basins because NHM')
            addtl_config_list$n_train_val_basins = addtl_config_list$n_train_val_basins * 2
        }

        if(! run_prediction_routine_across_sites){
            train_validate_basins <- sample(train_validate_basins,
                                            addtl_config_list$n_train_val_basins)
        }
    }

    if(addtl_config_list$include_finetune_basins_in_pretrain){

        if(! is.na(addtl_config_list$nhm_lineup)) stop('not set up for nhm')

        replace_inds <- sample(1:length(train_validate_basins),
                               length(addtl_config_list$finetune_basins))
        train_validate_basins[replace_inds] <- addtl_config_list$finetune_basins
        if(any(duplicated(train_validate_basins))) stop('*\n*\n*dupes\n*\n*')
    }

} else if(addtl_config_list$basin_set == 'custom'){

    # stop('out of order (unmaintained)')

    basin_data <- read_delim(file.path(camels_dir,
                                       'camels_attributes_v2.0',
                                       'camels_topo.txt'),
                             delim = ';') %>%
        filter(gauge_id %in% custom_basin_list,
               ! gauge_id %in% addtl_config_list$camels_test_basins)

    area_quants <- quantile(
        basin_data$area_gages2,
        probs = seq(0, 1, length.out = addtl_config_list$n_train_val_basins),
        names = FALSE) %>%
        tibble(join_area = .)

    #nearest join
    data.table::setDT(basin_data)
    data.table::setDT(area_quants)
    train_validate_basins <- basin_data[area_quants,
                                        on = c('area_gages2' = 'join_area'),
                                        roll = 'nearest'] %>%
        as_tibble() %>%
        pull(gauge_id) %>%
        unique()

    if(addtl_config_list$include_ms_data){

        ms_trainval_basins <- ms_site_data %>%
            filter(site_type == 'stream_gauge') %>%
            select(domain, site_code) %>%
            filter(domain %in% addtl_config_list$ms_train_domains,
                   site_code %in% dateranges_ms$basin_id,
                   site_code %in% attr_data$site_code) %>%
            pull(site_code)

        ms_trainval_basins <- c(ms_trainval_basins, addtl_config_list$ms_train_basins)
        train_validate_basins <- c(train_validate_basins, ms_trainval_basins)
    }
    trainval_basins_unmod <- train_validate_basins

    if(addtl_config_list$n_train_val_basins <= length(train_validate_basins)){
        train_validate_basins <- sample(train_validate_basins,
                                        addtl_config_list$n_train_val_basins)
    }

} else { # full

    stop('are ms sites being included properly?')

    #make sure a few big basins end up in there
    n_large_basins <- ceiling(addtl_config_list$n_train_val_basins *
                                  (addtl_config_list$pct_large_basins / 100))
    n_small_basins <- addtl_config_list$n_train_val_basins - n_large_basins

    train_validate_basinsA <- sample(setdiff(omitted_basins, addtl_config_list$camels_test_basins),
                                     n_large_basins)

    if(addtl_config_list$include_ms_data){

        train_validate_basinsB <- sample(c(setdiff(kratzert_basin_list,
                                                   addtl_config_list$camels_test_basins),
                                           setdiff(ms_trainval_basins, unlist(addtl_config_list$ms_test_basins,
                                                                              use.names = FALSE))),
                                         n_small_basins)
    } else {

        train_validate_basinsB <- sample(setdiff(kratzert_basin_list, addtl_config_list$camels_test_basins),
                                         n_small_basins)
    }

    train_validate_basins <- c(train_validate_basinsA, train_validate_basinsB)
    trainval_basins_unmod <- train_validate_basins

    if(addtl_config_list$n_train_val_basins <= length(train_validate_basins)){
        train_validate_basins <- sample(train_validate_basins,
                                        addtl_config_list$n_train_val_basins)
    }
}

all_test_basins <- addtl_config_list$camels_test_basins
if(addtl_config_list$include_ms_data){
    if(new_generalist){
        all_test_basins <- c(all_test_basins, sites_for_routine)
    } else {
        all_test_basins <- c(all_test_basins,
                             unlist(addtl_config_list$ms_test_basins,
                                    use.names = FALSE))
    }
}

if(addtl_config_list$learn_from_test_basins){

    if(! is.na(addtl_config_list$nhm_lineup)) stop('not set up for nhm')
    train_validate_basins[sample(seq_along(train_validate_basins),
                                 length(all_test_basins))] <- all_test_basins
} else {

    still_available_basins <- setdiff(trainval_basins_unmod, all_test_basins)
    replace_these <- train_validate_basins %in% all_test_basins
    train_validate_basins[replace_these] <-
        sample(still_available_basins, sum(replace_these))
}

#make sure no trainval dupes
if(any(duplicated(train_validate_basins))){
    dupe_bool <- duplicated(train_validate_basins)
    still_available_trainbasins <- setdiff(trainval_basins_unmod, train_validate_basins)
    train_validate_basins[dupe_bool] <- sample(
        still_available_trainbasins, sum(dupe_bool))
}

if(any(duplicated(all_test_basins))) stop('*\n*\n*test dupes\n*\n*')

if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 1){

    is_nhm_basin <- substr(train_validate_basins, 1, 4) == 'NHM_'

    write_lines(train_validate_basins[is_nhm_basin], nhm_cfg$train_basin_file)
    write_lines(train_validate_basins[! is_nhm_basin], cfg$train_basin_file)
} else {
    write_lines(train_validate_basins, cfg$train_basin_file)
}



if(! run_prediction_routine_across_sites){
    write_lines(all_test_basins, cfg$test_basin_file)
} else {
    for(ii in seq_len(n_configs)){
        i <- ix * n_configs + ii - n_configs

        write_lines(all_test_basins,
                    gsub('test.txt',
                         file.path(paste0('run', run_seq[i]),
                                   'test.txt'),
                         cfg$test_basin_file))

        if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2){

            all_test_basins0 <- all_test_basins
            extrap_bsns <- grepl('_extrapolate$', all_test_basins0)
            all_test_basins0[extrap_bsns] <- str_match(all_test_basins0[extrap_bsns], '(.*)?_extrapolate')[, 2, drop = TRUE]
            write_lines(paste0('NHM_', all_test_basins0),
                        gsub('test.txt',
                             file.path(paste0('run', run_seq[i]),
                                       'finetune.txt'),
                             cfg$test_basin_file))
        }
    }
}

#NH needs some config elements to be quoted or unquoted, against the will of the yaml package.
#also resolve list -> mapping write issue
for(ii in seq_len(n_configs)){
    i <- ix * n_configs + ii - n_configs

    if(doing_param_search) runid <- run_seq[i]

    cfg_writepath <- file.path(config_dir,
                               glue('{ex}run{cr}/{type}{cr}.yml',
                                    ex = path_extra,
                                    type = ifelse(run_prediction_routine_across_sites,
                                                  'continue', 'run'),
                                    cr = runid))

    cfg_str <- read_file(cfg_writepath) %>% str_replace_all("'", '')
    cfg_str <- gsub('- type', '  type', cfg_str)
    write_file(cfg_str, cfg_writepath)

    if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 1){
        nhm_cfg_writepath <- file.path(config_dir,
                                   glue('{ex}run{cr}/run{cr}nhm.yml',
                                        ex = path_extra,
                                        cr = runid))

        cfg_str <- read_file(nhm_cfg_writepath) %>% str_replace_all("'", '')
        cfg_str <- gsub('- type', '  type', cfg_str)
        write_file(cfg_str, nhm_cfg_writepath)
    }
}


## 3d config fine-tune (REQUIRES INPUT; ignored if finetune_to_test_basins==F and finetune_basins is NA) ####

finetuning <- addtl_config_list$finetune_to_test_basins || ! is.na(addtl_config_list$finetune_basins[1])
# finetuning <- addtl_config_list$finetune_to_test_basins || (! is.na(addtl_config_list$finetune_basins[1]) && ! addtl_config_list$configure_generalist)

if(finetuning){

    if(run_prediction_routine_across_sites || doing_param_search) nfin <- n_configs else nfin <- 1
    for(ii in seq_len(nfin)){
        i <- ix * nfin + ii - nfin

        #create second config file with fine-tuning specs
        cfg2 <- list()

        if(! doing_param_search){
            finetune_basin_file <- glue('{cdr}/run{ri}/run{ri}_finetune.txt',
                                        cdr = config_dir,
                                        ri = runid)
            pbtpf <- glue('{cdr}/run{ri}/finetune_train_ranges.pkl',
                          cdr = config_dir,
                          ri = runid)
        } else {

            runid <- paste0('run', run_seq[i])
            pathextra <- ifelse(run_prediction_routine_across_sites,
                                paste0(runid, '/'),
                                '')

            finetune_basin_file <- glue('{rpd}/{m}{xx}.txt',
                                        m = pathextra,
                                        xx = ifelse(run_prediction_routine_across_sites,
                                                    'continue', 'finetune'),
                                        rpd = runset_parent_dir)
            pbtpf <- glue('{rpd}/{m}finetune_train_ranges.pkl',
                          m = pathextra,
                          rpd = runset_parent_dir)
        }

        # cfg2$experiment_name <- glue('finetune', runid)

        if(length(addtl_config_list$finetune_basins) > 1 || ! is.na(addtl_config_list$finetune_basins)){

            if(run_prediction_routine_across_sites){
                if((! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2) ||
                   addtl_config_list$configure_specialist){
                    finetune_basin_file <- sub('continue', 'finetune', finetune_basin_file)
                } else {
                    finetune_basin_file <- sub('continue', 'test', finetune_basin_file)
                }
            }

            cfg2$train_basin_file <- finetune_basin_file
            cfg2$validation_basin_file <- finetune_basin_file
            if(! addtl_config_list$separate_pretrain_finetune_test){
                cfg2$test_basin_file <- finetune_basin_file
            } else {
                cfg2$test_basin_file <- sub('(finetune|continue)\\.txt', 'test.txt', finetune_basin_file)
            }
            cfg2$per_basin_train_periods_file <- pbtpf
            cfg2$per_basin_validation_periods_file <- gsub('finetune_train_ranges',
                                                           'finetune_validation_ranges',
                                                           cfg2$per_basin_train_periods_file)
            # if(! addtl_config_list$separate_pretrain_finetune_test){
            if(run_prediction_routine_across_sites){
                cfg2$per_basin_test_periods_file <- gsub('finetune_train_ranges',
                                                         'test_ranges',
                                                         cfg2$per_basin_train_periods_file)
            } else {
                cfg2$per_basin_test_periods_file <- gsub('finetune_train_ranges',
                                                         'finetune_test_ranges',
                                                         cfg2$per_basin_train_periods_file)
            }
            # }

            if(! addtl_config_list$include_test_basins_in_finetune){

                to_drop_bsns <- intersect(addtl_config_list$finetune_basins, all_test_basins)

                if(length(to_drop_bsns)){
                    warning(paste('Dropping', paste(to_drop_bsns, collapse = ', '),
                                  '\nfrom finetune because include_test_basins_in_finetune is FALSE'))
                }

                addtl_config_list$finetune_basins <-
                    addtl_config_list$finetune_basins[! addtl_config_list$finetune_basins %in% all_test_basins]
            }

            if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 1){
                addtl_config_list$finetune_basins <- paste0('NHM_', to_drop_bsns)
            }

            if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2){
                addtl_config_list$finetune_basins <- c(addtl_config_list$finetune_basins, nhm_ms_ids)
            }
            addtl_config_list$finetune_basins <- addtl_config_list$finetune_basins[! duplicated(addtl_config_list$finetune_basins)]

            write_lines(addtl_config_list$finetune_basins, ifelse(run_prediction_routine_across_sites,
                                                                  sub('finetune|test', 'continue', finetune_basin_file),
                                                                  finetune_basin_file))

            # if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2){
            #     write_lines(paste0('NHM_', all_test_basins), finetune_basin_file)
            # }

            # if(addtl_config_list$include_finetune_basins_in_pretrain){
            #     stop('not implemented')

            dateranges0 <- dateranges00 <- daterangesnhm2 <- dateranges %>%
                mutate(start_dates = case_when(mindate_q < mindate_f ~ mindate_q, TRUE ~ mindate_f),
                       end_dates = case_when(maxdate_q > maxdate_f ~ maxdate_q, TRUE ~ maxdate_f)) %>%
                select(basin_id, start_dates, end_dates)

                # select(dateranges,
                #                   basin_id,
                #                   start_dates = mindate_q,
                #                   end_dates = maxdate_q)

            train_end_pct <- addtl_config_list$train_split_pct

            if(addtl_config_list$configure_generalist){
                val_end_pct <- 100
            } else if(! addtl_config_list$separate_pretrain_finetune_test || run_prediction_routine_across_sites){
                val_end_pct <- train_end_pct + ((100 - train_end_pct) %/% 2)
            } else {
                val_end_pct <- 100
            }

            sfse0 = if(is.na(specialist_finetune_site_elect)) NULL else specialist_finetune_site_elect
            dateranges0 <- dateranges0 %>%
                filter(basin_id %in% c(addtl_config_list$finetune_basins, sfse0)) %>%
                mutate(train_ends = as.Date(map2_dbl(.x = start_dates,
                                                     .y = end_dates,
                                                     ~quantile(as.numeric(c(.x, .y)),
                                                               probs = train_end_pct / 100,
                                                               names = FALSE)),
                                            origin = as.Date('1970-01-01')),
                       validation_starts = train_ends + 1,
                       validation_ends = as.Date(map2_dbl(.x = start_dates,
                                                          .y = end_dates,
                                                          ~quantile(as.numeric(c(.x, .y)),
                                                                    probs = val_end_pct / 100,
                                                                    names = FALSE)),
                                                 origin = as.Date('1970-01-01')),
                       test_starts = validation_ends + 1)

            # basins00 <- if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2){
            #     paste0('NHM_', all_test_basins)
            # } else {
            #     all_test_basins
            # }

            dateranges00 <- dateranges00 %>%
                filter(basin_id %in% all_test_basins) %>%
                mutate(train_ends = as.Date(map2_dbl(.x = start_dates,
                                                     .y = end_dates,
                                                     ~quantile(as.numeric(c(.x, .y)),
                                                               probs = train_end_pct / 100,
                                                               names = FALSE)),
                                            origin = as.Date('1970-01-01')),
                       validation_starts = train_ends + 1,
                       validation_ends = as.Date(map2_dbl(.x = start_dates,
                                                          .y = end_dates,
                                                          ~quantile(as.numeric(c(.x, .y)),
                                                                    probs = val_end_pct / 100,
                                                                    names = FALSE)),
                                                 origin = as.Date('1970-01-01')),
                       test_starts = validation_ends + 1)

            if(ugly_fix_12000){
                assign(paste0('dateranges00_', addtl_config_list$ms_test_basins$neon), dateranges00, envir = .GlobalEnv)
            }

            if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2){
                daterangesnhm2 <- daterangesnhm2 %>%
                    filter(basin_id %in% paste0('NHM_', sub('_MANUALQ', '', all_test_basins0))) %>%
                    mutate(train_ends = as.Date(map2_dbl(.x = start_dates,
                                                         .y = end_dates,
                                                         ~quantile(as.numeric(c(.x, .y)),
                                                                   probs = train_end_pct / 100,
                                                                   names = FALSE)),
                                                origin = as.Date('1970-01-01')),
                           validation_starts = train_ends + 1,
                           validation_ends = as.Date(map2_dbl(.x = start_dates,
                                                              .y = end_dates,
                                                              ~quantile(as.numeric(c(.x, .y)),
                                                                        probs = 1,
                                                                        names = FALSE)),
                                                     origin = as.Date('1970-01-01')))
            }

            # if(addtl_config_list$configure_specialist){
            #     dateranges_spec <- dateranges_spec %>%
            #         filter(basin_id %in% paste0('NHM_', all_test_basins0)) %>%
            #         mutate(train_ends = as.Date(map2_dbl(.x = start_dates,
            #                                              .y = end_dates,
            #                                              ~quantile(as.numeric(c(.x, .y)),
            #                                                        probs = train_end_pct / 100,
            #                                                        names = FALSE)),
            #                                     origin = as.Date('1970-01-01')),
            #                validation_starts = train_ends + 1,
            #                validation_ends = as.Date(map2_dbl(.x = start_dates,
            #                                                   .y = end_dates,
            #                                                   ~quantile(as.numeric(c(.x, .y)),
            #                                                             probs = 1,
            #                                                             names = FALSE)),
            #                                          origin = as.Date('1970-01-01')))
            # }

            if(! addtl_config_list$separate_pretrain_finetune_test){

                test_ranges <- select(dateranges0,
                                      basin_id,
                                      start_dates = test_starts,
                                      end_dates)

                dateranges0$end_dates = NULL
                dateranges0$test_starts = NULL
            }

            # if(run_prediction_routine_across_sites){
            if((! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2) ||
               addtl_config_list$configure_specialist){

                #may be overridden below
                test_ranges <- dateranges00 %>%
                    select(basin_id,
                           start_dates,
                           end_dates) %>%
                    filter(basin_id %in% all_test_basins)

                train_ranges <- select(dateranges0,
                                       basin_id,
                                       start_dates,
                                       end_dates = train_ends)

                validation_ranges <- select(dateranges0,
                                            basin_id,
                                            start_dates = validation_starts,
                                            end_dates)
            } else if(! addtl_config_list$configure_generalist){

                test_ranges <- dateranges00 %>%
                    select(basin_id,
                           start_dates = test_starts,
                           end_dates) %>%
                    filter(basin_id %in% all_test_basins)

                train_ranges <- select(dateranges0,
                                       basin_id,
                                       start_dates,
                                       end_dates = train_ends)

                if(addtl_config_list$configure_specialist){

                    validation_ranges <- select(dateranges0,
                                                basin_id,
                                                start_dates = validation_starts,
                                                end_dates)
                } else {

                    validation_ranges <- select(dateranges0,
                                                basin_id,
                                                start_dates = validation_starts,
                                                end_dates = validation_ends)
                }
            } else {

                train_ranges <- select(dateranges0,
                                       basin_id,
                                       start_dates,
                                       end_dates = train_ends)

                validation_ranges <- select(dateranges0,
                                            basin_id,
                                            start_dates = validation_starts,
                                            end_dates)

                test_ranges <- dateranges00 %>%
                    select(basin_id,
                           start_dates,
                           end_dates) %>%
                    filter(basin_id %in% all_test_basins)
            }
            # }


            if(addtl_config_list$reverse_finetune_val_test){

                if(addtl_config_list$separate_pretrain_finetune_test){
                    stop('separate_pretrain_finetune_test incompatible with reverse_finetune_val_test')
                }
                xx <- validation_ranges
                validation_ranges <- test_ranges
                test_ranges <- xx
            }

            if(! doing_param_search){
                rangefile_base <- file.path(config_dir,
                                            glue('run', runid))
            } else {
                if(run_prediction_routine_across_sites){
                    rangefile_base <- file.path(runset_parent_dir, runid)
                } else {
                    rangefile_base <- runset_parent_dir
                }
            }

            write_csv(train_ranges,
                      file.path(rangefile_base,
                                glue('{xx}_train_ranges.csv',
                                     xx = ifelse(run_prediction_routine_across_sites,
                                                 'continue', 'finetune'))))
            write_csv(validation_ranges,
                      file.path(rangefile_base,
                                glue('{xx}_validation_ranges.csv',
                                     xx = ifelse(run_prediction_routine_across_sites,
                                                 'continue', 'finetune'))))

            if(addtl_config_list$finetune_predict_on_trainseq){

                full_ranges <- test_ranges
                full_ranges$start_dates <- train_ranges$start_dates

                if(! addtl_config_list$separate_pretrain_finetune_test){
                    write_csv(full_ranges,
                              file.path(rangefile_base,
                                        'finetune_test_ranges.csv'))
                }

            } else if(! addtl_config_list$separate_pretrain_finetune_test ||
                      run_prediction_routine_across_sites){
                write_csv(test_ranges,
                          file.path(rangefile_base,
                                    glue('{xx}test_ranges.csv',
                                         xx = ifelse(run_prediction_routine_across_sites,
                                                '', 'finetune_'))))
            }

        } else {
            cfg2$train_basin_file <- cfg0$test_basin_file
            cfg2$validation_basin_file <- cfg0$test_basin_file
        }

        #IF THE NUMBER OF PARAMS EVER CHANGES, update the filtering of run_params
        #   below and in the training config specification section
        #ALSO: if ever doing a param search WITHOUT finetune, more engineering needed
        #ALSO ALSO: this gets ignored if param searching
        cfg2$learning_rate <- list(`0` = '5e-3',
                                   `5` = '1e-4')
        # cfg2$learning_rate <- list(`0` = '5e-4',
        #                            `5` = '1e-5')
        cfg2$epochs <- 10L
        cfg2$finetune_modules <- list('head', 'lstm')

        # for(i in seq_len(n_configs)){

        cfg2$experiment_name <- paste0('finetune', run_seq[i])

        if(doing_param_search){

            runid <- run_seq[i]
            run_params <- run_param_list[[i]]
            run_params <- run_params[names(run_params) %in%
                                         c('learning_rate2', 'epochs2', 'finetune_modules')]

            for(j in seq_along(run_params)){

                #finalize learning rate steps based on the number of training epochs
                eps <- run_params$epochs2
                names(run_params$learning_rate2[[1]]) <- floor(seq(0, eps, length.out = 3))[1:2]

                p <- run_params[j]
                paramj <- names(p)
                paramj_officialname <- str_match(paramj, '^([^2]+)2?$')[, 2]
                if(paramj %in% c('learning_rate2', 'finetune_modules')){
                    cfg2[[paramj_officialname]] <- p[[1]][[1]]
                } else {
                    cfg2[[paramj_officialname]] <- p[[1]]
                }
            }
        }

        cfg2_writepath <- file.path(config_dir,
                                    glue('{ex}run{cr}/finetune{cr}.yml',
                                         ex = path_extra,
                                         cr = runid))
        write_yaml(cfg2, cfg2_writepath)

        #NH needs some config elements to be quoted or unquoted, against the will of the yaml package
        cfg2_str <- read_file(cfg2_writepath) %>% str_replace_all("'", '')
        write_file(cfg2_str, cfg2_writepath)
        # }
    }

}

if(exists('cfg2') && any(as.numeric(names(cfg2$learning_rate)) >= cfg2$epochs)){
    stop('learning_rate specified for final epoch (or greater)')
}

## 4a write train/validation/test periods CSVs ####

# dateranges <- select(dateranges,
#                      basin_id,
#                      start_dates = mindate_q,
#                      end_dates = maxdate_q)
dateranges <- dateranges %>%
    mutate(start_dates = case_when(mindate_q < mindate_f ~ mindate_q, TRUE ~ mindate_f),
           end_dates = case_when(maxdate_q > maxdate_f ~ maxdate_q, TRUE ~ maxdate_f)) %>%
    select(basin_id, start_dates, end_dates)

if(! is.na(addtl_config_list$nhm_lineup)){

    dateranges_nhm <- dateranges_nhm %>%
        mutate(start_dates = case_when(mindate_q < mindate_f ~ mindate_q, TRUE ~ mindate_f),
               end_dates = case_when(maxdate_q > maxdate_f ~ maxdate_q, TRUE ~ maxdate_f)) %>%
        select(basin_id, start_dates, end_dates)
    # dateranges_nhm <- select(dateranges_nhm,
    #                          basin_id,
    #                          start_dates = mindate_q,
    #                          end_dates = maxdate_q)
}

if(addtl_config_list$learn_from_test_basins){


    if(! is.na(addtl_config_list$nhm_lineup)) stop('not set up for nhm')

    #training is performed on train_split_pct of each basin's dates; the rest are
    #   evenly divided between validation and testing
    train_end_pct <- addtl_config_list$train_split_pct
    val_end_pct <- train_end_pct + ((100 - train_end_pct) %/% 2)
    dateranges <- dateranges %>%
        mutate(train_ends = as.Date(map2_dbl(.x = start_dates,
                                             .y = end_dates,
                                             ~quantile(as.numeric(c(.x, .y)),
                                                       probs = train_end_pct / 100,
                                                       names = FALSE)),
                                    origin = as.Date('1970-01-01')),
               validation_starts = train_ends + 1,
               validation_ends = as.Date(map2_dbl(.x = start_dates,
                                                  .y = end_dates,
                                                  ~quantile(as.numeric(c(.x, .y)),
                                                            probs = val_end_pct / 100,
                                                            names = FALSE)),
                                         origin = as.Date('1970-01-01')),
               test_starts = validation_ends + 1)

    train_ranges <- select(dateranges,
                           basin_id,
                           start_dates,
                           end_dates = train_ends)

    validation_ranges <- select(dateranges,
                                basin_id,
                                start_dates = validation_starts,
                                end_dates = validation_ends)

    test_ranges <- select(dateranges,
                          basin_id,
                          start_dates = test_starts,
                          end_dates)

} else {

    basins00 <- if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2){
        paste0('NHM_', all_test_basins)
    } else {
        all_test_basins
    }

    #testing occurs on the full available temporal range of a subset of basins
    test_ranges <- filter(dateranges,
                          basin_id %in% basins00)

    dateranges <- filter(dateranges,
                         basin_id %in% train_validate_basins)

    #training is performed on train_split_pct of the non-test basins; validation is
    #   performed on the remaining
    dateranges <- dateranges %>%
        mutate(train_ends = as.Date(map2_dbl(.x = start_dates,
                                             .y = end_dates,
                                             ~quantile(as.numeric(c(.x, .y)),
                                                       probs = addtl_config_list$train_split_pct / 100,
                                                       names = FALSE)),
                                    origin = as.Date('1970-01-01')),
               validation_starts = train_ends + 1)

    if(run_prediction_routine_across_sites){

        if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2){

            train_ranges <- select(daterangesnhm2,
                                   basin_id,
                                   start_dates,
                                   end_dates = train_ends)

            validation_ranges <- select(daterangesnhm2,
                                        basin_id,
                                        start_dates = validation_starts,
                                        end_dates = validation_ends)

        # } else if(addtl_config_list$configure_specialist){
        #
        #     train_ranges <- select(dateranges_spec,
        #                            basin_id,
        #                            start_dates,
        #                            end_dates = train_ends)
        #
        #     validation_ranges <- select(dateranges_spec,
        #                                 basin_id,
        #                                 start_dates = validation_starts,
        #                                 end_dates = validation_ends)
        } else {

            train_ranges <- select(dateranges00,
                                   basin_id,
                                   start_dates,
                                   end_dates = train_ends)

            validation_ranges <- select(dateranges00,
                                        basin_id,
                                        start_dates = validation_starts,
                                        end_dates = validation_ends)
        }

        train_ranges_pre <- select(dateranges,
                                   basin_id,
                                   start_dates,
                                   end_dates = train_ends)

        validation_ranges_pre <- select(dateranges,
                                        basin_id,
                                        start_dates = validation_starts,
                                        end_dates)
        write_csv(train_ranges_pre,
                  file.path(runset_parent_dir,
                            'train_ranges.csv'))
        write_csv(validation_ranges_pre,
                  file.path(runset_parent_dir,
                            'validation_ranges.csv'))

    } else {
        train_ranges <- select(dateranges,
                               basin_id,
                               start_dates,
                               end_dates = train_ends)

        validation_ranges <- select(dateranges,
                                    basin_id,
                                    start_dates = validation_starts,
                                    end_dates)
    }

    if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 1){

        dateranges_nhm <- filter(dateranges_nhm,
                                 basin_id %in% train_validate_basins)

        dateranges_nhm <- dateranges_nhm %>%
            mutate(train_ends = as.Date(map2_dbl(.x = start_dates,
                                                 .y = end_dates,
                                                 ~quantile(as.numeric(c(.x, .y)),
                                                           probs = addtl_config_list$train_split_pct / 100,
                                                           names = FALSE)),
                                        origin = as.Date('1970-01-01')),
                   validation_starts = train_ends + 1)

        train_ranges_nhm <- select(dateranges_nhm,
                                   basin_id,
                                   start_dates,
                                   end_dates = train_ends)

        validation_ranges_nhm <- select(dateranges_nhm,
                                        basin_id,
                                        start_dates = validation_starts,
                                        end_dates)
    }

    if((! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2) ||
       addtl_config_list$configure_specialist){

        #may be overridden again below
        dateranges00 %>%
            filter(basin_id == all_test_basins) %>%
            select(basin_id, start_dates, end_dates) %>%
            write_csv(file.path(rangefile_base,
                                'test_ranges.csv'))
    }

}

if(! doing_param_search){
    rangefile_base <- file.path(config_dir,
                                glue('run', runid))
} else {
    rangefile_base <- runset_parent_dir
}

if(run_prediction_routine_across_sites){
    for(ii in seq_len(n_configs)){
        i <- ix * n_configs + ii - n_configs
        runid <- run_seq[i]
        write_csv(train_ranges,
                  file.path(runset_parent_dir,
                            paste0('run', runid),
                            'finetune_train_ranges.csv'))
        write_csv(validation_ranges,
                  file.path(runset_parent_dir,
                            paste0('run', runid),
                            'finetune_validation_ranges.csv'))
    }
} else {
    write_csv(train_ranges,
              file.path(rangefile_base,
                        'train_ranges.csv'))
    write_csv(validation_ranges,
              file.path(rangefile_base,
                        'validation_ranges.csv'))
}

if(addtl_config_list$separate_pretrain_finetune_test && finetuning){
    if(! run_prediction_routine_across_sites){
        write_csv(test_ranges,
                  file.path(rangefile_base,
                            'finetune_test_ranges.csv'))
    }
} else {
    write_csv(test_ranges,
              file.path(rangefile_base,
                        'test_ranges.csv'))
}

if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 1){

    write_csv(test_ranges,
              file.path(rangefile_base,
                        'test_ranges.csv'))

    write_csv(train_ranges_nhm,
              file.path(rangefile_base,
                        'nhm_train_ranges.csv'))
    write_csv(validation_ranges_nhm,
              file.path(rangefile_base,
                        'nhm_validation_ranges.csv'))
}

## 4aA write pretrain_model_loc.txt ####
    if(! is.na(addtl_config_list$pretrained_model)){
        write_lines(addtl_config_list$pretrained_model,
                    file.path(rangefile_base, 'pretrained_model_loc.txt'))
    }

## 4aB force-modify test ranges if necessary ####

    # if(addtl_config_list$guarantee_full_test_daterange){
    #
    # }
## 4aC handle specialist_finetune_site_elect ####

    if(! is.na(specialist_finetune_site_elect)){

        runds = list.dirs(runset_parent_dir, recursive = FALSE)
        runds = str_match(runds, 'run([0-9]+)$')[, 2]
        ixst = n_configs * (ix - 1) + 1
        ixnd = n_configs * (ix)
        runds = runds[ixst:ixnd]

        for(runid in runds){

            sfse = specialist_finetune_site_elect
            rf = file.path(runset_parent_dir, paste0('run', runid))

            #update train
            d1 = read_csv(file.path(rf, 'continue_train_ranges.csv'))
            d2 = read_csv(file.path(rf, 'finetune_train_ranges.csv'))

            if(grepl('^NHM_', sfse)){
                d2[1, ] = filter(dateranges_nhm_ms, basin_id == !!sfse) %>%
                    select(basin_id, start_dates = mindate_q, end_dates = maxdate_q) %>%
                    mutate(train_ends = as.Date(map2_dbl(.x = start_dates,
                                                         .y = end_dates,
                                                         ~quantile(as.numeric(c(.x, .y)),
                                                                   probs = addtl_config_list$train_split_pct / 100,
                                                                   names = FALSE)),
                                                origin = as.Date('1970-01-01')),
                           validation_starts = train_ends + 1) %>%
                    select(basin_id, start_dates, end_dates = train_ends)
            } else {
                d2[1, ] = filter(d1, basin_id == !!sfse)
            }

            d1 = filter(d1, basin_id != !!sfse)
            write_csv(d1, file.path(rf, 'continue_train_ranges.csv'))
            write_csv(d2, file.path(rf, 'finetune_train_ranges.csv'))

            #update val
            d1 = read_csv(file.path(rf, 'continue_validation_ranges.csv'))
            d2 = read_csv(file.path(rf, 'finetune_validation_ranges.csv'))

            if(grepl('^NHM_', sfse)){
                d2[1, ] = filter(dateranges_nhm_ms, basin_id == !!sfse) %>%
                    select(basin_id, start_dates = mindate_q, end_dates = maxdate_q) %>%
                    mutate(train_ends = as.Date(map2_dbl(.x = start_dates,
                                                         .y = end_dates,
                                                         ~quantile(as.numeric(c(.x, .y)),
                                                                   probs = addtl_config_list$train_split_pct / 100,
                                                                   names = FALSE)),
                                                origin = as.Date('1970-01-01')),
                           validation_starts = train_ends + 1) %>%
                    select(basin_id, start_dates = validation_starts, end_dates)
            } else {
                d2[1, ] = filter(d1, basin_id == !!sfse)
            }
            d1 = filter(d1, basin_id != !!sfse)
            write_csv(d1, file.path(rf, 'continue_validation_ranges.csv'))
            write_csv(d2, file.path(rf, 'finetune_validation_ranges.csv'))

            d3 = read_lines(file.path(rf, 'continue.txt'))
            d3 = d3[! d3 == specialist_finetune_site_elect]
            write_lines(d3, file.path(rf, 'continue.txt'))

            file.create(file.path(rf, 'finetune_special.txt'))
            write_lines(specialist_finetune_site_elect, file.path(rf, 'finetune_special.txt'))

            yml = read_yaml(file.path(rf, paste0('finetune', runid, '.yml')))
            yml$train_basin_file = gsub('finetune.txt$', 'finetune_special.txt', yml$train_basin_file)
            yml$validation_basin_file = gsub('finetune.txt$', 'finetune_special.txt', yml$validation_basin_file)
            write_yaml(yml, file.path(rf, paste0('finetune', runid, '.yml')))
            cfgstr <- read_file(file.path(rf, paste0('finetune', runid, '.yml'))) %>% str_replace_all("'", '')
            write_file(cfgstr, file.path(rf, paste0('finetune', runid, '.yml')))
        }
    }


} #close mega outer loop

### 4ZZZ manually create NHM generalist ####

#1. create dir and subdirs
#2. copy useful template files into nhm_goodies
#3. ensure that template train ranges cover full range:
    # xx = read_csv(file.path(goodies, 'continue_train_ranges.csv'))
    # yy = read_csv(file.path(goodies, 'continue_validation_ranges.csv'))
    # xx$end_dates = yy$end_dates
    # write_csv(xx, file.path(goodies, 'continue_train_ranges.csv'))

# goodies = '~/git/macrosheds/qa_experimentation/imputation/src/nh_methods/nhm_goodies'
# runset_parent_dir = "/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs/runs_3813-3842/"
# run_seq = 3813:3842
# regular_generalist = 3003:3032
# rgd_ = '/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs/runs_3003-3032'
# for(i in seq_along(run_seq)){
#
#     r = run_seq[i]
#     cfgd_ <- paste0('run', r)
#     for(f in c('continue_train_ranges.csv', 'continue.txt', 'test.txt', 'test_ranges.csv')){
#         file.copy(file.path(goodies, f), file.path(runset_parent_dir, cfgd_, f))
#     }
#
#     cfg_ = file.path(runset_parent_dir, cfgd_, paste0('continue', r, '.yml'))
#     r0 = regular_generalist[i]
#     file.copy(file.path(rgd_, paste0('run', r0), paste0('continue', r0, '.yml')),  cfg_)
#
#     xx = read_lines(cfg_)
#     old = paste0('runs_', paste(range(regular_generalist), collapse = '-'))
#     new = paste0('runs_', paste(range(run_seq), collapse = '-'))
#     xx = str_replace(xx, old, new)
#     xx = str_replace(xx, paste0('run', r0), cfgd_)
#     xx = str_replace(xx, 'run_configs', 'lstm_configs')
#     xx = str_replace(xx, 'seed: [0-9]+', paste('seed:', sample(99999, 1)))
#     xx = str_replace(xx, 'generalist_specialist_pretrain', 'generalist_specialist_pretrain_NHM_2304_170417')
#     write_lines(xx, cfg_)
# }

### 4aD remove camels target sites from camels pretrain ####
if(length(addtl_config_list$camels_test_basins)){
    rm_bsns = gsub('_GAPPED$', '', addtl_config_list$camels_test_basins)

    trainval_basin_file = file.path(runset_parent_dir, 'train_validate.txt')
    tv_bsns = read_lines(trainval_basin_file)
    tv_bsns = tv_bsns[! tv_bsns %in% rm_bsns]
    write_lines(tv_bsns, trainval_basin_file)

    trainrng_file = file.path(runset_parent_dir, 'train_ranges.csv')
    read_csv(trainrng_file) %>% filter(! basin_id %in% rm_bsns) %>% write_csv(trainrng_file)

    valrng_file = file.path(runset_parent_dir, 'validation_ranges.csv')
    read_csv(valrng_file) %>% filter(! basin_id %in% rm_bsns) %>% write_csv(valrng_file)
}

### 4aE remove files for finetune 2 if this is a generalist ####
if(strategy == 'generalist'){
    rundirs = list.dirs(runset_parent_dir, recursive = FALSE)
    for(rd in rundirs){
        ffs = list.files(rd, pattern = '^finetune', full.names = TRUE)
        file.remove(ffs)
    }
}

### 4aF write strategy courier file for nh_fiddling.py ####

strategy <- ifelse(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2, 'pgdl2', strategy)
write_lines(strategy, '~/git/macrosheds/qa_experimentation/imputation/data/nh_methods/strategy_courier.txt')

### 4aG if features vary across run, split pretrain.yml into several files ####

if(exists('variable_features')){
    for(rd in rundirs){

        runid = str_match(rd, '([0-9]+)$')[, 2]

        file.copy(file.path(rangefile_base, 'pretrain.yml'), rd)
        cont_conf = read_yaml(file.path(rd, paste0('continue', runid, '.yml')))
        pret_conf = read_yaml(file.path(rd, paste0('pretrain.yml')))
        pret_conf$dynamic_inputs = cont_conf$dynamic_inputs
        pret_conf$metrics = list(pret_conf$metrics)
        pret_conf$target_variables = list(pret_conf$target_variables)
        pret_conf$clip_targets_to_zero = list(pret_conf$clip_targets_to_zero)
        try({pret_conf$statics_embedding$hiddens = list(pret_conf$statics_embedding$hiddens)})
        try({pret_conf$dynamics_embedding$hiddens = list(pret_conf$dynamics_embedding$hiddens)})
        write_yaml(pret_conf, file.path(rd, paste0('pretrain.yml')))

        cfgstr <- read_file(file.path(rd, paste0('pretrain.yml'))) %>% str_replace_all("'", '')
        cfgstr <- gsub('- type', '  type', cfgstr)
        cfgstr <- gsub(' yes', ' True', cfgstr)
        cfgstr <- gsub(' no', ' False', cfgstr)

        write_file(cfgstr, file.path(rd, paste0('pretrain.yml')))
    }
    file.remove(file.path(rangefile_base, 'pretrain.yml')) #pretty sure moving this before the next bracket was right
}

### 4aH [OPTIONAL, uncomment] ensure TOMB isn't sneaking into the test set ####

#only set up for generalist atm

# rundirs = list.files(runset_parent_dir, pattern = 'run', full.names = FALSE)
# for(rd in rundirs){
#     testrng_path = file.path(runset_parent_dir, rd, 'test.txt')
#     try({
#         read_lines(testrng_path) %>%
#             grep('TOMB', ., invert = TRUE, value = TRUE) %>%
#             write_lines(testrng_path)
#     }, silent = TRUE)
# }

###!. 4aI ensure random seeds for each run ####

#might want to do this for pretrain.yml too

rundirs = list.files(runset_parent_dir, pattern = 'run', full.names = FALSE)
for(rd in rundirs){
    contfile = list.files(file.path(runset_parent_dir, rd),
                          pattern = 'continue[0-9]+',
                          full.names = TRUE)
    read_lines(contfile) %>%
        c(., paste('seed:', sample(1:99999, size = 1))) %>%
        write_lines(contfile)

    try({
        contfile = list.files(file.path(runset_parent_dir, rd),
                              pattern = 'finetune[0-9]+',
                              full.names = TRUE)
        read_lines(contfile) %>%
            c(., paste('seed:', sample(1:99999, size = 1))) %>%
            write_lines(contfile)
    }, silent = TRUE)
}

### 4aJ [OPTIONAL, uncomment] protect last x years of eval data as a hold-out set ####

# testrng_path = file.path(runset_parent_dir, 'finetune_test_ranges.csv')
# if(file.exists(testrng_path)){
#
#     read_csv(testrng_path) %>%
#         mutate(end_dates = '2019-12-31') %>% #2020 and 2021 can now be holdout
#         write_csv(testrng_path)
# } else {
#     rundirs = list.files(runset_parent_dir, pattern = 'run', full.names = FALSE)
#     for(rd in rundirs){
#         testrng_path = file.path(runset_parent_dir, rd, 'test_ranges.csv')
#         read_csv(testrng_path) %>%
#             mutate(end_dates = '2019-12-31') %>% #2020 and 2021 can now be holdout
#             write_csv(testrng_path)
#     }
# }


### 4aK [NHM] handle NHM option 2 finetune site_code error ####

if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2){

    runds = list.dirs(runset_parent_dir, recursive = FALSE)
    runds = str_match(runds, 'run([0-9]+)$')[, 2]

    for(runid in runds){
        read_lines(file.path(runset_parent_dir, paste0('run', runid), 'finetune.txt')) %>%
            stringr::str_replace('_MANUALQ', '') %>%
            write_lines(file.path(runset_parent_dir, paste0('run', runid), 'finetune.txt'))
    }
}

### 4aL [NHM] ensure that finetune params inherit properly from pretrain for NHM runs ####

if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2){

    runds = list.dirs(runset_parent_dir, recursive = FALSE)
    runds = str_match(runds, 'run([0-9]+)$')[, 2]

    inherited_size <- read_lines(file.path(runset_parent_dir, 'pretrain.yml')) %>%
        str_subset('hidden_size') %>%
        str_match('hidden_size: ([0-9]+)') %>%
        {.[, 2]}

    for(runid in runds){
        cfgx <- read_lines(file.path(runset_parent_dir, paste0('run', runid), paste0('continue', runid, '.yml')))
        cfgx[str_which(cfgx, 'hidden_size')] <- paste('hidden_size:', inherited_size)
        write_lines(cfgx, file.path(runset_parent_dir, paste0('run', runid), paste0('continue', runid, '.yml')))
    }
}

### 4aM [AS-NEEDED; uncomment] add test metrics to configs in run directories ####

# runs_to_modify <- 1468:1520 #generalist
# runs_to_modify <- c(1548:1627, 1748:1937) #specialist
# runs_to_modify <- 2028:2117 #pgdl
#
# rundir_ <- file.path(config_dir, '..', 'runs')
# rundirs_ <- list.files(rundir_)
#
# for(r in runs_to_modify){
#     rund_ <- grep(paste0('^finetune', r, '_'), rundirs_, value = TRUE)
#     if(! length(rund_)){
#         rund_ <- grep(paste0('^run', r, '_'), rundirs_, value = TRUE)
#     }
#     xx <- try(read_lines(file.path(rundir_, rund_, 'config.yml')), silent = TRUE)
#     if(inherits(xx, 'try-error')) next
#     keepinds <- grep('^- KGE|pbias$', xx, invert = TRUE)
#     xx <- xx[keepinds]
#     xxind <- grep('^- NSE$', xx)
#     xx <- c(xx[1:xxind], c('- KGE', '- pbias'), xx[(xxind + 1):length(xx)])
#     write_lines(xx, file.path(rundir_, rund_, 'config.yml'))
# }

###!. 4aN add test metrics to configs in config directories ####

for(r in run_seq){
    cfgd_ <- paste0('run', r)
    xx <- read_lines(file.path(runset_parent_dir, cfgd_, paste0('continue', r, '.yml')))
    xxind <- grep('^- NSE$', xx)
    xx <- c(xx[1:xxind], c('- KGE', '- pbias'), xx[(xxind + 1):length(xx)])
    write_lines(xx, file.path(runset_parent_dir, cfgd_, paste0('continue', r, '.yml')))
}

###! 4aN2 drop any _MANUALQ sites from continue configs ####

for(r in run_seq){

    cfgd_ <- paste0('run', r)

    read_csv(file.path(runset_parent_dir, cfgd_, paste0('continue_train_ranges.csv'))) %>%
        filter(! grepl('_MANUALQ$', basin_id)) %>%
        write_csv(file.path(runset_parent_dir, cfgd_, paste0('continue_train_ranges.csv')))

    try({
        read_csv(file.path(runset_parent_dir, cfgd_, paste0('continue_validation_ranges.csv'))) %>%
            filter(! grepl('_MANUALQ$', basin_id)) %>%
            write_csv(file.path(runset_parent_dir, cfgd_, paste0('continue_validation_ranges.csv')))
    }, silent = TRUE)

    zz = read_lines(file.path(runset_parent_dir, cfgd_, paste0('continue.txt')))
    zz = grep('_MANUALQ$', zz, value = TRUE, invert = TRUE)
    write_lines(zz, file.path(runset_parent_dir, cfgd_, paste0('continue.txt')))
}

### 4aN3 [DISABLED] modify some continue train/val dateranges ####

# for(r in run_seq){
#
#     cfgd_ <- paste0('run', r)
#
#     qqqq = read_csv(file.path(runset_parent_dir, cfgd_, paste0('continue_train_ranges.csv')))
#     qqqq$end_dates[qqqq$basin_id %in% c('MCRA', 'TECR', 'OKSR')] = as.Date('2019-03-01')
#     write_csv(qqqq, file.path(runset_parent_dir, cfgd_, paste0('continue_train_ranges.csv')))
#
#     qqqq = read_csv(file.path(runset_parent_dir, cfgd_, paste0('continue_validation_ranges.csv')))
#     qqqq$start_dates[qqqq$basin_id %in% c('MCRA', 'TECR', 'OKSR')] = as.Date('2019-03-02')
#     write_csv(qqqq, file.path(runset_parent_dir, cfgd_, paste0('continue_validation_ranges.csv')))
#
#     # qqqq = read_csv(file.path(runset_parent_dir, cfgd_, paste0('test_ranges.csv')))
#     # qqqq$start_dates[qqqq$basin_id %in% c('TECR_MANUALQ')] = as.Date('2017-11-26')
#     # qqqq$end_dates[qqqq$basin_id %in% c('TECR_MANUALQ')] = as.Date('2019-12-31')
#     # write_csv(qqqq, file.path(runset_parent_dir, cfgd_, paste0('test_ranges.csv')))
# }

###!. 4aN4 remove NHM_upper_ipswich from pgdl run (doesn't have pet and can't add) ####

if(! is.na(addtl_config_list$nhm_lineup) && addtl_config_list$nhm_lineup == 2){
    for(r in run_seq){

        cfgd_ <- paste0('run', r)

        read_csv(file.path(runset_parent_dir, cfgd_, paste0('continue_train_ranges.csv'))) %>%
            filter(! grepl('NHM_upper_ipswich', basin_id)) %>%
            write_csv(file.path(runset_parent_dir, cfgd_, paste0('continue_train_ranges.csv')))

        read_csv(file.path(runset_parent_dir, cfgd_, paste0('continue_validation_ranges.csv'))) %>%
            filter(! grepl('NHM_upper_ipswich', basin_id)) %>%
            write_csv(file.path(runset_parent_dir, cfgd_, paste0('continue_validation_ranges.csv')))

        zz = read_lines(file.path(runset_parent_dir, cfgd_, paste0('continue.txt')))
        zz = grep('NHM_upper_ipswich', zz, value = TRUE, invert = TRUE)
        write_lines(zz, file.path(runset_parent_dir, cfgd_, paste0('continue.txt')))
    }
}

### 4aO [REQUIRES INPUT, NON-GENERALIST, uncomment] inherit params from best models by site ####

# runs_to_search = c('remainder of process guided runs (no pet)', 'one more process guided attempt? (just finetuning on NHM)') #pgdl param search for paper
# runs_to_search = c("specialist for un-lmable, ungeneralistable sites", "remainder of specialist runs") #specialist param search for paper
#
# for(r in run_seq){
#
#     cfgd_ <- paste0('run', r)
#     target_site <- read_lines(file.path(runset_parent_dir, cfgd_, paste0('test.txt')))
#
#     bestrun = suppressWarnings(read_csv(file.path(config_dir, '..', '..', 'accumulated_results.csv')) %>%
#         filter(run_note %in% runs_to_search,
#                site_code == !!target_site) %>%
#         arrange(desc(fine_NSE)) %>%
#         slice(1) %>%
#         pull(experiment_name))
#
#     runs = list.files(file.path(config_dir, '..', 'runs'))
#     rundir = try(grep(paste0(bestrun, '_'), runs, value = TRUE), silent = TRUE)
#     if(inherits(rundir, 'try-error') || ! length(rundir)){
#         bestrun = sub('run', 'finetune', bestrun)
#         rundir = grep(paste0(bestrun, '_'), runs, value = TRUE)
#     }
#     # best_cfg = read_lines(file.path(config_dir, '..', 'runs', rundir, 'config.yml'))
#
#     best_cfg_ = list.files(config_dir, recursive = TRUE, full.names = TRUE,
#                            pattern = paste0('continue', sub('finetune|run', '', bestrun), '\\.yml'))
#     best_cfg = read_lines(best_cfg_)
#
#     xx <- read_lines(file.path(runset_parent_dir, cfgd_, paste0('continue', r, '.yml')))
#     xx[grep('^output_dropout', xx)] = best_cfg[grep('^output_dropout', best_cfg)]
#     xx[grep('^batch_size', xx)] = best_cfg[grep('^batch_size', best_cfg)]
#     xx[grep('^epochs', xx)] = best_cfg[grep('^epochs', best_cfg)]
#     xx[grep('^  [0-9]+:', xx)] = best_cfg[grep('^  [0-9]+:', best_cfg)]
#
#     best_cfg2_ = list.files(config_dir, recursive = TRUE, full.names = TRUE,
#                             pattern = paste0('finetune', sub('finetune|run', '', bestrun), '\\.yml'))
#     best_cfg2 = read_lines(best_cfg2_)
#
#     xx2 <- read_lines(file.path(runset_parent_dir, cfgd_, paste0('finetune', r, '.yml')))
#     xx2[grep('^  [0-9]+:', xx2)] = best_cfg2[grep('^  [0-9]+:', best_cfg2)]
#     xx2[grep('^epochs', xx2)] = best_cfg2[grep('^epochs', best_cfg2)]
#
#     # print(best_cfg2[grep('head|lstm$', best_cfg2)])
#     zz = sub('finetune_modules: ', '', best_cfg2[grep('head|lstm$', best_cfg2)])
#     zz = sub('- ', '', zz)
#     zz = paste('-', zz)
#     z_ = grep('^- head|lstm', xx2)
#     xx2 = c(xx2[1:(z_[1] - 1)], zz, xx2[(z_[length(z_)] + 1):length(xx2)])
#
#     write_lines(xx, file.path(runset_parent_dir, cfgd_, paste0('continue', r, '.yml')))
#     write_lines(xx2, file.path(runset_parent_dir, cfgd_, paste0('finetune', r, '.yml')))
# }

### 4aP [DISABLED] specialist ensure that finetune train sites are MANUALQ ####

# for(r in run_seq){
#
#     cfgd_ <- paste0('run', r)
#
#     read_csv(file.path(runset_parent_dir, cfgd_, paste0('finetune_train_ranges.csv'))) %>%
#         mutate(basin_id = ifelse(! grepl('_MANUALQ$', basin_id), paste0(basin_id, '_MANUALQ'), basin_id)) %>%
#         write_csv(file.path(runset_parent_dir, cfgd_, paste0('finetune_train_ranges.csv')))
#
#     read_csv(file.path(runset_parent_dir, cfgd_, paste0('finetune_validation_ranges.csv'))) %>%
#         mutate(basin_id = ifelse(! grepl('_MANUALQ$', basin_id), paste0(basin_id, '_MANUALQ'), basin_id)) %>%
#         write_csv(file.path(runset_parent_dir, cfgd_, paste0('finetune_validation_ranges.csv')))
#
#     zz = read_lines(file.path(runset_parent_dir, cfgd_, paste0('finetune_special.txt')))
#     if(! grepl('_MANUALQ$', zz)) zz = paste0(zz, '_MANUALQ')
#     write_lines(zz, file.path(runset_parent_dir, cfgd_, paste0('finetune_special.txt')))
# }

### 4aQ [DISABLED] specialist set train/val/test ranges for e.g. MANQ train MANQ test ####

# for(r in run_seq){
#
#     if(nrow(dateranges00) != 1) stop()
#
#     cfgd_ <- paste0('run', r)
#
#     site <- read_lines(file.path(runset_parent_dir, cfgd_, paste0('test.txt')))
#     dr00 <- get(paste0('dateranges00_', site))
#
#     read_csv(file.path(runset_parent_dir, cfgd_, paste0('finetune_train_ranges.csv'))) %>%
#         mutate(start_dates = dr00$start_dates,
#                end_dates = dr00$train_ends) %>%
#         write_csv(file.path(runset_parent_dir, cfgd_, paste0('finetune_train_ranges.csv')))
#
#     read_csv(file.path(runset_parent_dir, cfgd_, paste0('finetune_validation_ranges.csv'))) %>%
#         mutate(start_dates = dr00$validation_starts,
#                end_dates = dr00$validation_ends) %>%
#         write_csv(file.path(runset_parent_dir, cfgd_, paste0('finetune_validation_ranges.csv')))
#
#     read_csv(file.path(runset_parent_dir, cfgd_, paste0('test_ranges.csv'))) %>%
#         mutate(start_dates = dr00$test_starts,
#                end_dates = dr00$end_dates) %>%
#         write_csv(file.path(runset_parent_dir, cfgd_, paste0('test_ranges.csv')))
# }

### 4aR [AS-NEEDED, SPECIALIST, uncomment] piggyback off generalist ensemble and just run finetune2 ####

# if(length(run_seq) != length(run_range)) stop('oi')
# for(i in seq_along(run_seq)){
#
#     r <- run_seq[i]
#     r_gen <- run_range[i]
#
#     cfgd_ <- paste0('run', r)
#     cfgd_gen <- paste0('run', r_gen)
#
#     target_site <- read_lines(file.path(runset_parent_dir, cfgd_, paste0('test.txt')))
#
#     rundir <- list.files('~/git/macrosheds/qa_experimentation/imputation/src/nh_methods/runs',
#                          pattern = cfgd_gen, full.names = TRUE)
#     if(length(rundir) > 1) stop('oi')
#
#     cfgpath <- file.path(runset_parent_dir, cfgd_, paste0('finetune', r, '.yml'))
#     read_lines(cfgpath) %>%
#         append(paste('base_run_dir:', rundir)) %>%
#         write_lines(cfgpath)
# }

###!. 4aS1 set up generalist param search. clean more junkfiles  ####

runref = read_lines(file.path(runset_parent_dir, 'pretrained_model_loc.txt'))
for(r in run_seq){
    cfgd_ <- paste0('run', r)
    write_lines(paste('base_run_dir:', runref),
                file.path(runset_parent_dir, cfgd_, paste0('continue', r, '.yml')),
                append = TRUE)
}

for(r in run_seq){
    cfgd_ <- paste0('run', r)
    rm_files = list.files(file.path(runset_parent_dir, cfgd_),
                          pattern = 'validation', full.names = TRUE)
    file.remove(rm_files)
}

rm_files = list.files(runset_parent_dir,
                      pattern = 'pretrain|train|valid', full.names = TRUE)
file.remove(rm_files)

### 4aS2 set up specialist param search. clean more junkfiles ####

for(r in run_seq){
    cfgd_ <- paste0('run', r)
    rm_files = list.files(file.path(runset_parent_dir, cfgd_),
                          pattern = 'continue', full.names = TRUE)
    file.remove(rm_files)
}

runref = read_lines('~/git/macrosheds/qa_experimentation/imputation/src/nh_methods/runlists/generalist_search.txt') %>%
    sort()
# grps_tst = c()
for(i in seq_along(run_seq)){

    r = run_seq[i]
    cfgd_ <- paste0('run', r)
    grp = (i - 1) %/% 30 + 1
    # grps_tst = c(grps_tst, grp)
    write_lines(paste('base_run_dir:', file.path('/hpc/home/mjv22/q_sim/lstm_runs', runref[grp])),
                file.path(runset_parent_dir, cfgd_, paste0('finetune', r, '.yml')),
                append = TRUE)
}

rm_files = list.files(runset_parent_dir,
                      pattern = 'pretrain|train|valid', full.names = TRUE)
file.remove(rm_files)

###!. 4aT eliminate validation period! ####

for(r in run_seq){

    cfgd_ <- paste0('run', r)
    phase = if(new_generalist) 'continue' else 'finetune'

    x_train = read_csv(file.path(runset_parent_dir, cfgd_, paste0(phase, '_train_ranges.csv')))
    x_val = read_csv(file.path(runset_parent_dir, cfgd_, paste0(phase, '_validation_ranges.csv')))
    x_train$end_dates = x_val$end_dates
    write_csv(x_train, file.path(runset_parent_dir, cfgd_, paste0(phase, '_train_ranges.csv')))

    read_lines(file.path(runset_parent_dir, cfgd_, paste0(phase, r, '.yml'))) %>%
        str_replace('validate_every: [0-9]+', 'validate_every:') %>%
        write_lines(file.path(runset_parent_dir, cfgd_, paste0(phase, r, '.yml')))
}

for(r in run_seq){
    cfgd_ <- paste0('run', r)
    rm_files = list.files(file.path(runset_parent_dir, cfgd_),
                          pattern = 'validation', full.names = TRUE)
    file.remove(rm_files)
}
rm_files = list.files(file.path(runset_parent_dir),
                      pattern = 'validation', full.names = TRUE)
file.remove(rm_files)

###!. 4aU set paths for DCC ####

newpth = '/home/mike/git/macrosheds/papers/q_sim/in'
oldpth = '/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods'
system(glue("find {runset_parent_dir} -name '*.yml' | xargs sed -e 's|{oldpth}|/hpc/home/mjv22/q_sim|g' -i"))
system(glue("find {runset_parent_dir} -name '*.yml' | xargs sed -e 's|{newpth}|/hpc/home/mjv22/q_sim|g' -i"))
system(glue("find {runset_parent_dir} -name '*.yml' | xargs sed -e 's|run_configs|lstm_configs|g' -i"))
system(glue("find {runset_parent_dir} -name '*.yml' | xargs sed -e 's|num_workers: 48|num_workers: 8|g' -i"))
system(glue("find {runset_parent_dir} -name 'pretrained_model_loc.txt' | xargs sed -e 's|{oldpth}|/hpc/home/mjv22/q_sim|g' -i"))
# find . -name '*.yml' | xargs sed -e 's|/home/mike/git/macrosheds/papers/q_sim/in|/hpc/home/mjv22/q_sim|g' -i

### 4aV [DONE] modify old pretrain configs for DCC ####

## generalist/specialist pretrain (no nhm)

#COPY PRETRAIN FILES FROM 1468 into run_configs/generalist_specialist_pretrain
#COPY PRETRAIN FILES FROM 1628 into run_configs/generalist_specialist_pretrain_NHM

# pretrain_dir = '/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs/generalist_specialist_pretrain'
# setwd(pretrain_dir)
#
# newpth = '/home/mike/git/macrosheds/papers/q_sim/in'
# oldpth = '/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods'
# system(glue("find . -name '*.yml' | xargs sed -e 's|{oldpth}|/hpc/home/mjv22/q_sim|g' -i"))
# system(glue("find . -name '*.yml' | xargs sed -e 's|{newpth}|/hpc/home/mjv22/q_sim|g' -i"))
# system(glue("find . -name '*.yml' | xargs sed -E 's|runs_[0-9]+-[0-9]+|generalist_specialist_pretrain|g' -i"))
#
# x_train = read_csv('train_ranges.csv')
# x_val = read_csv('validation_ranges.csv')
# x_train$end_dates = x_val$end_dates
# write_csv(x_train, 'train_ranges.csv')
#
# read_lines('pretrain.yml') %>%
#     str_replace('validate_every: [0-9]+', 'validate_every:') %>%
#     write_lines('pretrain.yml')

## pgdl generalist/specialist pretrain (nhm)

# 1. copy train and val rangefiles (CSVs) from runs_2633-2662
# 2. copy pretrain.yml from generalist_specialist_pretrain
# 3. run the code below to fix date ranges, repickle, modify config

pretrain_dir = '/home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs/generalist_specialist_pretrain_NHM'
setwd(pretrain_dir)

x = read_csv('train_ranges.csv')
x$end_dates = read_csv('validation_ranges.csv')$end_dates
write_csv(x, 'train_ranges.csv')
file.remove('validation_ranges.csv')

use_condaenv('nh2')
py_run_file(pickle_script)
system('/home/mike/anaconda3/envs/nh2/bin/python /home/mike/git/macrosheds/papers/q_sim/src/lstm_dungeon/dungeon_lvl_2/pickle_pretrain.py /home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs/generalist_specialist_pretrain_NHM')

system(glue("find . -name '*.yml' | xargs sed -E 's|generalist_specialist_pretrain|generalist_specialist_pretrain_NHM|g' -i"))

#also manually modify learning rates if necessary (if any 5eX)
#and num workers set to 8
#and validate_every (should be done)
#and the experiment name (should be good)

###! 4b convert periods CSVs to python dicts ####

# pickle_script <- file.path(config_dir, '../../nh_methods/pickle_pretrain.py')

# pickle_script <- file.path(config_dir, '../../nh_methods/pickle_train_periods.py')
# use_condaenv('nh2')
# py_run_file(pickle_script)

# pickle_script <- '~/git/macrosheds/papers/q_sim/src/lstm_dungeon/dungeon_lvl_2/nh_methods/pickle_all.py'
system('/home/mike/anaconda3/envs/nh2/bin/python /home/mike/git/macrosheds/papers/q_sim/src/lstm_dungeon/dungeon_lvl_2/pickle_all.py /home/mike/git/macrosheds/qa_experimentation/imputation/src/nh_methods/run_configs/runs_3813-3842')

### 5 [OPTIONAL] run LSTM ####

# py_run_file(file.path(config_dir, '../../nh_methods/nh_fiddling.py'))


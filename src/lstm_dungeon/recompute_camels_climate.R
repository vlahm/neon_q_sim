# Spencer Rhea, Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-12

#map Daymet data
folds <- stringr::str_pad(1:18, 2, 'left', '0')
site_map <- list()
for(i in seq_along(folds)){

    site_map[[folds[i]]] <-
        list.files(paste0('in/CAMELS/basin_dataset_public_v1p2/basin_mean_forcing/daymet/',
                          folds[i])) %>%
        str_extract('^([0-9]+)_', group = 1)
}

#load data from  camels, daymet, Aschonitis et al. 2017
camels_sites <- read_delim('in/CAMELS/camels_attributes_v2.0/camels_topo.txt',
                           col_types = cols('c')) %>%
    select(site_code = gauge_id, latitude = gauge_lat,
           longitude = gauge_lon, elev = elev_mean)

camels_shapes <- st_read('in/CAMELS/HCDN_nhru_final_671.shp')

daymet_files <- list.files('in/CAMELS/basin_dataset_public_v1p2/basin_mean_forcing/daymet/',
                           full.names = TRUE, recursive = TRUE)

a_coef <- raster('in/CAMELS/aptt1_30s/aptt1')

dir.create('in/CAMELS/camels_pet_isolate', showWarnings = FALSE)

#recompute climate indices
daymet_with_pet <- tibble()
for(i in seq_along(site_map)){

    print(paste('recomputing CAMELS climate indices; chunk', i, 'of', length(site_map)))

    for(s in site_map[[i]]){

        ## Read in daymet data

        daymet_f <- try({
            read_file(
                glue('in/CAMELS/basin_dataset_public_v1p2/basin_mean_forcing/daymet/{f}/{s}_lump_cida_forcing_leap.txt',
                     f = names(site_map[i]))
            )
        })

        if(inherits(daymet_f, 'try-error')){
            daymet_f <- read_file(grep(s, daymet_files, value = TRUE))
        }

        daymet_f <- str_replace_all(daymet_f, '\t', ' ')

        daymet_d <- try({
            read_delim(daymet_f, delim = ' ', skip = 3) %>%
                mutate(date = ymd(paste(Year, Mnth, Day, sep = '-'))) %>%
                select(date, dayl = `dayl(s)`, prcp = `prcp(mm/day)`, srad = `srad(W/m2)`,
                       swe = `swe(mm)`, tmax = `tmax(C)`, tmin = `tmin(C)`, vp = `vp(Pa)`) %>%
                mutate(site_code = !!s)
        }, silent = TRUE)

        if(inherits(daymet_d, 'try-error')) next

        daymet_d <- left_join(daymet_d, camels_sites, by = 'site_code')

        ## Get alpha

        site_shp <- camels_shapes %>%
            filter(hru_id == !!s) %>%
            sf::st_transform(., sf::st_crs(a_coef))

        if(nrow(site_shp) == 0){

            site_shp <- camels_shapes %>%
                filter(hru_id == !!substr(s, 2, nchar(s))) %>%
                sf::st_transform(., sf::st_crs(a_coef))

            if(nrow(site_shp) == 0){

                site_shp <- camels_shapes %>%
                    filter(hru_id == !!paste0(0, s)) %>%
                    sf::st_transform(., sf::st_crs(a_coef))

                if(nrow(site_shp) == 0) next
            }
        }

        daymet_d$alpha <- terra::crop(a_coef, site_shp) %>%
            terra::mask(site_shp) %>%
            terra::values() %>%
            mean(na.rm = TRUE) %>%
            {. / 100}

        ## compute pet and climate indices

        daymet_d_supp <- calc_daymet_pet(daymet_d) %>%
            mutate(temp = (t_min + t_max)/2) %>%
            mutate(pet = ifelse(pet < 0, 0, pet)) %>%
            mutate(pet = ifelse(is.na(pet), 0, pet))

        write_csv(select(daymet_d_supp, date, site_code, pet),
                  glue('in/CAMELS/camels_pet_isolate/{s}.csv'))

        look <- daymet_d_supp %>%
            filter(date >= ymd('1989-10-01'),
                   date <= ymd('2009-09-30'))

        look <- compute_clim_indices_camels(
            temp = look$temp,
            prec = look$prcp,
            pet = look$pet,
            day = look$date,
            tol = 0.1
        ) %>%
            mutate(site_code = !!s)

        daymet_with_pet <- rbind(daymet_with_pet, look)
    }
}

as_tibble(daymet_with_pet) %>%
    select(site_code, pet_mean, aridity, frac_snow) %>%
    write_csv('in/CAMELS/recomputed_attributes/clim.csv')

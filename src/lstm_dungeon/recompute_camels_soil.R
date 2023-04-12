# Spencer Rhea, Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-12

sm <- suppressMessages
sw <- suppressWarnings

dir.create('in/CAMELS/recomputed_attributes/soil', showWarnings = FALSE)

# boundaries <- st_read('in/NEON/NEONAquaticWatershed/NEON_Aquatic_Watershed.shp')
boundaries <- st_read('in/CAMELS/HCDN_nhru_final_671.shp')
get_nrcs_soils <- function(network, domain, nrcs_var_name, site, ws_boundaries){

    # Use soilDB to download  soil map unit key (mukey) calssification raster
    site_boundary <- ws_boundaries %>%
        filter(hru_id == !!site)

    bb <- sf::st_bbox(site_boundary)

    soil <- try(sw(soilDB::mukey.wcs(aoi = bb, db = 'gssurgo', quiet = TRUE)))

    # should build a chunking method for this
    if(class(soil) == 'try-error'){

        this_var_tib <- tibble(site_code = site,
                               year = NA,
                               var =  names(nrcs_var_name),
                               val = NA)

        return(this_var_tib)
    }

    mukey_values <- unique(soil@data@values)

    #### Grab soil vars

    # Query Soil Data Acess (SDA) to get the component key (cokey) in each mukey.
    #### each map unit made up of componets but components do not have
    #### spatial informaiton associted with them, but they do include information
    #### on the percentage of each mukey that is made up of each component.
    #### This informaiton on compositon is held in the component table in the
    #### comppct_r column, givin in a percent

    mukey_sql <- soilDB::format_SQL_in_statement(mukey_values)
    component_sql <- sprintf("SELECT cokey, mukey, compname, comppct_r, majcompflag FROM component WHERE mukey IN %s", mukey_sql)
    component <- sm(soilDB::SDA_query(component_sql))

    # Check is soil data is available
    if(length(unique(component$compname)) == 1 &&
       unique(component$compname) == 'NOTCOM'){

        this_var_tib <- tibble(site_code = site, year = NA,
                               var =  names(nrcs_var_name), val = NA)

        return(this_var_tib)
    }

    cokey <- unique(component[,1])
    cokey <- data.frame(ID=cokey)

    # Query SDA for the componets to get infromation on their horizons from the
    #### chorizon table. Each componet is made up of soil horizons (layer of soil
    #### vertically) identified by a chkey.
    #### Informaiton about the depth of each horizon is needed to calculate weighted
    #### averges of any parameter for the whole soil column
    #### hzdept_r = depth to top of horizon
    #### hzdepb_r = depth to bottom of horizon
    #### om_r = percent organic matter
    cokey_sql <- soilDB::format_SQL_in_statement(cokey$ID)
    chorizon_sql <- sprintf(paste0('SELECT cokey, chkey, hzname, desgnmaster, hzdept_r, hzdepb_r, ',
                                   paste(nrcs_var_name, collapse = ', '),
                                   ' FROM chorizon WHERE cokey IN %s'), cokey_sql)

    full_soil_data <- sm(soilDB::SDA_query(chorizon_sql))

    # Calculate weighted average for the entire soil column. This involves 3 steps.
    #### First the component's weighted average of all horizones. Second, the
    #### weighted avergae of each component in a map unit. And third, the weighted
    #### average of all mukeys in the watershed (weighted by their area)

    # cokey weighted average of all horizones
    soil_data_joined <- full_join(full_soil_data, component, by = c("cokey"))

    all_soil_vars <- tibble()
    for(s in 1:length(nrcs_var_name)){

        this_var <- unname(nrcs_var_name[s])

        soil_data_one_var <- soil_data_joined %>%
            select(cokey, chkey, hzname, desgnmaster, hzdept_r, hzdepb_r,
                   !!this_var, mukey, compname, comppct_r, majcompflag)

        cokey_size <- soil_data_one_var %>%
            filter(!is.na(.data[[this_var]])) %>%
            mutate(layer_size = hzdepb_r-hzdept_r) %>%
            group_by(cokey) %>%
            summarise(cokey_size = sum(layer_size)) %>%
            ungroup()

        cokey_weighted_av <- soil_data_one_var %>%
            filter(!is.na(.data[[this_var]])) %>%
            mutate(layer_size = hzdepb_r-hzdept_r) %>%
            left_join(., cokey_size, by = 'cokey') %>%
            mutate(layer_prop = layer_size/cokey_size) %>%
            mutate(value_mat_weith = .data[[this_var]] * layer_prop)  %>%
            group_by(mukey, cokey) %>%
            summarise(value_comp = sum(value_mat_weith),
                      comppct_r = unique(comppct_r)) %>%
            ungroup()

        # mukey weighted average of all compenets
        mukey_weighted_av <- cokey_weighted_av %>%
            group_by(mukey) %>%
            mutate(comppct_r_sum = sum(comppct_r)) %>%
            summarise(value_mukey = sum(value_comp*(comppct_r/comppct_r_sum)))

        # Watershed weighted average
        site_boundary_p <- sf::st_transform(site_boundary, crs = sf::st_crs(soil))

        soil_masked <- sw(terra::mask(soil, site_boundary_p))

        watershed_mukey_values <- soil_masked@data@values %>%
            as_tibble() %>%
            filter(!is.na(value)) %>%
            group_by(value) %>%
            summarise(n = n()) %>%
            rename(mukey = value) %>%
            left_join(mukey_weighted_av, by = 'mukey')


        # Info on soil data
        # query SDA's Mapunit Aggregated Attribute table (muaggatt) by mukey.
        # table information found here: https://sdmdataaccess.sc.egov.usda.gov/documents/TableColumnDescriptionsReport.pdf
        # https://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/survey/geo/?cid=nrcs142p2_053631
        # how the tables relate using mukey, cokey, and chkey is here: https://www.nrcs.usda.gov/Internet/FSE_DOCUMENTS/nrcs142p2_050900.pdf
        if(all(is.na(watershed_mukey_values$value_mukey))){

            na_prop <- 100
            watershed_value <- NA
        } else{

            total_cells <- sum(watershed_mukey_values$n)
            na_cells <- watershed_mukey_values %>%
                filter(is.na(value_mukey))
            na_cells <- sum(na_cells$n)

            watershed_mukey_values_weighted <- watershed_mukey_values %>%
                filter(!is.na(value_mukey)) %>%
                mutate(sum = sum(n, na.rm = TRUE)) %>%
                mutate(prop = n/sum)

            watershed_mukey_values_weighted <- watershed_mukey_values_weighted %>%
                mutate(weighted_av = prop*value_mukey)

            watershed_value <- sum(watershed_mukey_values_weighted$weighted_av,
                                   na.rm = TRUE)

            watershed_value <- round(watershed_value, 2)
            na_prop <- round((100*(na_cells/total_cells)), 2)
        }

        this_var_tib <- tibble(year = NA,
                               site_code = site,
                               var =  names(nrcs_var_name[s]),
                               val = watershed_value,
                               pctCellErr = na_prop)

        all_soil_vars <- rbind(all_soil_vars, this_var_tib)
    }

    return(all_soil_vars)

    # #Used to visualize raster
    # for(i in 1:nrow(watershed_mukey_values)){
    #   soil_masked@data@values[soil_masked@data@values == pull(watershed_mukey_values[i,1])] <- pull(watershed_mukey_values[i,3])
    # }
    # soil_masked@data@isfactor <- FALSE
    #
    # mapview::mapview(soil_masked)
    # raster::plot(soil_masked)

    # Other table in SDA system
    # # component table
    # compnent_sql <- soilDB::format_SQL_in_statement(cokey$ID)
    # component_sql <- sprintf("SELECT cokey, runoff, compname, compkind, comppct_r, majcompflag, otherph, localphase, slope_r, hydricrating, taxorder, taxsuborder, taxsubgrp, taxpartsize FROM component WHERE cokey IN %s", compnent_sql)
    # component_return <- soilDB::SDA_query(component_sql)
    #
    # chtext
    # compnent_sql <- soilDB::format_SQL_in_statement(full_soil_data$chkey)
    # component_sql <- sprintf("SELECT texture, stratextsflag, rvindicator, texdesc, chtgkey FROM chtexturegrp WHERE chkey IN %s", compnent_sql)
    # component_return <- soilDB::SDA_query(component_sql)

    # component_sql <- sprintf("SELECT musym, brockdepmin, wtdepannmin, wtdepaprjunmin, niccdcd FROM muaggatt WHERE mukey IN %s", mukey_sql)
    # component_return <- soilDB::SDA_query(component_sql)
    #
    # # corestrictions table
    # corestrictions_sql <- sprintf("SELECT cokey, reskind, resdept_r, resdepb_r, resthk_r FROM corestrictions WHERE cokey IN %s", compnent_sql)
    # corestrictions_return <- soilDB::SDA_query(corestrictions_sql)
    #
    # # cosoilmoist table
    # cosoilmoist_sql <- sprintf("SELECT cokey,  FROM cosoilmoist WHERE cokey IN %s", compnent_sql)
    # corestrictions_return <- soilDB::SDA_query(corestrictions_sql)
    #
    # # pores
    # chkey_sql <- soilDB::format_SQL_in_statement(unique(full_soil_data$chkey))
    # cokey_to_chkey_sql_pores <- sprintf(paste0('SELECT chkey, poresize FROM chpores WHERE chkey IN %s'), chkey_sql)
    #
    # soil_pores <- soilDB::SDA_query(cokey_to_chkey_sql_pores)

}

clst_type <- ifelse(.Platform$OS.type == 'windows', 'PSOCK', 'FORK')
ncores <- parallel::detectCores() %/% 2
clst <- parallel::makeCluster(spec = ncores, type = clst_type)
doParallel::registerDoParallel(clst)

print(paste('recomputing camels soil metrics using', ncores, 'cores. may take a while.'))

sites <- boundaries$hru_id
catch <- foreach(i = seq_along(sites)) %dopar% {

    s <- sites[i]

    soil_tib <- try({
        sm(get_nrcs_soils(
            nrcs_var_name = c(
                'soil_org' = 'om_r', # percent
                'soil_sand' = 'sandtotal_r', # percent
                'soil_silt' = 'silttotal_r', # percent
                'soil_clay' = 'claytotal_r', # percent
                # mass/volume
                # 'soil_partical_density' = 'partdensity',
                # micrometers per second
                # 'soil_ksat' = 'ksat_r', centimeters of water per centimeter of
                # soil, quantity of water that the soil is capable of storing
                # for use by plants
                'soil_awc' = 'awc_r',
                # Water content, 1/10 bar, is the amount of soil water retained
                # at a tension of 15 bars, expressed as a volumetric percentage
                # of the whole soil material.
                'soil_water_0.1bar' = 'wtenthbar_r',
                # Water content, 1/3 bar, is the amount of soil water retained
                # at a tension of 15 bars, expressed as a volumetric percentage
                # of the whole soil material. 15 bar = wilting point
                'soil_water_0.33bar' = 'wthirdbar_r',
                # Water content, 15 bar, is the amount of soil water retained at
                # a tension of 15 bars, expressed as a volumetric percentage of
                # the whole soil material. 15 bar = field capacity
                'soil_water_15bar' = 'wfifteenbar_r',
                # Water content, 0 bar, is the amount of soil water retained at
                # a tension of 15 bars, expressed as a volumetric percentage of
                # the whole soil material.
                'soil_water_0bar' = 'wsatiated_r',
                # percent of carbonates, by weight, in the fraction of the soil
                # less than 2 millimeters in size. 'pf_soil_carbonate' =
                # 'caco3_r', percent, by weight hydrated calcium sulfates in the
                # fraction of the soil less than 20 millimeters in size
                # 'pf_soil_gypsum' = 'gypsum_r', Cation-exchange capacity
                # (CEC-7) is the total amount of extractable cations that can be
                # held by the soil, expressed in terms of milliequivalents per
                # 100 grams of soil at neutrality (pH 7.0)
                'soil_cat_exchange_7' = 'cec7_r',
                # Effective cation-exchange capacity refers to the sum of
                # extractable cations plus aluminum expressed in terms of
                # milliequivalents per 100 grams of soil
                # 'pf_soil_cat_exchange_eff' = 'ecec_r', Electrical conductivity
                # (EC) is the electrolytic conductivity of an extract from
                # saturated soil paste, expressed as decisiemens per meter at 25
                # degrees C. 'pf_soil_elec_cond' = 'ec_r', Sodium adsorption
                # ratio is a measure of the amount of sodium (Na) relative to
                # calcium (Ca) and magnesium (Mg) in the water extract from
                # saturated soil paste. It is the ratio of the Na concentration
                # divided by the square root of one-half of the Ca + Mg
                # concentration. Soils that have SAR values of 13 or more may be
                # characterized by an increased dispersion of organic matter and
                # clay particles, reduced saturated hydraulic conductivity
                # (Ksat) and aeration, and a general degradation of soil
                # structure. 'pf_soil_SAR' = 'sar_r', pH is the 1:1 water
                # method. A crushed soil sample is mixed with an equal amount of
                # water, and a measurement is made of the suspension.
                'soil_ph' = 'ph1to1h2o_r'),
            #Bulk density, one-third bar, is the ovendry weight of the soil
            #material less than 2 millimeters in size per unit volume of soil at
            #water tension of 1/3 bar, expressed in grams per cubic centimeter
            #'pf_soil_bulk_density' = 'dbthirdbar_r'), Linear extensibility
            #refers to the change in length of an unconfined clod as moisture
            #content is decreased from a moist to a dry state. It is an
            #expression of the volume change between the water content of the
            #clod at 1/3- or 1/10-bar tension (33kPa or 10kPa tension) and oven
            #dryness. The volume change is reported as percent change for the
            #whole soil 'soil_linear_extend' = 'lep_r', Liquid limit (LL) is one
            #of the standard Atterberg limits used to indicate the plasticity
            #characteristics of a soil. It is the water content, on a percent by
            #weight basis, of the soil (passing #40 sieve) at which the soil
            #changes from a plastic to a liquid state 'soil_liquid_limit' =
            #'ll_r', Plasticity index (PI) is one of the standard Atterberg
            #limits used to indicate the plasticity characteristics of a soil.
            #It is defined as the numerical difference between the liquid limit
            #and plastic limit of the soil. It is the range of water content in
            #which a soil exhibits the characteristics of a plastic solid.
            #'soil_plasticity_index' = 'pi_r'
            site = s,
            ws_boundaries = boundaries))
    })

    write_csv(soil_tib, glue('in/CAMELS/recomputed_attributes/soil/{s}.csv'))
}

parallel::stopCluster(clst)

all_soil <- list.files('in/CAMELS/recomputed_attributes/soil',
                       full.names = TRUE) %>%
    map_dfr(read_csv) %>%
    filter(! var %in% c('soil_cat_exchange_7', 'soil_ph')) %>%
    mutate(n = nchar(site_code)) %>%
    mutate(site_code = ifelse(n == 7, paste0(0, site_code), site_code)) %>%
    select(-n, -year)

write_csv(all_soil, 'in/CAMELS/recomputed_attributes/soil.csv')

HERE. FIGURE OUT WHICH SOIL.FEATHER IS IN USE AND ONLY WRITE THAT ONE
#is this supposed to be the same file as the line above?
camels_soil <- read_feather('neon_camels_attr/data/large/camels_soil/soil.feather')
camels_soil <- all_soil

camels_soil <- camels_soil %>%
    select(-pctCellErr) %>%
    filter(var %in% c('soil_org', 'soil_sand', 'soil_silt', 'soil_clay')) %>%
    pivot_wider(names_from = 'var', values_from = 'val') %>%
    rename(sand_frac = soil_sand,
           silt_frac = soil_silt,
           clay_frac = soil_clay,
           organic_frac = soil_org)

dir.create('neon_camels_attr/data/ms_attributes/camels/')
write_feather(camels_soil, 'neon_camels_attr/data/ms_attributes/camels/soil.feather')

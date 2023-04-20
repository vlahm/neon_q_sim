# Spencer Rhea, Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-12

neon_daymet <- read_csv('in/NEON/neon_forcings.csv')

neon_sites <- read_csv('in/NEON/neon_site_info.csv') %>%
    filter(SiteID != 'TOOK') %>%
    select(site_code = SiteID, latitude = Latitude, longitude = Longitude)

neon_elevs <- read_csv('in/NEON/neon_site_info2.csv') %>%
    filter(field_site_id %in% neon_sites$site_code) %>%
    select(site_code = field_site_id, elev = field_mean_elevation_m)

neon_boundaries <- st_read('in/NEON/NEONAquaticWatershed/NEON_Aquatic_Watershed.shp') %>%
    filter(SiteID %in% neon_sites$site_code) %>%
    select(site_code = SiteID, geometry)

df <- left_join(neon_daymet, neon_sites, by = 'site_code') %>%
    left_join(., neon_elevs, by = 'site_code')

a_cof <- terra::rast('in/CAMELS/aptt1_30s/aptt1')

neon_alpha <- tibble()
for(i in seq_len(nrow(neon_sites))){

    site_i <- neon_boundaries[i, ] %>%
        sf::st_transform(., sf::st_crs(a_cof)) %>%
        terra::vect(.)

    vals <- terra::crop(a_cof, site_i) %>%
        terra::mask(site_i) %>%
        terra::values() %>%
        {.[, 1]}

    neon_alpha <- tibble(site_code = site_i$site_code,
                       alpha = mean(vals, na.rm = TRUE) / 100) %>%
        bind_rows(neon_alpha)
}

df <- left_join(df, neon_alpha, by = 'site_code')

calc_daymet_pet(df) %>%
    select(date, site_code, dayl, prcp, srad = s_rad, swe, tmax = t_max,
           tmin = t_min, vp, pet) %>%
    # srad is converted to different units in function from camels, converting back
    mutate(srad = srad * 11.57407) %>%
    write_csv('in/NEON/neon_forcings.csv')

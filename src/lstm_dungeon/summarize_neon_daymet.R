# Spencer Rhea, Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-13

library(googledrive)
library(geojsonio)
library(rgee)

boundaries <- st_read('in/NEON/NEONAquaticWatershed/NEON_Aquatic_Watershed.shp') %>%
    filter(SiteID %in% neon_sites) %>%
    select(site_code = SiteID, geometry)

# load authorization email for Google Earth Engine and Google Drive.
# same account must have access to both.
gee_login <- read_lines('cfg/google_auth.cfg', skip_empty_rows = TRUE, n_max = 1) %>%
    str_extract('^([^ #])+')

googledrive::drive_auth(email = gee_login)
rgee::ee_Initialize(user = gee_login, drive = TRUE)
user_info <- rgee::ee_user_info(quiet = TRUE)

#create remote folder to hold watershed boundaries and daymet summaries
asset_folder <- glue('{a}/neon_daymet', a = user_info$asset_home)
rgee::ee_manage_create(asset_folder)

print('uploading NEON watershed boundaries to Google Earth Engine')

#upload boundaries one at a time (too large altogether)
for(i in 1:nrow(boundaries)){
    one_boundary <- boundaries[i, ]
    asset_path <- file.path(asset_folder, one_boundary$site_code)
    sf_as_ee(one_boundary,
             via = 'getInfo_to_asset',
             assetId = asset_path,
             overwrite = TRUE,
             quiet = TRUE)
}

## set up daymet summary job

asset_path <- rgee::ee_manage_assetlist(asset_folder)

ws_boundary_asset <- ee$FeatureCollection(asset_path$ID[1])
for(i in 2:nrow(asset_path)){
    one_ws <- ee$FeatureCollection(asset_path$ID[i])
    ws_boundary_asset <- ws_boundary_asset$merge(one_ws)
}

imgcol <- ee$ImageCollection('NASA/ORNL/DAYMET_V4')$
    filterBounds(ws_boundary_asset)$
    select('dayl','prcp', 'tmin', 'tmax', 'srad', 'swe', 'vp')

results <- ws_boundary_asset$map(function(f){
    fin = imgcol$map(function(i){
        mean = i$reduceRegion(
            geometry = f$geometry(),
            reducer = ee$Reducer$mean(),
            scale = 1000
        )
        f$setMulti(mean)$set(list(date = i$date()))
    })
})$flatten()

gee <- results$select(
    propertySelectors = c('site_code', 'date', 'dayl', 'prcp', 'srad', 'swe',
                          'tmax', 'tmin', 'vp'),
    retainGeometry = FALSE
)

ee_task <- ee$batch$Export$table$toDrive(
    collection = gee,
    description = 'neon_daymet_out',
    fileFormat = 'CSV',
    folder = 'GEE',
    fileNamePrefix = 'rgee'
)

print('summarizing Daymet layers for NEON watersheds. may take hours or even days depending on demand.')

#run job
ee_task$start()
ee_monitoring(ee_task, max_attempts = Inf)

#retrieve results from google drive
temp_rgee <- tempfile(fileext = '.csv')
googledrive::drive_download(file = 'GEE/rgee.csv', temp_rgee)
read_csv(temp_rgee) %>%
    select(date, site_code, dayl, prcp, srad, swe, tmax, tmin, vp) %>%
    write_csv('in/NEON/neon_forcings.csv')

#remove remote file if desired
# googledrive::drive_rm('GEE/rgee.csv')

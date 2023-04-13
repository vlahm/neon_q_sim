# Spencer Rhea, Mike Vlah
# vlahm13@gmail.com
# last edit: 2023-04-12

library(googledrive)
library(geojsonio)
library(rgee)

boundaries <- st_read('in/NEON/NEONAquaticWatershed/NEON_Aquatic_Watershed.shp') %>%
    filter(SiteID %in% neon_sites) %>%
    select(site_code = SiteID, geometry)

conf <- jsonlite::fromJSON('../../data_acquisition/config.json',
                           simplifyDataFrame = FALSE)

#connect rgee to earth engine and python
gee_login <- conf$gee_login_mike

#load authorization file for macrosheds google sheets and drive
#same account must have GEE and GDrive access
# googlesheets4::gs4_auth(email = gee_login)
googledrive::drive_auth(email = gee_login)

#initialize and authorize GEE account
rgee::ee_Initialize(user = gee_login, drive = TRUE)

user_info <- rgee::ee_user_info(quiet = TRUE)

asset_folder <- glue('{a}/neon_daymet', a = user_info$asset_home)
rgee::ee_manage_create(asset_folder)

print('uploading NEON watershed boundaries to Google Earth Engine')

# asset_path <- file.path(asset_folder, 'ws_boundaries')
# rgee::ee_manage_create(asset_path)

# ee_shape <- sf_as_ee(boundaries,
#                      via = 'getInfo_to_asset',
#                      assetId = asset_path,
#                      overwrite = TRUE,
#                      quiet = TRUE)
#
# if('try-error' %in% class(ee_shape)){

for(i in 1:nrow(boundaries)){
    one_boundary <- boundaries[i, ]
    asset_path <- file.path(asset_folder, one_boundary$site_code)
    sf_as_ee(one_boundary,
             via = 'getInfo_to_asset',
             assetId = asset_path,
             overwrite = TRUE,
             quiet = TRUE)
}
# }



asset_path <- rgee::ee_manage_assetlist(asset_folder)

# if(nrow(asset_path) > 1){
#   for(i in 1:nrow(asset_path)){
#
#     if(i == 1){
#       ws_boundary_asset <- ee$FeatureCollection(asset_path$ID[i])
#     }
#     if(i > 1){
#       one_ws <- ee$FeatureCollection(asset_path$ID[i])
#
#       ws_boundary_asset <- ws_boundary_asset$merge(one_ws)
#     }
#   }
# } else{
ws_boundary_asset <- ee$FeatureCollection(asset_path$ID)
# }

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

ee_description <- glue('{n}_{d}_{p}', d = domain, n = network, p = prodname_ms)

ee_task <- ee$batch$Export$table$toDrive(
    collection = gee,
    description = ee_description,
    fileFormat = 'CSV',
    folder = 'GEE',
    fileNamePrefix = 'rgee'
)

ee_task$start()
ee_monitoring(ee_task, quiet = TRUE, max_attempts = Inf)

temp_rgee <- tempfile(fileext = '.csv')

# expo_backoff(
#     expr = {
googledrive::drive_download(file = 'GEE/rgee.csv',
                            temp_rgee,
                            verbose = FALSE)
#     },
#     max_attempts = 5
# ) %>% invisible()

fin_table <- read_csv(temp_rgee)

googledrive::drive_rm('GEE/rgee.csv', verbose = FALSE)

fin_table <- fin_table %>%
    select(date, site_code, dayl, prcp, srad, swe, tmax, tmin, vp)

# if(nrow(fin_table) == 0){
#     return(generate_ms_exception(glue('No data was retrived for {s}',
#                                       s = site_code)))
# }

# dir.create(glue('data/{n}/{d}/ws_traits/daymet/',
#                 n = network,
#                 d = domain),
#            showWarnings = FALSE)

file_path <- glue('data/{n}/{d}/ws_traits/daymet/domain_climate.feather',
                  n = network,
                  d = domain)

write_feather(fin_table, file_path)

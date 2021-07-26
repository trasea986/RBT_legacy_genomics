#Script for extracting environmental variables at sample location

library(tidyverse)
library(terra)

#define projection for this project
projection <- "EPSG: 4326" ## NEED TO DOUBLE CHECK THIS

#set up points and extract
legacy_df <- read.csv('./data/legacy_data_formatted.csv')

#convert to spatial
legacy_spatial <- vect(legacy_df, geom=c("Longitude", "Latitude"), crs=projection)

#create list of env data for ind bioclim files 
bio_files <- list.files(path = './data/wc2.1_30s_bio', pattern = '*.tif', all.files = TRUE, full.names = TRUE)

#load in the bioclim rasters
bio_layers <- rast(bio_files)

#going to crop to save space prior to dealing with correlated variables for bioclim
x_min <- xmin(legacy_spatial)
x_max <- xmax(legacy_spatial)
y_min <- ymin(legacy_spatial)
y_max <- ymax(legacy_spatial)

#create extent object by slightly increasing the min/max values from data
ext <- rast(xmin= (x_min - 0.5), xmax =(x_max + 0.5), 
            ymin= (y_min -0.5), ymax = (y_max + 0.5))

#define crs of this raster
crs(ext) <- crs(bio_layers$wc2.1_30s_bio_1)

#crop all layers in the stack
bio_layers <- crop(bio_layers, ext)

#plot to double check alignment
plot(bio_layers$wc2.1_30s_bio_1)
points(legacy_spatial)


#... and now NorWeST. multiple basins were combined in ArcGis Pro prior to this step
#collaborators not happy with interpolation issues, so ended up doing this manually but keeping the code here for future reference.
# nw_lines <- vect('./data/NorWeST_all_basins.shp')
# 
# #reproject one bioclim layer for cropping and rasterization template
# bio1 <- bio_layers$wc2.1_30s_bio_1
# bio1 <- project(bio1, "+proj=aea +lat_0=30 +lon_0=-114 +lat_1=43 +lat_2=47 +x_0=1500000 +y_0=0 +datum=NAD83 +units=m +no_defs")
# 
# nw_lines_crop <- crop(nw_lines, bio1)
# 
# 


# #extract values
# #method bilinear used to interpolate based on location within the cell and not just cell value
bio_values <- terra::extract(bio_layers, legacy_spatial, xy = TRUE, method = 'bilinear')

# 
# #for norwest, need list of the columns wanted
# nw_list <- list('ELEV',
#                 'CANOPY',
#                 'SLOPE',
#                 'PRECIP',
#                 'CUMDRAINAG',
#                 'Flow_Aug',
#                 'S1_93_11')
# 
# nw_lines_crop_buffer <- buffer(nw_lines_crop, 1000)
# 
# #rasterize the shapefile, reproject, and then extract
# norwest_sites <- lapply(nw_list, function(i){
#   tmp <- rasterize(nw_lines_crop_buffer, bio1, field=paste(i), background = NA)
#   tmp <- project(tmp, "+proj=longlat +datum=WGS84 +no_defs")
#   tmp[tmp<(-20)] <- NA #remove missing data
#   terra::extract(tmp, legacy_spatial, xy = TRUE, method = 'simple', touches = TRUE, na.rm=TRUE)
# })
# 
# #combine all the individual dataframes together
# norwest_df <- do.call("cbind", norwest_sites)
# 
# #remove duplicates (id, x, and y repeat many times from cbind, maybe expore left_join in the future here)
# norwest_df <- norwest_df[!duplicated(as.list(norwest_df))]

norwest_df <- read.csv('./data/legacy_data_nw_man_screened.csv')

#combine this with any other columns of env data
#add back in the population and location information
full_env <- cbind(norwest_df, bio_values)

#x/y from bio_values and long/lat from legacy_df match (which is good) so going to delete the x and y
#full_env <- full_env[!duplicated(as.list(full_env))] only needed if manual/interpolation used

full_env$x <- NULL
full_env$y <- NULL
full_env$ID <- NULL

#export
write.csv(full_env, file = './outputs/full_env.csv', row.names = FALSE)

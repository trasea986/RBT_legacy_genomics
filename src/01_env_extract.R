library(tidyverse)
library(terra)

#define projection for this project
projection <- "EPSG: 4326" ## NEED TO DOUBLE CHECK THIS

#set up points and extract
legacy_df <- read.csv('./data/legacy_data_formatted.csv')

#convert to spatial
legacy_spatial <- vect(legacy_df, geom=c("Longitude", "Latitude"), crs=projection)

#create list of env data
bio_files <- list.files(path = './data/wc2.1_30s_bio', pattern = '*.tif', all.files = TRUE, full.names = TRUE)

#load in the rasters
bio_layers <- rast(bio_files)

#going to crop and export to save space prior to dealing with correlated variables
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

#extract values
bio_values <- raster::extract(bio_layers, legacy_spatial)

#convert to dataframe
bio_values <- as.data.frame(bio_values)

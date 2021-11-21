#will use predicted WorldClim variables
#will need to make a new lfmm2 object that is just uses worldclim
#devtools::install_github("bcm-uga/LEA") need most recent version
library(tidyverse)
library(LEA)
library(terra)
library(gstat)
library(raster)
library(auto)

library(ggmap)
library(maps)
library(mapdata)
library(cowplot)
library(ggthemes)
library(sf)
library(RColorBrewer)
library(maptools)
library(FRK)
library(ggrepel)

#Load the allele frequency matrix with 20 simulated individuals per population
rbeta_data <- readRDS("../outputs/sim_data.rds")
rbeta_data$.id <- NULL #remove .id column that gets generated when reading the data

#pull env data in from rda analysis
env_lfmm <- read.csv('../outputs/env_rda.csv')

#duplicate for each of the 20 individuals
# Create an index of the rows you want with duplications
env_lfmm <- env_lfmm[rep(seq_len(nrow(env_lfmm)), each = 20), ]

pop <- env_lfmm[,1]

#remove columns that do not have the env variables, and remove the variables for anything that is not from worldclim
env_lfmm <- env_lfmm[,-c(1:6)]
env_lfmm <- env_lfmm[,-c(6)]

#older versions of LEA might have needed you to run lfmm2 before, and defined with "object" parameter".
#run lfmm2, generally needs to be run on a cluster
#lfmm2_model <- lfmm2(input = rbeta_data, env = env_lfmm, K = 9)
#saveRDS(lfmm2_model, '../outputs/lfmm2_model_offset.RDS')
#lfmm2_model <- readRDS('../outputs/lfmm2_model_offset.RDS')

#next, need matrix for the two climate change scenarios. This is going to be very similar to the 01 script with bringing in the legacy data and processing to extract the WorldClim rasters
#set up points and extract
#define projection for this project
projection <- "EPSG: 4326"
env_final <- read.csv('../outputs/env_final_pca.csv', stringsAsFactors = FALSE)
legacy_df <- read.csv('../data/legacy_data_formatted.csv')

#remove duplicate rows
legacy_df <- legacy_df[!duplicated(legacy_df), ]

#also remove the duplicated little jack's, keithley, and johnson
legacy_df <- subset(legacy_df, Stream!="Little Jacks Creek -> Jacks Creek - HydroID: 6440" & Stream!="Keithly Creek -> Weiser River - HydroID: 10032" & Stream!="Johnson Creek -> North Fork Boise River - HydroID: 4458")

#remove the hatchery and remove lower dry creek
legacy_df <- legacy_df[-c(7, 14),]

#convert to spatial
legacy_spatial <- vect(legacy_df, geom=c("Longitude", "Latitude"), crs=projection)

#create list of env data for ind bioclim files 
bio_files245 <- list.files(path = '../data/ssp245/2081-2100/', pattern = '*.tif', all.files = TRUE, full.names = TRUE, recursive = TRUE)

#load in the bioclim rasters
bio_layers245 <- rast(bio_files245)

#remove the rasters that we do not need
bio_layers245 <- bio_layers245[[c(15, 2, 3, 6, 7)]]

#extract values
bio_values245 <- terra::extract(bio_layers245, legacy_spatial, xy = TRUE, method = 'bilinear')

bio_values245$x <- NULL
bio_values245$y <- NULL
bio_values245$ID <- NULL

colnames(bio_values245) = c(colnames(env_final[,7:11]))

#add row for each individual
bio_values245 <- bio_values245[rep(seq_len(nrow(bio_values245)), each = 20), ]

#export
write.csv(bio_values245, file = '../outputs/ssp245_env.csv', row.names = FALSE)

#now to do 585/high carbon emission
#create list of env data for ind bioclim files 
bio_files585 <- list.files(path = '../data/ssp585/2081-2100/', pattern = '*.tif', all.files = TRUE, full.names = TRUE, recursive = TRUE)
bio_layers585 <- rast(bio_files585)
bio_layers585 <- bio_layers585[[c(15, 2, 3, 6, 7)]]
bio_values585 <- terra::extract(bio_layers585, legacy_spatial, xy = TRUE, method = 'bilinear')
bio_values585$x <- NULL
bio_values585$y <- NULL
bio_values585$ID <- NULL
colnames(bio_values585) = c(colnames(env_final[,7:11]))
bio_values585 <- bio_values585[rep(seq_len(nrow(bio_values585)), each = 20), ]

write.csv(bio_values585, file = '../outputs/ssp585_env.csv', row.names = FALSE)


# Computing genetic offset must run on the cluster

g_offset245 <- genetic.offset(input = rbeta_data, 
                           env = env_lfmm, new.env = bio_values245, 
                           pop.labels = pop, K = 9)
saveRDS(g_offset245, '../outputs/g_offset245.RDS')

g_offset585 <- genetic.offset(input = rbeta_data, 
                              env = env_lfmm, new.env = bio_values585, 
                              pop.labels = pop, K = 9)
saveRDS(g_offset585, '../outputs/g_offset585.RDS')


g_offset245 <- readRDS('../outputs/g_offset245.RDS')
g_offset585 <- readRDS('../outputs/g_offset585.RDS')

#create raster by interpolating the offset values
#going to crop to save space prior to dealing with correlated variables for bioclim
x_min <- xmin(legacy_spatial)
x_max <- xmax(legacy_spatial)
y_min <- ymin(legacy_spatial)
y_max <- ymax(legacy_spatial)

#create extent object by slightly increasing the min/max values from data
ext <- rast(xmin= (x_min - 0.5), xmax =(x_max + 0.5), 
            ymin= (y_min -0.5), ymax = (y_max + 0.5))

#define crs of this raster
crs(ext) <- crs(bio_layers245$wc2_15)

#set up offset files
offset245 <- as.data.frame(g_offset245)
offset245 <- cbind(offset245, legacy_df)
offset245 <- offset245 %>% select('Latitude', 'Longitude', 'g_offset245')
offset245$SSP <- c('245')
offset245$offset <- offset245$g_offset245
offset245$g_offset245 <- NULL

#set up offset files
offset585 <- as.data.frame(g_offset585)
offset585 <- cbind(offset585, legacy_df)
offset585 <- offset585 %>% select('Latitude', 'Longitude', 'g_offset585')
offset585$SSP <- c('585')
offset585$offset <- offset585$g_offset585
offset585$g_offset585 <- NULL

offset <- rbind(offset245, offset585)

#convert to spatial to change the projection, then go back to dataframe for ggplot
offset_sp <- coordinates


states_plot <- c("idaho")

dmap <- map("state", regions=states_plot, col="transparent", plot=FALSE, fill = TRUE)

area_poly <- map2SpatialPolygons(dmap, IDs=dmap$names, , proj4string=CRS("+proj=longlat +datum=WGS84"))

counties <- map_data("county")

county_sub <- subset(counties, region %in% c("idaho"))

offset_plot <- ggplot(offset, aes(x = Longitude, y = Latitude, color = offset)) +
  geom_point(size = 4, alpha = 0.75) +
  scale_color_gradient("Genetic Offset", low = "blue", high = "red", limits = c(0.75, 1)) +
  geom_polygon(data = county_sub, mapping = aes(x = long, y = lat, group = group), fill = NA, color = "darkgrey") + #darkgrey county lines
  geom_polygon(data = area_poly, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) + #black lines for the states
  facet_grid(cols = vars(SSP)) +
  theme_classic(base_size = 16)

ggsave('../outputs/offset_plot.tiff', plot = offset_plot, height = 8, width = 10, units = "in")

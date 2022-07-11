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
rbeta_data <- readRDS("../../outputs/no_koot/sim_data.rds")
rbeta_data$.id <- NULL #remove .id column that gets generated when reading the data

#pull env data in from rda analysis
env_lfmm <- read.csv('../../outputs/no_koot/env_rda.csv')

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
#saveRDS(lfmm2_model, '../../outputs/no_koot/lfmm2_model_offset.RDS')
#lfmm2_model <- readRDS('../../outputs/no_koot/lfmm2_model_offset.RDS')

#next, need matrix for the two climate change scenarios. This is going to be very similar to the 01 script with bringing in the legacy data and processing to extract the WorldClim rasters
#set up points and extract
#define projection for this project
projection <- "EPSG: 4326"
env_final <- read.csv('../../outputs/no_koot/env_final_pca.csv', stringsAsFactors = FALSE)
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
write.csv(bio_values245, file = '../../outputs/no_koot/ssp245_env.csv', row.names = FALSE)

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

write.csv(bio_values585, file = '../../outputs/no_koot/ssp585_env.csv', row.names = FALSE)


# Computing genetic offset must run on the cluster

g_offset245 <- genetic.offset(input = rbeta_data, 
                           env = env_lfmm, new.env = bio_values245, 
                           pop.labels = pop, K = 9)
saveRDS(g_offset245, '../../outputs/no_koot/g_offset245.RDS')

g_offset585 <- genetic.offset(input = rbeta_data, 
                              env = env_lfmm, new.env = bio_values585, 
                              pop.labels = pop, K = 9)
saveRDS(g_offset585, '../../outputs/no_koot/g_offset585.RDS')


g_offset245 <- readRDS('../../outputs/no_koot/g_offset245.RDS')
g_offset585 <- readRDS('../../outputs/no_koot/g_offset585.RDS')

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
ext <- setValues(ext, 1)

writeRaster(ext, filename = "../../outputs/no_koot/extent.tif")

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

ggsave('../../outputs/no_koot/figures/offset_plot.tiff', plot = offset_plot, height = 8, width = 10, units = "in")


#one question here is what the maps may look like for the amount the environment is changing
bio_files_present <- list.files(path = '../data/wc2.1_30s_bio/', pattern = '*.tif', all.files = TRUE, full.names = TRUE, recursive = TRUE)
bio_layers_present <- rast(bio_files_present)
#note that order here is slightly different because of the different names conventions of present vs. future data
bio_layers_present <- bio_layers_present[[c(7, 12, 13, 16, 17)]]

#rename the layers to match future scenarios
names(bio_layers_present) <- c("wc2_15", "wc2_2",  "wc2_3",  "wc2_6",  "wc2_7")

#for resolution shift, going to crop layers first to save on space
ext <- rast('../../outputs/no_koot/extent.tif')
bio_layers_present <- terra::crop(bio_layers_present, ext)
bio_layers245 <- terra::crop(bio_layers245, ext)
bio_layers585 <- terra::crop(bio_layers585, ext)

#make resolution of present day match:
bio_layers_present <- terra::resample(bio_layers_present, bio_layers585, method = 'bilinear') # resample output

dif_245 <- bio_layers245 - bio_layers_present
dif_585 <- bio_layers585 - bio_layers_present

#now to turn the stack into a data frame, add the ssp to the dataframe, then rbind. facet rows by SSP and columns by environmental variable

dif_245df <- as.data.frame(dif_245, xy = TRUE)
dif_245df$ssp <- c('245')
dif_585df <- as.data.frame(dif_585, xy = TRUE)
dif_585df$ssp <- c('585')

dif_comb <- rbind(dif_585df, dif_245df)

#for plotting, will need to adjust so that the axis on 245 and 585 is the same and easily comparable
# > max(dif_comb$wc2_15)
# [1] 6.655718
# > max(dif_comb$wc2_2)
# [1] 0.2751064
# > max(dif_comb$wc2_3)
# [1] -1.325596
# > max(dif_comb$wc2_6)
# [1] 9.97856
# > max(dif_comb$wc2_7)
# [1] 5.951464
# 
# > min(dif_comb$wc2_15)
# [1] -9.011786
# > min(dif_comb$wc2_2)
# [1] -1.171982
# > min(dif_comb$wc2_3)
# [1] -6.760864
# > min(dif_comb$wc2_6)
# [1] 4.47056
# > min(dif_comb$wc2_7)
# [1] -0.4212608

#going to want to pull in other states, because the cropped area around the points isn't just focused on Idaho

states_plot <- c("idaho", "washington", "oregon", "montana")
dmap <- map("state", regions=states_plot, col="transparent", plot=FALSE, fill = TRUE)
area_poly <- map2SpatialPolygons(dmap, IDs=dmap$names, , proj4string=CRS("+proj=longlat +datum=WGS84"))
counties <- map_data("county")
county_sub <- subset(counties, region %in% c("idaho", "washington", "oregon", "montana"))

cc15 <- ggplot()+
  geom_raster(data = dif_comb, aes(x = x, y = y, fill = wc2_15)) + 
  scale_fill_gradient2(name = "Env Change Value \n (mm/m)", 
                      low = "blue",high = "red", mid = "white",
                      midpoint = 0,
                      guide = "colourbar",
                      limits = c(-9.1, 6.7)) +
  geom_polygon(data = county_sub, mapping = aes(x = long, y = lat, group = group), fill = NA, color = "darkgray", alpha = 0.5) + #darkgrey county lines
  geom_polygon(data = area_poly, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) + #black lines for the states
  geom_point(data = legacy_df, mapping = aes(x = Longitude, y = Latitude)) +
  facet_grid(cols = vars(ssp)) +
  coord_cartesian(xlim = c(-120, -113), ylim = c(42, 49)) +
  ggtitle("Precipitation Seasonality (Coefficient of Variation)") +
  theme_classic(base_size = 16)
ggsave('../../outputs/no_koot/figures/climate_change_plot_15.png', plot = cc15, height = 8, width = 10, units = "in")

cc2 <- ggplot()+
  geom_raster(data = dif_comb, aes(x = x, y = y, fill = wc2_2)) + 
  scale_fill_gradient2(name = "Env Change Value \n (1/10 degree Celsius)", 
                       low = "blue",high = "red", mid = "white",
                       midpoint= 0,
                       guide = "colourbar",
                      limits = c(-1.2, 0.3)) +
  geom_polygon(data = county_sub, mapping = aes(x = long, y = lat, group = group), fill = NA, color = "darkgray", alpha = 0.5) + #darkgrey county lines
  geom_polygon(data = area_poly, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) + #black lines for the states
  geom_point(data = legacy_df, mapping = aes(x = Longitude, y = Latitude)) +
  facet_grid(cols = vars(ssp)) +
  coord_cartesian(xlim = c(-120, -113), ylim = c(42, 49)) +
  ggtitle("Mean Diurnal Range") +
  theme_classic(base_size = 16)
ggsave('../../outputs/no_koot/figures/climate_change_plot_2.png', plot = cc2, height = 8, width = 10, units = "in")

cc3 <- ggplot()+
  geom_raster(data = dif_comb, aes(x = x, y = y, fill = wc2_3)) + 
  scale_fill_gradient2(name = "Env Change Value \n (1/10 degree Celsius)", 
                       low = "blue",high = "red", mid = "white",
                       midpoint = 0,
                       guide = "colourbar",
                      limits = c(-7, 0)) +
  geom_polygon(data = county_sub, mapping = aes(x = long, y = lat, group = group), fill = NA, color = "darkgray", alpha = 0.5) + #darkgrey county lines
  geom_polygon(data = area_poly, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) + #black lines for the states
  geom_point(data = legacy_df, mapping = aes(x = Longitude, y = Latitude)) +
  facet_grid(cols = vars(ssp)) +
  coord_cartesian(xlim = c(-120, -113), ylim = c(42, 49)) +
  ggtitle("Isothermality") +
  theme_classic(base_size = 16)
ggsave('../../outputs/no_koot/figures/climate_change_plot_3.png', plot = cc3, height = 8, width = 10, units = "in")

cc6 <- ggplot()+
  geom_raster(data = dif_comb, aes(x = x, y = y, fill = wc2_6)) + 
  scale_fill_gradient2(name = "Env Change Value \n (1/10 degree Celsius)", 
                       low = "white ",high = "red",
                       guide = "colourbar",
                      limits = c(0, 10)) +
  geom_polygon(data = county_sub, mapping = aes(x = long, y = lat, group = group), fill = NA, color = "darkgray", alpha = 0.5) + #darkgrey county lines
  geom_polygon(data = area_poly, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) + #black lines for the states
  geom_point(data = legacy_df, mapping = aes(x = Longitude, y = Latitude)) +
  facet_grid(cols = vars(ssp)) +
  coord_cartesian(xlim = c(-120, -113), ylim = c(42, 49)) +
  ggtitle("Min Temperature of Coldest Month") +
  theme_classic(base_size = 16)
ggsave('../../outputs/no_koot/figures/climate_change_plot_6.png', plot = cc6, height = 8, width = 10, units = "in")

cc7 <- ggplot()+
  geom_raster(data = dif_comb, aes(x = x, y = y, fill = wc2_7)) + 
  scale_fill_gradient2(name = "Env Change Value  \n (1/10 degree Celsius)", 
                       low = "blue",high = "red", mid = "white",
                       midpoint = 0,
                       guide = "colourbar",
                      limits = c(-0.5, 6)) +
  geom_polygon(data = county_sub, mapping = aes(x = long, y = lat, group = group), fill = NA, color = "darkgray", alpha = 0.5) + #darkgrey county lines
  geom_polygon(data = area_poly, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) + #black lines for the states
  geom_point(data = legacy_df, mapping = aes(x = Longitude, y = Latitude)) +
  facet_grid(cols = vars(ssp)) +
  coord_cartesian(xlim = c(-120, -113), ylim = c(42, 49)) +
  ggtitle("Temperature Annual Range") +
  theme_classic(base_size = 16)
ggsave('../../outputs/no_koot/figures/climate_change_plot_7.png', plot = cc7, height = 8, width = 10, units = "in")

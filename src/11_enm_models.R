#warning: some of the code after the enm is made requires script 10 to be run to pull in the offset values and plot the legacy locations

library(readxl)
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
library(SDMtune)
library(ENMeval)
library(dismo)

#pull in template
#create list of env data for ind bioclim files 
bio_files245 <- list.files(path = '../data/ssp245/2081-2100/', pattern = '*.tif', all.files = TRUE, full.names = TRUE, recursive = TRUE)
bio_layers245 <- rast(bio_files245)
bio_layers245 <- bio_layers245[[c(15, 2, 3, 6, 7)]]

projection <- "EPSG: 4326"

IDFG_Redband_RBT_Summary_PanelInfo_SRN <- read_excel("D:/OneDrive/GEM3_PostDoc/Redband Legacy Genomics/RBT_legacy_genomics/data/IDFG Redband_RBT Summary_PanelInfo-SRN.xlsx")

points <- IDFG_Redband_RBT_Summary_PanelInfo_SRN %>%
  dplyr::select(Longitude, Latitude)

#remove rows with missing data
points <- na.omit(points)

#some problematic points here. probably missing negative sign, but also some way south and in the middle of the ocean
points <- points %>%
  filter(Longitude < -110) %>%
  filter(Latitude > 39)

#quick check of points on the raster... requires script 10 to be donr
#shows these points much greater extent than legacy genetic samples
# plot(bio_layers_present[[1]])
# points(points, col = "blue", cex = 1.5)

#this data set of points is much larger than the previous, so going to repeat some of the steps in 10 for cropping bioclim to match this new extent
#create raster by interpolating the offset values
#going to crop to save space prior to dealing with correlated variables for bioclim
points_spatial <-  vect(points, geom=c("Longitude", "Latitude"), crs=projection)
x_min <- xmin(points_spatial)
x_max <- xmax(points_spatial)
y_min <- ymin(points_spatial)
y_max <- ymax(points_spatial)

#create extent object by slightly increasing the min/max values from data
ext_enm <- rast(xmin= (x_min - 0.5), xmax =(x_max + 0.5), 
            ymin= (y_min -0.5), ymax = (y_max + 0.5))

#define crs of this raster
crs(ext_enm) <- crs(bio_layers245)
ext_enm <- setValues(ext_enm, 1)

writeRaster(ext_enm, filename = "../outputs/extent_enm.tif", overwrite = TRUE)

#now to bring in env data
bio_files_present <- list.files(path = '../data/wc2.1_30s_bio/', pattern = '*.tif', all.files = TRUE, full.names = TRUE, recursive = TRUE)
bio_layers_present <- rast(bio_files_present)
#note that order here is slightly different between present and future because of the different names conventions of present vs. future data
bio_layers_present <- bio_layers_present[[c(7, 12, 13, 16, 17)]]

bio_files245 <- list.files(path = '../data/ssp245/2081-2100/', pattern = '*.tif', all.files = TRUE, full.names = TRUE, recursive = TRUE)
bio_layers245 <- rast(bio_files245)
bio_layers245 <- bio_layers245[[c(15, 2, 3, 6, 7)]]

bio_files585 <- list.files(path = '../data/ssp585/2081-2100/', pattern = '*.tif', all.files = TRUE, full.names = TRUE, recursive = TRUE)
bio_layers585 <- rast(bio_files585)
bio_layers585 <- bio_layers585[[c(15, 2, 3, 6, 7)]]

#rename the layers to match future scenarios. if worried about matching can double check with 'names'
names(bio_layers_present) <- c("wc2_15", "wc2_2",  "wc2_3",  "wc2_6",  "wc2_7")

#for resolution shift, going to crop layers first to save on space
#ext_enm <- rast('../outputs/extent_enm.tif')
bio_layers_present <- terra::crop(bio_layers_present, ext_enm)
bio_layers245 <- terra::crop(bio_layers245, ext_enm)
bio_layers585 <- terra::crop(bio_layers585, ext_enm)

#make resolution of present day match:
bio_layers_present <- terra::resample(bio_layers_present, bio_layers585, method = 'bilinear')

#quick check of points on the raster
plot(bio_layers_present[[2]])
points(points, col = "blue", cex = 1.5)

#give java extra memory
options(java.parameters = "-Xmx6g" )

#Next is to run ENMeval to check regularization multiplier values
#set up model list to test
tune_args_list  <- list(fc = c("L","Q","LQ","LQH", "H"), rm = c(0.5,1:4))

#running enmeval independently to make watching progress easier
enmeval_results <- ENMevaluate(points, bio_layers_present, n.bg=10000, tune.args = tune_args_list, partitions='checkerboard2', algorithm='maxnet')

eval <- eval.results(enmeval_results)
write.csv(eval, "../outputs/enm_eval_results.csv")

#based on the results from enm eval, the feature classes to be used are LQH or just H (AICc < 2), with AUC slightly higher for LQH, so going with that. The regularization multiplier is 0.5

#also need to use raster stacks and not spat raster
bio_layers_present <- stack(bio_layers_present)
bio_layers245 <- stack(bio_layers245)
bio_layers585 <- stack(bio_layers585)

#look and export plots
#little convoluted with how I am doing names here, because the name change was requested in writing steps, and don't want to mess up anything downstream
bio_layers585_plot <- bio_layers585
bio_layers245_plot <- bio_layers245
bio_layers_present_plot <-bio_layers_present

names(bio_layers585_plot) <- c("Precipitation Seasonality", "Mean Diurnal Temp Range",  "Isothermality",  "Min Temp of Coldest Month",  "Temp Annual Range")
names(bio_layers245_plot) <- c("Precipitation Seasonality", "Mean Diurnal Temp Range",  "Isothermality",  "Min Temp of Coldest Month",  "Temp Annual Range")
names(bio_layers_present_plot) <- c("Precipitation Seasonality", "Mean Diurnal Temp Range", "Isothermality", "Min Temp of Coldest Month",  "Temp Annual Range")

plot(bio_layers585_plot)
plot(bio_layers245_plot)
plot(bio_layers_present_plot)

enm_model <- 
  maxent(x=bio_layers_present,
         p=as.data.frame(points),
         removeDuplicates=TRUE, args=c(
           'maximumbackground=10000',
           'defaultprevalence=1.00',
           'betamultiplier=0.5',
           'plots=true',
           'pictures=true',
           'linear=true',
           'quadratic=true',
           'product=false',
           'threshold=false',
           'hinge=true',
           'threads=6',
           'responsecurves=true',
           'jackknife=true',
           'askoverwrite=false',
           'replicates=10',
           'replicatetype=crossvalidate'),
         path = '../outputs/enm')

#determine which model had the best AUC, note that first in sequence is 0
colnames(as.data.frame(enm_model@results))[max.col(as.data.frame(enm_model@results)[c("Test.AUC"),],ties.method="first")]

#species_7 model is best
#extract the model
model_7 <- enm_model@models[[8]]
enm_results <- enm_model@models[[8]]@results

imp <- cbind(row.names(enm_results), enm_results)
cont <- as.data.frame(imp[16:20,])
colnames(cont) <- c("Variable", "Permuatation Importance")
cont$Variable_Readable <- c('Precipitation Seasonality (Coefficient of Variation)', 'Mean Diurnal Range (Mean of monthly (max temp - min temp))', 'Isothermality (BIO2/BIO7) (×100)', 'Min Temperature of Coldest Month', 'Temperature Annual Range (BIO5-BIO6)')
write.csv(cont, '../outputs/enm_contribution_table.csv', row.names = FALSE)


#how to predict distribution across the landscape (not sure if i need this to answer my question)
predict_present <- predict(bio_layers_present, model_7, progress = 'text')

#view map
plot(predict_present)

#now to predict with the future variables
predict_245 <- predict(bio_layers245, model_7, progress = 'text')
predict_585 <- predict(bio_layers585, model_7, progress = 'text')

#next steps: calculate and visualize the amount of change
change_245 <- predict_245 - predict_present
change_585 <- predict_585 - predict_present

#positive values are areas of gain in the future, negative values are loss
#next, make the two plots

present_df <- as.data.frame(predict_present, xy = TRUE)
predict_585_df <- as.data.frame(predict_585, xy = TRUE)
predict_245_df <- as.data.frame(predict_245, xy = TRUE)
change_245_df <- as.data.frame(change_245, xy=TRUE)
change_585_df <- as.data.frame(change_585, xy=TRUE)

#now to make the plots
#bring in original point locations
legacy_df <- read.csv('../data/legacy_data_formatted.csv')
legacy_df <- legacy_df[!duplicated(legacy_df), ]
legacy_df <- subset(legacy_df, Stream!="Little Jacks Creek -> Jacks Creek - HydroID: 6440" & Stream!="Keithly Creek -> Weiser River - HydroID: 10032" & Stream!="Johnson Creek -> North Fork Boise River - HydroID: 4458")
legacy_df <- legacy_df[-c(7, 14),]

#bring in state lines
states_plot <- c("idaho", "washington", "oregon", "montana")
dmap <- map("state", regions=states_plot, col="transparent", plot=FALSE, fill = TRUE)
area_poly <- map2SpatialPolygons(dmap, IDs=dmap$names, , proj4string=CRS("+proj=longlat +datum=WGS84"))
counties <- map_data("county")
county_sub <- subset(counties, region %in% c("idaho", "washington", "oregon", "montana"))

#plot of present
present_enm <- ggplot() + 
  geom_raster(data = present_df, aes(x = x, y = y, fill = layer)) + 
  scale_fill_gradient2("ENM Value",
                       low = 'white', high = 'blue4',
                       na.value = NA) +
  geom_polygon(data = county_sub, mapping = aes(x = long, y = lat, group = group), fill = NA, color = "darkgray", alpha = 0.5) + #darkgrey county lines
  geom_polygon(data = area_poly, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) + #black lines for the states
  geom_point(data = legacy_df, mapping = aes(x = Longitude, y = Latitude), size = 3, color = "yellow") +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Present") +
  coord_cartesian(xlim = c(-120, -113), ylim = c(42, 49)) +
  theme_classic(base_size = 15) +
  theme(legend.position = "right")

ggsave('../outputs/figures/present_enm.png', plot = present_enm, height = 8, width = 8, units = "in")

#predictions
plot_245 <- ggplot() + 
  geom_raster(data = predict_245_df, aes(x = x, y = y, fill = layer)) + 
  scale_fill_gradient2("ENM Value",
                       low = 'white', high = 'blue4',
                       na.value = NA, limits=c(0, 1)) +
  geom_polygon(data = county_sub, mapping = aes(x = long, y = lat, group = group), fill = NA, color = "darkgray", alpha = 0.5) + #darkgrey county lines
  geom_polygon(data = area_poly, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) + #black lines for the states
  geom_point(data = legacy_df, mapping = aes(x = Longitude, y = Latitude), size = 3, color = "yellow") +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("SSP245") +
  coord_cartesian(xlim = c(-120, -113), ylim = c(42, 49)) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none")

plot_585 <- ggplot() + 
  geom_raster(data = predict_585_df, aes(x = x, y = y, fill = layer)) + 
  scale_fill_gradient2("ENM Value",
                       low = 'white', high = 'blue4',
                       na.value = NA, limits=c(0, 1)) +
  geom_polygon(data = county_sub, mapping = aes(x = long, y = lat, group = group), fill = NA, color = "darkgray", alpha = 0.5) + #darkgrey county lines
  geom_polygon(data = area_poly, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) + #black lines for the states
  geom_point(data = legacy_df, mapping = aes(x = Longitude, y = Latitude), size = 3, color = "yellow") +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("SSP585") +
  coord_cartesian(xlim = c(-120, -113), ylim = c(42, 49)) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none")

#get legend from above
legend_enm <- get_legend(ggplot() + 
                           geom_raster(data = predict_585_df, aes(x = x, y = y, fill = layer)) + 
                           scale_fill_gradient2("ENM \nValue",
                                                low = 'white', high = 'blue4',
                                                na.value = NA, limits=c(0, 1)) +
                           geom_polygon(data = county_sub, mapping = aes(x = long, y = lat, group = group), fill = NA, color = "darkgray", alpha = 0.5) + #darkgrey county lines
                           geom_polygon(data = area_poly, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) + #black lines for the states
                           geom_point(data = legacy_df, mapping = aes(x = Longitude, y = Latitude), size = 3, color = "yellow") +
                           xlab("Longitude") +
                           ylab("Latitude") +
                           ggtitle("SSP585") +
                           coord_cartesian(xlim = c(-120, -113), ylim = c(42, 49)) +
                           theme_classic(base_size = 12) +
                           theme(legend.position = "right"))

#change plots
plot_change_245 <- ggplot() + 
  geom_raster(data = change_245_df, aes(x = x, y = y, fill = layer)) + 
  scale_fill_gradient2("ENM Value",
                       low = 'red4', mid = 'white', high = 'blue4',
                       na.value = NA, midpoint = 0, limits = c(-1,1)) +
  geom_polygon(data = county_sub, mapping = aes(x = long, y = lat, group = group), fill = NA, color = "darkgray", alpha = 0.5) + #darkgrey county lines
  geom_polygon(data = area_poly, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) + #black lines for the states
  geom_point(data = legacy_df, mapping = aes(x = Longitude, y = Latitude), size = 2, color = 'yellow') +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("SSP245") +
  coord_cartesian(xlim = c(-120, -113), ylim = c(42, 49)) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none")

plot_change_585 <- ggplot() + 
  geom_raster(data = change_585_df, aes(x = x, y = y, fill = layer)) + 
  scale_fill_gradient2("ENM Value",
                       low = 'red4', mid = 'white', high = 'blue4',
                       na.value = NA, midpoint = 0, limits = c(-1,1)) +
  geom_polygon(data = county_sub, mapping = aes(x = long, y = lat, group = group), fill = NA, color = "darkgray", alpha = 0.5) + #darkgrey county lines
  geom_polygon(data = area_poly, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) + #black lines for the states
  geom_point(data = legacy_df, mapping = aes(x = Longitude, y = Latitude), size = 2, color = 'yellow') +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("SSP585") +
  coord_cartesian(xlim = c(-120, -113), ylim = c(42, 49)) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none")

legend_change <- get_legend(ggplot() + 
                              geom_raster(data = change_585_df, aes(x = x, y = y, fill = layer)) + 
                              scale_fill_gradient2("ENM \nChange",
                                                   low = 'red4', mid = 'white', high = 'blue4',
                                                   na.value = NA, midpoint = 0, limits = c(-1,1)) +
                              geom_polygon(data = county_sub, mapping = aes(x = long, y = lat, group = group), fill = NA, color = "darkgray", alpha = 0.5) + #darkgrey county lines
                              geom_polygon(data = area_poly, mapping = aes(x = long, y = lat, group = group), color = "black", fill = NA) + #black lines for the states
                              geom_point(data = legacy_df, mapping = aes(x = Longitude, y = Latitude), size = 2, color = 'yellow') +
                              xlab("Longitude") +
                              ylab("Latitude") +
                              ggtitle("SSP585") +
                              coord_cartesian(xlim = c(-120, -113), ylim = c(42, 49)) +
                              theme_classic(base_size = 12) +
                              theme(legend.position = "right"))

plot_final <- plot_grid(plot_245, plot_585, legend_enm, plot_change_245, plot_change_585, legend_change, ncol = 3, rel_widths = c(1,1,.4,1,1,.4))


ggsave('../outputs/figures/enm_future.png', plot = plot_final, height = 8, width = 9, units = "in")


#next up is to extract the enm values to points, and then compare to the genetic offset values

g_offset245 <- readRDS('../outputs/g_offset245.RDS')
g_offset585 <- readRDS('../outputs/g_offset585.RDS')
#set up offset files
offset245 <- as.data.frame(g_offset245)
offset245 <- cbind(offset245, legacy_df)
offset245 <- offset245 %>% dplyr::select('Latitude', 'Longitude', 'g_offset245')
offset245$SSP <- c('245')
offset245$offset <- offset245$g_offset245
offset245$g_offset245 <- NULL
offset585 <- as.data.frame(g_offset585)
offset585 <- cbind(offset585, legacy_df)
offset585 <- offset585 %>% dplyr::select('Latitude', 'Longitude', 'g_offset585')
offset585$SSP <- c('585')
offset585$offset <- offset585$g_offset585
offset585$g_offset585 <- NULL

#next is get change values for the legacy_df
legacy_spatial <- vect(legacy_df, geom=c("Longitude", "Latitude"), crs=projection)
#for terra, need rast
change_245 <- rast(change_245)
legacy_245_change <- terra::extract(change_245, legacy_spatial, xy = TRUE, method = 'bilinear')
change_585 <- rast(change_585)
legacy_585_change <- terra::extract(change_585, legacy_spatial, xy = TRUE, method = 'bilinear')

offset245$enm_change <- legacy_245_change$layer
offset585$enm_change <- legacy_585_change$layer

off_enm <- rbind(offset245, offset585)

ENM_vs_Offset <- ggplot(off_enm, aes(x = enm_change, y = offset, color = SSP)) +
  geom_point(size = 3, aes(shape= SSP)) +
  scale_color_manual(values = c('blue', 'darkgreen'))+
  xlab('Change in ENM Value') +
  ylab('Genetic Offset') +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

#now to get a sense for the most imperiled based off of each metric
#keepin in mind higher offset is bad and lower enm is bad
rank_change <- off_enm %>%
  group_by(SSP) %>%
  mutate(rank_offset = rank(-offset)) %>%
  mutate(rank_enm = rank(enm_change))

rank_change <- cbind(rank_change, row.names(off_enm))

colnames(rank_change)[8] <- "Pop"

write.csv(rank_change, "../outputs/enm_offset_change_sites.csv", row.names = FALSE)

ENM_vs_Offset_rank <- ggplot(rank_change, aes(x = rank_offset, y = rank_enm, color = SSP)) +
  geom_point(size = 3, aes(shape= SSP)) +
  scale_color_manual(values = c('blue', 'darkgreen'))+
  xlab('Ranked Change in ENM Value') +
  ylab('Ranked Genetic Offset') +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

legend_compare <- get_legend(ggplot(rank_change, aes(x = rank_offset, y = rank_enm, color = SSP)) +
                              geom_point(size = 3, aes(shape= SSP)) +
                              scale_color_manual(values = c('blue', 'darkgreen'))+
                              xlab('Ranked Change in ENM Value') +
                              ylab('Ranked Genetic Offset') +
                              theme_classic(base_size = 14) +
                              theme(legend.position = "right"))

plot_final_compare <- plot_grid(ENM_vs_Offset, ENM_vs_Offset_rank, legend_compare, ncol = 3, rel_widths = c(1,1,.3))


ggsave('../outputs/figures/enm_offset_compare.png', plot = plot_final_compare, height = 6, width = 10, units = "in")

#other change stats
#average genetic offset

#average ENM change value across SSPs
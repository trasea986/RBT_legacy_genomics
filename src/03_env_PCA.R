#Script to pull in Ecotype data and do a PCA on the environmental variables

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(reshape2)
library(ggrepel)
library(cowplot)

#first bring in the final output from 02_corr_test.R

env_final <- read.csv('../outputs/env_final_pops_temp_swap.csv')


#next is bring in the sample request file to add in ecotypes for color coding the PCA
samples <- read.csv('../data/GEM3 Redband legacy sample request-090419-GIS.csv')
#View(samples)

#request doesn't perfectly match env_final which came from the actual samples that were ran/worked
env_final$Ecotype <- c(
  'Desert', #"Little Jacks Cr-Bruneau drainage"                       
  'Desert', #"Little Jacks Creek -> Jacks Creek - HydroID: 6440" 
  'Desert', #"Big Jacks Creek -> Jacks Creek - HydroID: 6439"
  'Desert', #"Big Jacks Creek -> Jacks Creek - HydroID: 6439"
  'Desert', #"Duncan Creek -> Big Jacks Creek - HydroID: 6563" 
  'Desert', #"Duncan Creek -> Big Jacks Creek - HydroID: 6563" 
  'Desert', #"Williams Creek -> Jordan Creek - HydroID: 11418" 
  'Cool Montane', # "Keithly Creek -> Weiser River - HydroID: 10032"         
  'Cool Montane', # "Keithley Cr-Weiser drainage"                            
  'Cold Montane', #"Johnson Creek -> North Fork Boise River - HydroID: 4458"
  'Cool Montane', #"Little Weiser River -> Weiser River - HydroID: 9840"    
  'Cool Montane', # "Dry Creek -> Boise River - HydroID: 8001"               
  'Cool Montane', #"Dry Creek -> Boise River - HydroID: 8001"               
  'Cold Montane', #"Fawn Cr-NF Payette"                                     
  'Cold Montane', #"Hitt Creek -> Mann Creek - HydroID: 11616"              
  'Cold Montane', #"Fourth Of July Creek -> Mann Creek - HydroID: 11629"    
  'Cool Montane', #"Trail Creek -> Deep Creek - HydroID: 8119"              
  'Cold Montane', #"South Callahan Creek -> Callahan Creek - HydroID: 6553" 
  'Cold Montane', #"South Callahan Creek -> Callahan Creek - HydroID: 6553" 
  'Coastal Hatchery') #"Hayspur Fish Hatchery")

env_final <- env_final %>% select(Stream, Ecotype, everything())

#remove the hatchery and remove lower dry creek
env_final <- env_final[-c(12, 20),]

write.csv(env_final, file = './outputs/env_final_pca.csv', row.names = FALSE)

#now to do the pca
#pulling heavily from this website due to not liking my previous PCA scripts: https://tem11010.github.io/Plotting-PCAs/

#also drop elevation due to issues in later GEA steps
matrix <- env_final[, 4:12]

#replace column names to match fig 6
colnames(matrix) <- c("Canopy","Slope","Drainage","Precipitation","Temp: Mean Diurnal Range","Isothermality","Temp: Min Coldest Month", "Temp: Annual Range", "Temp: Stream")

## center and scale the data
# for (i in 1:length(colnames(matrix))){
#   if (is.numeric(matrix[, i])==TRUE)
#     matrix[, i] <- as.numeric(scale(matrix[, i]))
#   else
#     matrix[, i] <- matrix[, i]
# }

pca1 <- PCA(matrix, graph = FALSE)
#pull out percentages of variance for first four axes
PCA1 <- 48.13
PCA2 <- 32.55
PCA3 <- 6.29
PCA4 <- 5.67


#pull in the first few PCAs
matrix$pc1 <- pca1$ind$coord[, 1] # indexing the first column
matrix$pc2 <- pca1$ind$coord[, 2] 
matrix$pc3 <- pca1$ind$coord[, 3] 
matrix$pc4 <- pca1$ind$coord[, 4] 

#next, extract variable contribution for the axes
pca.vars <- pca1$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")

matrix$Ecotype <- env_final$Ecotype
matrix$Stream <- env_final$Stream

matrix <- separate(matrix, Stream, into = c("Stream", "Hydro"), sep = '-')

#rename the Mann creek to "a" and "b"
matrix$Stream <- c("Little Jacks Cr", "Little Jacks Creek","Big Jacks Creek","Big Jacks Creek","Duncan Creek","Duncan Creek" ,"Williams Creek","Keithly Creek","Keithley Cr","Johnson Creek","Little Weiser River ","Dry Creek","Fawn Cr","Mann Creek (a)","Mann Creek (b)","Trail Creek","South Callahan Creek ","South Callahan Creek")
 

pc12_points <- ggplot(data = matrix, aes(x = pc1, y = pc2, color = Ecotype, shape = Ecotype)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_label_repel(aes(label = Stream),
                   box.padding   = .5, 
                   #point.padding = 5,
                   segment.color = 'darkgrey') +
  geom_point(size = 3) +
  stat_ellipse(geom="polygon", aes(fill = Ecotype),  alpha = 0.05, show.legend = FALSE, level = 0.95) +
  scale_color_manual(values = c('darkblue', 'darkgreen', 'darkred')) +
  scale_fill_manual(values = c('darkblue', 'darkgreen', 'darkred')) +
  scale_shape_manual(values=c(15, 19, 17))+
  xlab(paste("PC1 ", PCA1, "%", sep = '')) + 
  ylab(paste("PC2 ", PCA2, "%", sep = '')) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"),
        legend.position = "bottom")


pc34_points <- ggplot(data = matrix, aes(x = pc3, y = pc4, color = Ecotype, shape = Ecotype)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_label_repel(aes(label = Stream),
                   box.padding   = .5, 
                   #point.padding = 5,
                   segment.color = 'darkgrey') +
  geom_point(size = 3) +
  stat_ellipse(geom="polygon", aes(fill = Ecotype),  alpha = 0.05, show.legend = FALSE, level = 0.95) +
  scale_color_manual(values = c('darkblue', 'darkgreen', 'darkred')) +
  scale_fill_manual(values = c('darkblue', 'darkgreen', 'darkred')) +
  scale_shape_manual(values=c(15, 19, 17))+
  xlab(paste("PC3 ", PCA3, "%", sep = '')) + 
  ylab(paste("PC4 ", PCA4, "%", sep = '')) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"),
        legend.position = "bottom")


#next is to make plots of variable contribution
#one nice option is to include circle around 1
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
circ <- circleFun(c(0,0),2,npoints = 500)

vars12 <- ggplot() +
  geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_segment(data = pca.vars, aes(x = 0, xend = Dim.1, y = 0, yend = Dim.2),arrow = arrow(length = unit(0.025, "npc"), type = "open"),
               lwd = 1) + 
  geom_text(data = pca.vars, aes(x = Dim.1*1.15, y =  Dim.2*1.15,label = c(colnames(matrix[,1:9]))), check_overlap = F, size = 3) +
    xlab(paste("PC1 ", PCA1, "%", sep = '')) + 
    ylab(paste("PC2 ", PCA2, "%", sep = '')) +
    coord_equal() +
    theme_bw(base_size = 13) +
    theme(panel.grid = element_blank(), panel.border = element_rect(fill= "transparent"))

vars34 <- ggplot() +
  geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_segment(data = pca.vars, aes(x = 0, xend = Dim.3, y = 0, yend = Dim.4),arrow = arrow(length = unit(0.025, "npc"), type = "open"),
               lwd = 1) + 
  geom_text(data = pca.vars, aes(x = Dim.3*1.15, y =  Dim.4*1.15,label = c(colnames(matrix[,1:9]))), check_overlap = F, size = 3) +
  xlab(paste("PC3 ", PCA3, "%", sep = '')) + 
  ylab(paste("PC4 ", PCA4, "%", sep = '')) +
  coord_equal() +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(), panel.border = element_rect(fill= "transparent"))


plot_grid(pc12_points, vars12, pc34_points, vars34,
          labels = c("A)", "B)", "C)" ,"D)"), nrow = 2)


ggsave('../outputs/figures/pca_env.png', plot = last_plot(), width = 40, height = 40, units = "cm")

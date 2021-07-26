#Correlation script, created to give env options after 6/10/20 meeting

library(tidyverse)
library(caret)
library(corrplot)
library(psych)

full_env <- read.csv('./outputs/full_env.csv')

#package only wants the matrix of values, so drop first 5 columns (but save stream names, useful for later)
Stream <- full_env$Stream 
full_env <- full_env[,5:ncol(full_env)]

#next renaming the columns to just be the biolcim variable # / deal with other shorthand
#view current
colnames(full_env)
colnames(full_env) <- c('ELEV', 'CANOPY', 'SLOPE', 'PRECIP', 'CUMDRAINAG', 'FLOW_Aug', 'Stream_temp_93_11', 'Annual Mean Temperature', 'Mean Temperature of Warmest Quarter', 'Mean Temperature of Coldest Quarter', 'Annual Precipitation', 'Precipitation of Wettest Month', 'Precipitation of Driest Month', 'Precipitation Seasonality', 'Precipitation of Wettest Quarter', 'Precipitation of Driest Quarter', 'Precipitation of Warmest Quarter', 'Precipitation of Coldest Quarter', 'Mean Diurnal Range', 'Isothermality ', 'Temperature Seasonality', 'Max Temperature of Warmest Month', 'Min Temperature of Coldest Month', 'Temperature Annual Range', 'Mean Temperature of Wettest Quarter', 'Mean Temperature of Driest Quarter')

#main correlation plot, then save it
cor_full <- corr.test(full_env, method = "spearman", adjust = "none")
cor_full$r[cor_full$p > 0.05] <- 0 #this puts blanks in for not significant correlations


png('outputs/init_corr.png',   
    width     = 3.25,
    height    = 3.25,
    units     = "in",
    res       = 600,
    pointsize = 4)
cor_init_plot <- corrplot(cor_full$r, order = "FPC", method = "number", type = "lower", tl.cex = 0.7, tl.col = rgb(0, 0, 0))
dev.off()

#stepwise removal
#find attributes that are highly corrected
highlyCorrelated <- findCorrelation(cor_full$r, cutoff=0.80, exact = TRUE)
#ID names of highly correlated
highlyCorCol <- colnames(full_env)[highlyCorrelated]
#remove highly correlated variables
env_cl <- full_env[, -which(colnames(full_env) %in% highlyCorCol)]

#retest correlation and save matrix + plot
cl_cor <- corr.test(env_cl, method = "spearman")


png('outputs/final_corr.png',   
    width     = 3.25,
    height    = 3.25,
    units     = "in",
    res       = 600,
    pointsize = 4)
corrplot(cl_cor$r, order = "FPC", method = "number", type = "lower", tl.cex = 0.7, tl.col = rgb(0, 0, 0))
dev.off()

write.csv(env_cl, file = './outputs/env_final_matrix.csv', row.names = FALSE)

#if wanting a df with the stream names included with the environmental variables
env_final <- cbind(streams, env_cl)
write.csv(env_final, file = './outputs/env_final_pops.csv', row.names = FALSE)





# request to swap bio X with stream temp instead of air temp for obvious biological reasons
env_final <- read.csv('./outputs/env_final_pops.csv')

env_final$Mean.Temperature.of.Warmest.Quarter <- NULL
env_final <- cbind(env_final, norwest_df$S1_93_11)

env_final_matrix <- env_final
env_final_matrix$streams <- NULL

colnames(env_final_matrix)
colnames(env_final_matrix) <- c('ELEV', 'SLOPE', 'CUMDRAINAG', 'FLOW_Aug',  'Precipitation Seasonality', 'Mean Diurnal Temperature Range', 'Isothermality', 'Temperature Seasonality', 'Stream_temp_93_11')


#retest correlation and save matrix + plot
cl_cor <- corr.test(env_final_matrix, method = "spearman")


png('outputs/final_corr_strem_swap.png',   
    width     = 3.25,
    height    = 3.25,
    units     = "in",
    res       = 600,
    pointsize = 4)
corrplot(cl_cor$r, order = "FPC", method = "number", type = "lower", tl.cex = 0.7, tl.col = rgb(0, 0, 0))
dev.off()

env_final <- cbind(Stream, env_final_matrix)

#and not true final matrix
write.csv(env_final, file = './outputs/env_final_pops_temp_swap.csv', row.names = FALSE)

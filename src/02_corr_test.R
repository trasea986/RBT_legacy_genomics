#Correlation script, created to give env options after 6/10/20 meeting

library(tidyverse)
library(caret)
library(corrplot)
library(psych)

full_env <- read.csv('./outputs/full_env.csv')

#package only wants the matrix of values, so drop first 5 columns (but save stream names, useful for later)
streams <- full_env$Stream 
full_env <- full_env[,6:ncol(full_env)]

#next renaming the columns to just be the biolcim variable # / deal with other shorthand
colnames(full_env) <- c('bio1', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', "ELEV", "CANOPY", "SLOPE", "PRECIP", "CUMDRAINAG", "Flow_Aug", "Stream_temp_93_11")

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
highlyCorrelated <- findCorrelation(cor_full$r, cutoff=0.8, exact = TRUE)
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
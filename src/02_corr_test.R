#Correlation script, created to give env options after 6/10/20 meeting

library(tidyverse)
library(caret)
library(corrplot)
library(psych)

env_full <- read.csv("Steve_Files/env_final.csv")
pop_name <- env_full[,1] #moving out pop names
env_full <- env_full[,-1]

#### main correlation plot ####
cor_full <- corr.test(env_full, method = "spearman", adjust = "none")
cor_full$r[cor_full$p > 0.05] <- 0
corrplot(cor_full$r, order = "FPC", method = "number", type = "lower", tl.cex = 0.7, tl.col = rgb(0, 0, 0))

#stepwise removal:
# find attributes that are highly corrected
# going to first drop R2 to zero in cases of insignificant p values
highlyCorrelated <- findCorrelation(cor_full$r, cutoff=0.8, exact = TRUE)
# Indentifying Variable Names of Highly Correlated Variables
highlyCorCol <- colnames(env_full)[highlyCorrelated]
# Remove highly correlated variables and create a new dataset
dat_cl <- env_full[, -which(colnames(env_full) %in% highlyCorCol)]

#retest correlation
cl_cor <- corr.test(dat_cl, method = "spearman")
#this will make insignificant values 0
#cl_cor$r[cl_cor$p > 0.05] <- 0
corrplot(cl_cor$r, order = "FPC", method = "number", type = "lower", tl.cex = 0.7, tl.col = rgb(0, 0, 0))






library(tidyverse)
library(data.table)
#for rerunning without kootenai populations:

setwd("./src/no_koot/")

#read in env data
env_final <- read.csv('../../outputs/env_final_pca.csv', stringsAsFactors = FALSE)

#remove duplicate rows
env_final <- env_final[!duplicated(env_final), ]

#also remove the duplicated little jack's, keithley, and johnson, plus koot
env_final <- subset(env_final, Stream!="Little Jacks Creek -> Jacks Creek - HydroID: 6440" & Stream!="Keithly Creek -> Weiser River - HydroID: 10032" & Stream!="Johnson Creek -> North Fork Boise River - HydroID: 4458" & Stream!= "Trail Creek -> Deep Creek - HydroID: 8119" & Stream!="South Callahan Creek -> Callahan Creek - HydroID: 6553")

#next pull out mann and average if using only one Mann allele freq instead of two. saving this old code for later.
# mann <- subset(env_final, Stream == "Hitt Creek -> Mann Creek - HydroID: 11616" | Stream == "Fourth Of July Creek -> Mann Creek - HydroID: 11629")
# mann <- colMeans(mann[sapply(mann, is.numeric)])
# mann <- c('Mann Creek Average', 'Cold Montane', mann)
# env_final <- rbind(env_final, mann)
#remove old mann creek values
# env_final <- subset(env_final, Stream!="Hitt Creek -> Mann Creek - HydroID: 11616" & Stream != "Fourth Of July Creek -> Mann Creek - HydroID: 11629")
# next step is to reorder to match the environmental data
#env_final <- rbind(env_final[1,], env_final[2,], env_final[3,], env_final[4,], env_final[5,], env_final[6,], env_final[7,], env_final[8,], env_final[11,], env_final[9,], env_final[10,])

#read in the data
raw <- read.table("../../data/Legacy02_analyze02.fz", header = FALSE, sep = " ")

#remove blank col at the end. this needs to be changed to be number of pops + 1 after account for col 1-3.
raw <- raw[,-14]

#transpose
raw <- t(raw)

#create names from rows 1 to 3 anc reate column names
my_names <- paste(raw[1,], raw[2,], raw[3,], sep = '_')
my_names <- gsub(" ", "", my_names)
colnames(raw) <- my_names

#delete rows 1 to 3
raw <- as.data.frame(raw[-c(1, 2, 3),])

df <- cbind(pop = c('Little Jacks Cr',
            'Big Jacks Cr',
            'Duncan Cr',
            'Williams Creek',
            'Keithley Cr',
            'Little Weser Cr',
            'Dry Cr',
            'Fawn Cr',
            'Mann Cr1',
            'Mann Cr2'), raw)


saveRDS(df, file = "../../outputs/no_koot/allele_freq.RDS") 
saveRDS(env_final, file = "../../outputs/no_koot/env_rda.RDS") 

fwrite(df, '../../outputs/no_koot/allele_freq.csv', row.names = FALSE)
fwrite(env_final, '../../outputs/no_koot/env_rda.csv', row.names = FALSE)

library(tidyverse)
library(data.table)

#read in env data
env_final <- read.csv('../outputs/env_final_pca.csv', stringsAsFactors = FALSE)

#remove duplicate rows
env_final <- env_final[!duplicated(env_final), ]

#also remove the duplicated little jack's, keithley, and johnson
env_final <- subset(env_final, Stream!="Little Jacks Creek -> Jacks Creek - HydroID: 6440" & Stream!="Keithly Creek -> Weiser River - HydroID: 10032" & Stream!="Johnson Creek -> North Fork Boise River - HydroID: 4458")

#next pull out mann and average if using only one Mann allele freq
# mann <- subset(env_final, Stream == "Hitt Creek -> Mann Creek - HydroID: 11616" | Stream == "Fourth Of July Creek -> Mann Creek - HydroID: 11629")

# mann <- colMeans(mann[sapply(mann, is.numeric)])
# mann <- c('Mann Creek Average', 'Cold Montane', mann)
# env_final <- rbind(env_final, mann)

#remove old mann creek values
# env_final <- subset(env_final, Stream!="Hitt Creek -> Mann Creek - HydroID: 11616" & Stream != "Fourth Of July Creek -> Mann Creek - HydroID: 11629")

# next step is to reorder to match the environmental data
#env_final <- rbind(env_final[1,], env_final[2,], env_final[3,], env_final[4,], env_final[5,], env_final[6,], env_final[7,], env_final[8,], env_final[11,], env_final[9,], env_final[10,])


#read in the data
raw <- read.table("../data/Legacy02_analyze01_Mann2.fz", header = FALSE, sep = " ")
raw <- raw[,-15]

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
            'Mann Cr2',
            'Trail Cr',
            'S.F. Callahan Cr'), raw)


saveRDS(df, file = "../outputs/allele_freq.RDS") 
saveRDS(env_final, file = "../outputs/env_rda.RDS") 

fwrite(df, '../outputs/allele_freq.csv', row.names = FALSE)
fwrite(env_final, '../outputs/env_rda.csv', row.names = FALSE)

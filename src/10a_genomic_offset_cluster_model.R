print('Script starting')
time <- Sys.time()
print(time)

library(LEA)

print('Packages loaded')

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

bio_values245 <- read.csv('../outputs/ssp245_env.csv')
bio_values585 <- read.csv('../outputs/ssp585_env.csv')

print('Data loaded')
time <- Sys.time()
print(time)

g_offset245 <- genetic.offset(input = rbeta_data, 
                              env = env_lfmm, new.env = bio_values245, 
                              pop.labels = pop, K = 9)
print('Genetic offset ssp245 calculated')
time <- Sys.time()
print(time)

saveRDS(g_offset245, '../outputs/g_offset245.RDS')
print('Genetic offset ssp245 saved')
time <- Sys.time()
print(time)

g_offset585 <- genetic.offset(input = rbeta_data, 
                              env = env_lfmm, new.env = bio_values585, 
                              pop.labels = pop, K = 9)
print('Genetic offset ssp585 calculated')
time <- Sys.time()
print(time)

saveRDS(g_offset585, '../outputs/g_offset585.RDS')
print('Genetic offset ssp585 saved')
time <- Sys.time()
print(time)
print('Script complete')
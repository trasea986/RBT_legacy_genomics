library("tidyverse")
library("miceadds")
library("future.apply")
library("data.table")
print('packages loaded')

time_test_start = Sys.time()

#Load allele frequency matrix
allele_freq <- fread('../../outputs/no_koot/allele_freq.csv') #use this for subset of snps: ", select = c(1:100))"
pop_names <- allele_freq[,1]
allele_freq <- allele_freq[,-c(1)]

print('data loaded')

#Subset the maf2 dataframe based on population
LittleJacks <- allele_freq %>% slice(1)
BigJacks <- allele_freq %>% slice(2)
Duncan <- allele_freq %>% slice(3) 
Williams <- allele_freq %>% slice(4) 
Keithley <- allele_freq %>% slice(5) 
LittleWeser <- allele_freq %>% slice(6) 
Dry <- allele_freq %>% slice(7)
Fawn <- allele_freq %>% slice(8)
Mann1 <- allele_freq %>% slice(9)
Mann2 <- allele_freq %>% slice(10)


print('sliced')

plan(multisession, workers = 24) ## => parallelize on your local computer



print('starting loop')

#only simulating 20 ind here due to really large sample sizes

# start_time <- Sys.time()
# LJ <- list()
# LJ <- future_lapply(seq_along(LittleJacks), function(ii) {
  # j <- LittleJacks[[ii]]
  # tmp <- (rbeta(20, 
                # shape1= ((20 * as.numeric(LittleJacks[1, ii])) + 1), 
                # shape2= (20 * (1 - as.numeric(LittleJacks[1, ii])) + 1)))
  # tmp            ## same as return(tmp)
# }, future.seed = TRUE)
# print('LittleJacks loop complete')
# output_LJ <- as.data.frame(do.call(cbind, LJ))
# colnames(output_LJ) <- colnames(LittleJacks)
# print('LittleJacks bind and rename done')
# fwrite(output_LJ, "../../outputs/no_koot/LJ_sim.csv", row.names = FALSE)
# print('LittleJacks csv written')
# end_time <- Sys.time()
# time = end_time - start_time
# print(time)

# start_time <- Sys.time()
# BJ <- list()
# BJ <- future_lapply(seq_along(BigJacks), function(ii) {
  # j <- BigJacks[[ii]]
  # tmp <- (rbeta(20, 
                # shape1= ((20 * as.numeric(BigJacks[1, ii])) + 1), 
                # shape2= (20 * (1 - as.numeric(BigJacks[1, ii])) + 1)))
  # tmp            ## same as return(tmp)
# }, future.seed = TRUE)
# print('BigJacks loop complete')
# output_BJ <- as.data.frame(do.call(cbind, BJ))
# colnames(output_BJ) <- colnames(BigJacks)
# print('BigJacks bind and rename done')
# fwrite(output_BJ, "../../outputs/no_koot/BJ_sim.csv", row.names = FALSE)
# print('BigJacks csv written')
# end_time <- Sys.time()
# time = end_time - start_time
# print(time)

# start_time <- Sys.time()
# D <- list()
# D <- future_lapply(seq_along(Duncan), function(ii) {
  # j <- Duncan[[ii]]
  # tmp <- (rbeta(20, 
                # shape1= ((20 * as.numeric(Duncan[1, ii])) + 1), 
                # shape2= (20 * (1 - as.numeric(Duncan[1, ii])) + 1)))
  # tmp            ## same as return(tmp)
# }, future.seed = TRUE)
# print('Duncan loop complete')
# output_D <- as.data.frame(do.call(cbind, D))
# colnames(output_D) <- colnames(Duncan)
# print('Duncan bind and rename done')
# fwrite(output_D, "../../outputs/no_koot/D_sim.csv", row.names = FALSE)
# print('Duncan csv written')
# end_time <- Sys.time()
# time = end_time - start_time
# print(time)

# start_time <- Sys.time()
# W <- list()
# W <- future_lapply(seq_along(Williams), function(ii) {
  # j <- Williams[[ii]]
  # tmp <- (rbeta(20, 
                # shape1= ((20 * as.numeric(Williams[1, ii])) + 1), 
                # shape2= (20 * (1 - as.numeric(Williams[1, ii])) + 1)))
  # tmp            ## same as return(tmp)
# }, future.seed = TRUE)
# print('Williams loop complete')
# output_W <- as.data.frame(do.call(cbind, W))
# colnames(output_W) <- colnames(Williams)
# print('Method bind and rename done')
# fwrite(output_W, "../../outputs/no_koot/W_sim.csv", row.names = FALSE)
# print('Williams csv written')
# end_time <- Sys.time()
# time = end_time - start_time
# print(time)

# start_time <- Sys.time()
# K <- list()
# K <- future_lapply(seq_along(Keithley), function(ii) {
  # j <- Keithley[[ii]]
  # tmp <- (rbeta(20, 
                # shape1= ((20 * as.numeric(Keithley[1, ii])) + 1), 
                # shape2= (20 * (1 - as.numeric(Keithley[1, ii])) + 1)))
  # tmp            ## same as return(tmp)
# }, future.seed = TRUE)
# print('Keithley loop complete')
# output_K <- as.data.frame(do.call(cbind, K))
# colnames(output_K) <- colnames(Keithley)
# print('Keithley bind and rename done')
# fwrite(output_K, "../../outputs/no_koot/K_sim.csv", row.names = FALSE)
# print('Keithley csv written')
# end_time <- Sys.time()
# time = end_time - start_time
# print(time)

# start_time <- Sys.time()
# LW <- list()
# LW <- future_lapply(seq_along(LittleWeser), function(ii) {
  # j <- LittleWeser[[ii]]
  # tmp <- (rbeta(20, 
                # shape1= ((20 * as.numeric(LittleWeser[1, ii])) + 1), 
                # shape2= (20 * (1 - as.numeric(LittleWeser[1, ii])) + 1)))
  # tmp            ## same as return(tmp)
# }, future.seed = TRUE)
# print('LittleWeser loop complete')
# output_LW <- as.data.frame(do.call(cbind, LW))
# colnames(output_LW) <- colnames(LittleWeser)
# print('LittleWeser bind and rename done')
# fwrite(output_LW, "../../outputs/no_koot/LW_sim.csv", row.names = FALSE)
# print('LittleWeser csv written')
# end_time <- Sys.time()
# time = end_time - start_time
# print(time)

# start_time <- Sys.time()
# DR <- list()
# DR <- future_lapply(seq_along(Dry), function(ii) {
  # j <- Dry[[ii]]
  # tmp <- (rbeta(20, 
                # shape1= ((20 * as.numeric(Dry[1, ii])) + 1), 
                # shape2= (20 * (1 - as.numeric(Dry[1, ii])) + 1)))
  # tmp            ## same as return(tmp)
# }, future.seed = TRUE)
# print('Dry loop complete')
# output_DR <- as.data.frame(do.call(cbind, DR))
# colnames(output_DR) <- colnames(Dry)
# print('Dry bind and rename done')
# fwrite(output_DR, "../../outputs/no_koot/DR_sim.csv", row.names = FALSE)
# print('Dry csv written')
# end_time <- Sys.time()
# time = end_time - start_time
# print(time)

# start_time <- Sys.time()
# FA <- list()
# FA <- future_lapply(seq_along(Fawn), function(ii) {
  # j <- Fawn[[ii]]
  # tmp <- (rbeta(20, 
                # shape1= ((20 * as.numeric(Fawn[1, ii])) + 1), 
                # shape2= (20 * (1 - as.numeric(Fawn[1, ii])) + 1)))
  # tmp            ## same as return(tmp)
# }, future.seed = TRUE)
# print('Fawn loop complete')
# output_FA <- as.data.frame(do.call(cbind, FA))
# colnames(output_FA) <- colnames(Fawn)
# print('Fawn bind and rename done')
# fwrite(output_FA, "../../outputs/no_koot/FA_sim.csv", row.names = FALSE)
# print('Fawn csv written')
# end_time <- Sys.time()
# time = end_time - start_time
# print(time)

start_time <- Sys.time()
M1 <- list()
M1 <- future_lapply(seq_along(Mann1), function(ii) {
  j <- Mann1[[ii]]
  tmp <- (rbeta(20, 
                shape1= ((20 * as.numeric(Mann1[1, ii])) + 1), 
                shape2= (20 * (1 - as.numeric(Mann1[1, ii])) + 1)))
  tmp            ## same as return(tmp)
}, future.seed = TRUE)
print('Mann1 loop complete')
output_M1 <- as.data.frame(do.call(cbind, M1))
colnames(output_M1) <- colnames(Mann1)
print('Mann1 bind and rename done')
fwrite(output_M1, "../../outputs/no_koot/M1_sim.csv", row.names = FALSE)
print('Mann1 csv written')
end_time <- Sys.time()
time = end_time - start_time
print(time)

start_time <- Sys.time()
M2 <- list()
M2 <- future_lapply(seq_along(Mann2), function(ii) {
  j <- Mann2[[ii]]
  tmp <- (rbeta(20, 
                shape1= ((20 * as.numeric(Mann2[1, ii])) + 1), 
                shape2= (20 * (1 - as.numeric(Mann2[1, ii])) + 1)))
  tmp            ## same as return(tmp)
}, future.seed = TRUE)
print('Mann2 loop complete')
output_M2 <- as.data.frame(do.call(cbind, M2))
colnames(output_M2) <- colnames(Mann2)
print('Mann2 bind and rename done')
fwrite(output_M2, "../../outputs/no_koot/M2_sim.csv", row.names = FALSE)
print('Mann2 csv written')
end_time <- Sys.time()
time = end_time - start_time
print(time)

print('binding all together')
start_time <- Sys.time()
output_LJ <- fread("../../outputs/no_koot/LJ_sim.csv")
output_BJ <- fread("../../outputs/no_koot/BJ_sim.csv")
output_D <- fread("../../outputs/no_koot/D_sim.csv")
output_W <- fread("../../outputs/no_koot/W_sim.csv")
output_K <- fread("../../outputs/no_koot/K_sim.csv")
output_LW <- fread("../../outputs/no_koot/LW_sim.csv")
output_DR <- fread("../../outputs/no_koot/DR_sim.csv")
output_FA <- fread("../../outputs/no_koot/FA_sim.csv")
output_M1 <- fread("../../outputs/no_koot/M1_sim.csv")
output_M2 <- fread("../../outputs/no_koot/M2_sim.csv")

listofdf <- c('output_LJ', 'output_BJ', 'output_D', 'output_W', 'output_K', 'output_LW', 'output_DR', 'output_FA', 'output_M1', 'output_M2')

gen_all_pop <- rbindlist(mget(listofdf), idcol = TRUE)

print('done binding all pops')
end_time <- Sys.time()
time = end_time - start_time
print(time)

saveRDS(gen_all_pop, file = "../../outputs/no_koot/sim_data.rds")
print('RDS file saved')

gen_all_pop <- as.matrix(gen_all_pop)
print('converted to matrix')
fwrite(gen_all_pop, "../../outputs/no_koot/all_pop_sim.csv", row.names = FALSE)
print('output table with append complete')

exit <- function() {
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}

end_time_final <- Sys.time()
total_time <- end_time_final - time_test_start
print(total_time)

print("this is the last message")
exit()
print("you should not see this")

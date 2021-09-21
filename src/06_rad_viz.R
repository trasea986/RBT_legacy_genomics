#for viz
library(vegan)
loadRDS(rbt_rda, '../outputs/rbt_rda.RDS')

#Draft RDA analysis part 3
#Written by Yara Alshwairikh based on script from https://urldefense.com/v3/__http://popgen.nescent.org/2018-03-27_RDA_GEA.html__;!!JYXjzlvb!wZxPCa8dalpoF6uCmoDyuhjKOkrosRPJXJbN5Xo86BZ6G1epQVXfgdFzL05TNEWZgw$ 
#Last updated August 25th, 2020
#Script consists of 3 parts

###---Start script part 3---###

#A) Triplots
#read a "env" dataframe containing 1 column for "population", and 1 column for "type"
env_rda_plot <- read.csv("../outputs/env_rda.csv", stringsAsFactors = TRUE)

population <- env_rda_plot$Stream
type <- env_rda_plot$Ecotype

#set colors for each population
bg <- c("firebrick","brown3", "#orangered2", "lightcoral", "green3", "#mediumspringgreen", "#limegreen", "dodgerblue", "deepskyblue2", "lightgreen", "skyblue")

# axes 1 & 2
png(file="../outputs/RDA_plot1.png") 
plot(rbt_rda, type="n", scaling=3)
points(rbt_rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(rbt_rda, display="sites", pch=c(21, 24, 23) [type], cex=1.3, col="gray32", scaling=3, bg=bg[population]) # the populations (7 pops)
text(rbt_rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors (5)
legend("bottomright", legend=levels(population), bty="n", col="gray32", pch=c(21, 24, 23) [type], cex=1, pt.bg=bg)
#legend("bottomleft", legend=levels(type), bty="n", col="gray32", pch=c(1, 2), cex=1, pt.bg=type) 
dev.off()

#for viz


#Draft RDA analysis part 3
#Written by Yara Alshwairikh based on script from https://urldefense.com/v3/__http://popgen.nescent.org/2018-03-27_RDA_GEA.html__;!!JYXjzlvb!wZxPCa8dalpoF6uCmoDyuhjKOkrosRPJXJbN5Xo86BZ6G1epQVXfgdFzL05TNEWZgw$ 
#Last updated August 25th, 2020
#Script consists of 3 parts

###---Start script part 3---###

###---Plots---###
#Just like part 2, this part 3 depends on how many axes are significant and we want to plot
#Currently, this script assumes the first 3 axes are significant

#A) Triplots
#read a "env" dataframe containing 1 column for "population", and 1 column for "type"
env <- read.csv("env.csv")
#set levels
levels(env$population) <- c("Yakima","Wenatchee","Lyons Ferry","Methow","Deschutes","Clearwater", "Priest")
levels(env$type) <- c("Fall","Summer")

population <- env$population
type <- env$type

#set colors for each population
bg <- c("#A50026","#DC3D2D", "#F57D4A", "#354B99", "#FDB366", "#6DA5CC", "#FED98B")

# axes 1 & 2
png(file="RDA_plot1.png") 
plot(ots.rda, type="n", scaling=3)
points(ots.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)           # the SNPs
points(ots.rda, display="sites", pch=c(21, 24) [type], cex=1.3, col="gray32", scaling=3, bg=bg[population]) # the populations (7 pops)
text(ots.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors (5)
legend("bottomright", legend=levels(population), bty="n", col="gray32", pch=c(21, 24) [type], cex=1, pt.bg=bg)
legend("bottomleft", legend=levels(type), bty="n", col="gray32", pch=c(1, 2), cex=1, pt.bg=type) 
dev.off()
#note: the shapes in the legend by type are incorrect. Clearwater should be a triangle instead of Wenatchee

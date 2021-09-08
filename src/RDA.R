#Draft RDA analysis for LG2020 group project
#Written by Yara Alshwairikh based on script from https://urldefense.com/v3/__http://popgen.nescent.org/2018-03-27_RDA_GEA.html__;!!JYXjzlvb!xWdlbcC9VtEpswPuczv5tejdyniCdOyEP1f3Vqkr4Mzo369sWSVKn6Uwmk_2yHlaCA$ 
#Last updated June 29th, 2020

#Load libraries
library("tidyverse")
library("dplyr")
library("psych")
library("vegan")
library("miceadds")

#Import environmental data
env_combination <- read.csv("env_combination.csv")
env_site <- read.csv("env_site.csv")
env_migration <- read.csv("env_migration.csv")

load.Rdata(filename = "freqs.RData", "maf") #import the .Rdata file containing allele frequencies

sum(is.na(maf)) # confirm there are no NAs in the matrix

#maf <- maf[, sample(ncol(maf), 100000)] #Sample random 100,000 SNPs from the matrix
-------------
  
#Rearrange rows and add rownames to match the allele frequency matrix
#Original order from Steve is : Methow, Wentachee, Deschutes, Yakima, Priest, Clearwater, Lyons
#The maf file is: Yakima, Wenatchee, Lyons Ferry, Methow, Deschutes, Clearwater, Priest
#Create new column "Population" and assign each row a letter 
#rearrange alphabetrically to end up with the same order as the allele frequency matrix

env_combination <- env_combination %>% 
  mutate(Population = c("d", "b", "e", "a", "g", "f", "c")) %>% arrange(Population) 
rownames(env_combination)<- 
  c("Yakima","Wenatchee","Lyons Ferry","Methow",
    "Deschutes","Clearwater","Priest")

env_site <- env_site %>% 
  mutate(Population = c("d", "b", "e", "a", "g", "f", "c")) %>% arrange(Population)
rownames(env_site)<- 
  c("Yakima","Wenatchee","Lyons Ferry","Methow",
                       "Deschutes","Clearwater","Priest")

env_migration <- env_migration %>% 
  mutate(Population = c("d", "b", "e", "a", "g", "f", "c")) %>% arrange(Population)
rownames(env_migration)<- 
  c("Yakima","Wenatchee","Lyons Ferry","Methow",
    "Deschutes","Clearwater","Priest")

#Remove the population column
env_combination <- subset(env_combination, select = -Population) 
env_site <- subset(env_site, select = -Population) 
env_migration <- subset(env_migration, select = -Population) 

#Confirm that genotypes and environmental data are in the same order
identical(rownames(maf), rownames(env_combination))
identical(rownames(maf), rownames(env_site))
identical(rownames(maf), rownames(env_migration))

#convert all predictors to numeric for RDA to work
env_combination_num <- env_combination %>% mutate_all(as.numeric)
env_site_num <- env_site %>% mutate_all(as.numeric)
env_migration_num <- env_migration %>% mutate_all(as.numeric)

----------------------
#Run RDA for COMBINATION
ots.rda <- vegan::rda(maf ~ ., data=env_combination_num, Scale=T)
ots.rda

vegan::RsquareAdj(ots.rda)

summary(ots.rda)$concont 
screeplot(ots.rda)

anova.cca(ots.rda, by="axis") 
vegan::vif.cca(ots.rda)

signif.full <- anova.cca(ots.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full

signif.axis <- anova.cca(ots.rda, by="axis", parallel=getOption("mc.cores"))
signif.axis

vif.cca(ots.rda)

plot(ots.rda, scaling=3)          # default is axes 1 and 2
plot(ots.rda, choices = c(1,3), scaling=3)  # axes 1 and 3

load.rda <- summary(ots.rda)$species[,1:3]
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 
cand3 <- outliers(load.rda[,3],3) 

ncand <- length(cand1) + length(cand2) + length(cand3)
ncand 

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
 
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=11)  # 11 columns for 11 predictors
colnames(foo) <- c("Bio10_rang", "Slp_mean", "SRad_mean", "Wsp_range", 
                   "Wsp_site", "stOrd_mean", "stOrd_site", "HLI_Site",
                   "Bio3_site", "Bio5_site", "Bio5_range") 

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- maf[,nam]
  foo[i,] <- apply(env_combination_num,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$snp[duplicated(cand$snp)])  

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) 

table(foo[foo[,1]==2,2])
table(foo[foo[,1]==3,2]) 

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,15] <- names(which.max(abs(bar[3:14]))) # gives the variable
  cand[i,16] <- max(abs(bar[3:14]))              # gives the correlation
}
#the 3:14 are the columns where the environmental variables are
#can confirm through str(cand). MUST change cand[i,#] to specificy the column right
#after 14, which here is column 15, or else you get an error 
#do the same thing for the correlation cand[ ]

colnames(cand)[15] <- "predictor"
colnames(cand)[16] <- "correlation"

table(cand$predictor) 

#Bio10_rang  Bio3_site Bio5_range  Bio5_site   HLI_Site   Slp_mean  SRad_mean 
#13998      18844       2053       5453        577       3216       4915 
#stOrd_mean stOrd_site  Wsp_range   Wsp_site 
#10924       1945      19057       9165 

#most SNPs associated with Wsp_range, Bio3_site, Bio10_rang

sel <- cand$snp
env <- cand$predictor
env[env=="Bio10_rang"] <- '#1f78b4'
env[env=="Bio3_site"] <- '#a6cee3'
env[env=="Bio5_range"] <- '#6a3d9a'
env[env=="Bio3_site"] <- '#e31a1c'
env[env=="Bio5_site"] <- '#33a02c'
env[env=="HLI_Site"] <- '#ffff33'
env[env=="Slp_mean"] <- '#fb9a99'
env[env=="SRad_mean"] <- 'green'
env[env=="stOrd_mean"] <- 'blue'
env[env=="stOrd_site"] <- 'red'
env[env=="Wsp_range"] <- 'yellow'
env[env=="Wsp_site"] <- 'orange'

# color by predictor:
col.pred <- rownames(ots.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("chr",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33','#fb9a99', 'green', 'blue',
        'red', 'yellow', 'orange')

# axes 1 & 2
plot(ots.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(ots.rda, display="species", pch=21, cex=1, col="green", bg=col.pred, scaling=3)
points(ots.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(ots.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("Bio10_rang", "Slp_mean", "SRad_mean", "Wsp_range", 
                               "Wsp_site", "stOrd_mean", "stOrd_site", "HLI_Site"),
                                bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(ots.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(ots.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(ots.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(ots.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=c("Bio10_rang", "Slp_mean", "SRad_mean", "Wsp_range", 
                               "Wsp_site", "stOrd_mean", "stOrd_site", "HLI_Site"), 
                                bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
------------------
#Run RDA for SITE
ots.rda2 <- vegan::rda(maf ~ ., data=env_site_num, Scale=T)
ots.rda2

vegan::RsquareAdj(ots.rda2)

summary(ots.rda2)$concont 
screeplot(ots.rda2)

anova.cca(ots.rda2, by="axis") 
vegan::vif.cca(ots.rda2)

plot(ots.rda2, scaling=3)          # default is axes 1 and 2
plot(ots.rda2, choices = c(1,3), scaling=3)  # axes 1 and 3


load.rda2 <- summary(ots.rda2)$species[,1:3]
hist(load.rda2[,1], main="Loadings on RDA1")
hist(load.rda2[,2], main="Loadings on RDA2")
hist(load.rda2[,3], main="Loadings on RDA3") 


outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand11 <- outliers(load.rda2[,1],3) 
cand22 <- outliers(load.rda2[,2],3) 
cand33 <- outliers(load.rda2[,3],3)

ncand2 <- length(cand11) + length(cand22) + length(cand33)
ncand2 

cand11 <- cbind.data.frame(rep(1,times=length(cand11)), names(cand11), unname(cand11))
cand22 <- cbind.data.frame(rep(2,times=length(cand22)), names(cand22), unname(cand22))
cand33 <- cbind.data.frame(rep(3,times=length(cand33)), names(cand33), unname(cand33))

colnames(cand11) <- colnames(cand22) <- colnames(cand33) <- c("axis","snp","loading")

cand <- rbind(cand11, cand22, cand33)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=7)  # 7 columns for 7 predictors
colnames(foo) <- c("Bio11_site", "Wsp_site", "HLI_Site", "Rough_Site", "Bio15_site", "Bio3_site","Bio5_site") 

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- maf[,nam]
  foo[i,] <- apply(env_site_num,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$snp[duplicated(cand$snp)])  

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2])

table(foo[foo[,1]==2,2]) 
table(foo[foo[,1]==3,2]) 

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,11] <- names(which.max(abs(bar[4:10]))) # gives the variable
  cand[i,12] <- max(abs(bar[4:10]))              # gives the correlation
}

colnames(cand)[11] <- "predictor"
colnames(cand)[12] <- "correlation"

table(cand$predictor) 

#Bio11_site Bio15_site  Bio3_site  Bio5_site   HLI_Site Rough_Site   Wsp_site 
#1246         23       8426          1         50        328        127 
#most SNPs associated with Bio3_site, Wsp_site, Bio15_site 

sel <- cand$snp
env <- cand$predictor
env[env=="Bio11_site"] <- '#1f78b4'
env[env=="Bio15_site"] <- '#a6cee3'
env[env=="Bio3_site"] <- '#6a3d9a'
env[env=="Bio5_site"] <- '#e31a1c'
env[env=="HLI_Site"] <- '#33a02c'
env[env=="Rough_Site"] <- '#ffff33'
env[env=="Wsp_site"] <- '#fb9a99'

# color by predictor:
col.pred <- rownames(ots.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("chr",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33','#fb9a99')

# axes 1 & 2
plot(ots.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(ots.rda, display="species", pch=21, cex=1, col="green", bg=col.pred, scaling=3)
points(ots.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(ots.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("Bio10_rang", "Slp_mean", "SRad_mean", "Wsp_range", 
                               "Wsp_site", "stOrd_mean", "stOrd_site", "HLI_Site"),
       bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(ots.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(ots.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(ots.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(ots.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=c("Bio10_rang", "Slp_mean", "SRad_mean", "Wsp_range", 
                               "Wsp_site", "stOrd_mean", "stOrd_site", "HLI_Site"), 
       bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

-------------------
  
#Run RDA for MIGRATION
ots.rda3 <- vegan::rda(maf ~ ., data=env_migration_num, Scale=T)
ots.rda3

vegan::RsquareAdj(ots.rda3)

summary(ots.rda3)$concont 
screeplot(ots.rda3)

anova.cca(ots.rda3, by="axis") 
vegan::vif.cca(ots.rda3)

plot(ots.rda3, scaling=3)          # default is axes 1 and 2
plot(ots.rda3, choices = c(1,3), scaling=3)  # axes 1 and 3


load.rda3 <- summary(ots.rda3)$species[,1:3]
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 


outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand111 <- outliers(load.rda3[,1],3) 
cand222 <- outliers(load.rda3[,2],3) 
cand333 <- outliers(load.rda3[,3],3) 

ncand <- length(cand111) + length(cand222) + length(cand333)
ncand 

cand111 <- cbind.data.frame(rep(1,times=length(cand111)), names(cand111), unname(cand111))
cand222 <- cbind.data.frame(rep(2,times=length(cand222)), names(cand222), unname(cand222))
cand333 <- cbind.data.frame(rep(3,times=length(cand333)), names(cand333), unname(cand333))

colnames(cand111) <- colnames(cand222) <- colnames(cand333) <- c("axis","snp","loading")

cand <- rbind(cand111, cand222, cand333)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=8)  # 8 columns for 8 predictors
colnames(foo) <- c("DEM_mean", "DEM_range", "Bio1_mean", "Bio15_mean", "RvT_range", "Rough_mean",
                   "mig_distance", "Bio5_range") 

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- maf[,nam]
  foo[i,] <- apply(env_migration_num,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$snp[duplicated(cand$snp)])  

foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) 

table(foo[foo[,1]==2,2]) 
table(foo[foo[,1]==3,2]) 

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,12] <- names(which.max(abs(bar[4:11]))) # gives the variable
  cand[i,13] <- max(abs(bar[4:11]))              # gives the correlation
}

colnames(cand)[12] <- "predictor"
colnames(cand)[13] <- "correlation"

table(cand$predictor) 

#Bio1_mean   Bio15_mean   Bio5_range     DEM_mean    DEM_range mig_distance 
#8          113         8295           33         1110         5853 
#Rough_mean    RvT_range 
#1866          196 

#most SNPs associated with Bio5_range, mig_distance, Rough_mean


sel <- cand$snp
env <- cand$predictor
env[env=="Bio1_mean"] <- '#1f78b4'
env[env=="Bio15_mean"] <- '#a6cee3'
env[env=="Bio5_range"] <- '#6a3d9a'
env[env=="DEM_mean"] <- '#e31a1c'
env[env=="DEM_range"] <- '#33a02c'
env[env=="mig_distance"] <- '#ffff33'
env[env=="Rough_mean"] <- '#fb9a99'
env[env=="RvT_range"] <- 'green'

# color by predictor:
col.pred <- rownames(ots.rda$CCA$v) # pull the SNP names

for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("chr",col.pred)] <- '#f1eef6' # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#a6cee3','#6a3d9a','#e31a1c','#33a02c','#ffff33','#fb9a99', 'blue')

# axes 1 & 2
plot(ots.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(ots.rda, display="species", pch=21, cex=1, col="green", bg=col.pred, scaling=3)
points(ots.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(ots.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("Bio10_rang", "Slp_mean", "SRad_mean", "Wsp_range", 
                               "Wsp_site", "stOrd_mean", "stOrd_site", "HLI_Site"),
       bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

# axes 1 & 3
plot(ots.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(1,3))
points(ots.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3, choices=c(1,3))
points(ots.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3, choices=c(1,3))
text(ots.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(1,3))
legend("bottomright", legend=c("Bio10_rang", "Slp_mean", "SRad_mean", "Wsp_range", 
                               "Wsp_site", "stOrd_mean", "stOrd_site", "HLI_Site"), 
       bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

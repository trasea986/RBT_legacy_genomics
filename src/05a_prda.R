library("vegan")
library("miceadds")
library("data.table")
#library("psych") if checking correlations again

#Load allele frequency matrix
allele_freq <- fread('../outputs/allele_freq.csv')

allele_freq <- allele_freq[,-c(1)]
allele_freq <- as.matrix(allele_freq, ncolumns = ncol(allele_freq))

#use csv to delineate col as numeric
env_mat <- read.csv('../outputs/env_rda.csv')
env_mat$Stream <- as.character(env_mat$Stream)

#remove ecotype and pop column to make matching matrix. also remove elevation due to R2 issue.
env_mat <- env_mat[,-c(1:3)]

#pairs.panels(env_mat, Scale = T)

#add matching row names
rownames(allele_freq) <- rownames(env_mat)

allele_freq1 <- allele_freq[,1:500]

rbt_pca <- rda(allele_freq1, scale=T)
labels(rbt_pca)

rbt_pca_PC1 <- scores(rbt_pca, choices=1, display="sites", scaling=0)

rbt_rda <- vegan::rda(allele_freq ~ . + Condition(rbt_pca_PC1), data = env_mat, Scale = T)

rbt_rda

#write rda object
#saveRDS(rbt_rda, '../outputs/rbt_prda.RDS')

rbt_rsquare <- vegan::RsquareAdj(rbt_rda) #R squared values
rbt_rsquare
#write.csv(rbt_rsquare, '../outputs/prbt_rsquare.csv')

rbt_prop_variance <- summary(rbt_rda)$concont #proportion of variance explaines by each axis
rbt_prop_variance
#screeplot(ots.rda) #visualize the canconical eignevalues
#write.csv(rbt_prop_variance, '../outputs/prbt_prop_variance.csv')

#Check the FULL RDA model for significance
signif_full <- anova.cca(rbt_rda, parallel=getOption("mc.cores")) # default is permutation=999
signif_full
saveRDS(signif_full, file = "../outputs/psignif_full.RDS")

#Check each constrained axis for significance
signif_axis <- anova.cca(rbt_rda, by="axis", parallel=getOption("mc.cores"))
signif_axis #Check to see which axes have p = 0.05. Ideally this should match the results in the screeplots
saveRDS(signif_axis, file = "../outputs/psignif_axis.RDS")
write.csv(signif_axis, file = "../outputs/psignif_axis.csv")

#Check each term for significance
signif_term <- anova.cca(rbt_rda, by="terms", parallel=getOption("mc.cores"))
signif_term #Check to see which axes have p = 0.05. Ideally this should match the results in the screeplots
saveRDS(signif_term, file = "../outputs/psignif_term.RDS")
write.csv(signif_term, file = "../outputs/psignif_term.csv")

#Find candidate SNPs: this part relies on knowledge of which axes are significant. 
#Currently, the script is written assuming axes 1, 2, and 3 are significant
load_rda <- summary(rbt_rda)$species[,1:3]
#hist(load_rda[,1], main="Loadings on RDA1")
#hist(load.rda[,2], main="Loadings on RDA2")
#hist(load.rda[,3], main="Loadings on RDA3") 

saveRDS(load_rda, file = "../outputs/load_prda.RDS")

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
} #function to find the outlier SNPs


#apply the outliers() function to each axis
cand1 <- outliers(load_rda[,1],3.5) #second number is number of stdeviations
#cand2 <- outliers(load.rda[,2],3) 
#cand3 <- outliers(load.rda[,3],3) 

ncand <- length(cand1) #+ length(cand2) + length(cand3)
ncand #total number of candidate SNPs

#Organzie the results into a datadram with the axis, SNP, laoding, a correlation with each environmental variable
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
#cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
#cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

#this was originally chained where 2 pointed to 1, 3 pointed to 2, and then names
colnames(cand1) <- c("axis","snp","loading")#<- colnames(cand2) <- colnames(cand3) 

cand <- rbind(cand1) #, cand2, cand3)
cand$snp <- as.character(cand$snp)



###################### ID env associations

foo <- matrix(nrow=(ncand), ncol=9)  # create 5 columns for 5 predictors
colnames(foo) <- c("CANOPY","SLOPE","CUMDRAINAG","Precipitation.Seasonality","Mean.Diurnal.Temperature.Range","Isothermality","Min.Temperature.of.Coldest.Month","Temperature.Annual.Range","Stream_temp_93_11")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- allele_freq[,nam]
  foo[i,] <- apply(env_mat,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  

n_snps <- length(cand$snp[duplicated(cand$snp)]) #check for duplicate SNPs


#foo <- cbind(cand$axis, duplicated(cand$snp)) 
#table(foo[foo[,1]==1,2]) #check for duplicates on axis 1
#table(foo[foo[,1]==2,2]) #check for duplicates on axis 2
#table(foo[foo[,1]==3,2]) #check for duplicates on axis 3

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections

#To find which of the predictors each candidate SNP is most strongly correlated with:
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,13] <- names(which.max(abs(bar[4:12]))) # gives the variable
  cand[i,14] <- max(abs(bar[4:12]))              # gives the correlation
}
#note: the 4:8 are the columns where the environmental variables are, must manually adjust this for every dataset
#can confirm through str(cand). MUST change cand[i,#] to specify the column right after col 8, 
#which here is column 9, or else you get an error 
#do the same thing for the correlation cand[ ]

#assign column names
colnames(cand)[13] <- "predictor"
colnames(cand)[14] <- "correlation"

#table(cand$predictor) #lists top associations 

#Write the output into a .csv
write.csv(cand, "../outputs/pRDA_all_SNP_cor.csv", row.names = FALSE)

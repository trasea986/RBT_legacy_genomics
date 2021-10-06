#for viz
library(tidyverse)
library(reshape2)
library(vegan)
library(qqman)
library(IHW)
library(robust)
library(data.table)
library(qvalue) #note: bioconductor package BiocManager::install("qvalue")




# rda correlation plot ----------------------------------------------------
#load in output files from rda
df_cand <- read.csv("../outputs/RDA_cand_3.5.csv")
df_all <- read.csv("../outputs/RDA_all_SNP_cor.csv")

#only select desired columns
df_cand <- df_cand %>%
  select(snp, correlation)
df_all <- df_all %>%
  select(snp, correlation)

#this section is for testing purposes to sample random snps from df_all
#index <- sample(1:nrow(df_all), 4000)
#df_sub <- df_all[index,]
#add in the candidates

df_all <- rbind(df_all, df_cand)

#add in column to flag if candidate or not. not that you ahve to run duplicate both directions for both instances to be marked as true
df_all$cand <- duplicated(df_all[,1]) | duplicated(df_all[,1], fromLast=TRUE)

#delete the duplicates but retain the new column
df <- df_all %>% unique()

print('ID the outliers')

#for plotting, need to subset the strangely named chromosomes, format, and then recombine
sub2 <- df %>% filter(str_count(df$snp, "_") == 2)
sub3 <- df %>% filter(str_count(df$snp, "_") == 3)


#separate the snp column to get chr and bp
newCols1 <- colsplit(sub2$snp, "_", c("chr_name","BP", "junk"))
newCols1$junk <- NULL
newCols2 <- colsplit(newCols1$chr_name, "y", c("junk", "chr"))
newCols2$junk <- NULL
df_pval_sub2 <- cbind(sub2, newCols1, newCols2)

#same as above but for those with three _ in the name
newCols1 <- colsplit(sub3$snp, "_", c("chr_name","junk", "BP", "junk1"))
newCols1$junk <- NULL
newCols1$junk1 <- NULL
newCols2 <- colsplit(newCols1$chr_name, "y", c("junk", "chr"))
newCols2$junk <- NULL
df_pval_sub3 <- cbind(sub3, newCols1, newCols2)

#combine for plotting and format columns as numeric
df_pval_plot <-rbind(df_pval_sub2, df_pval_sub3)
df_pval_plot$BP <- as.numeric(df_pval_plot$BP)
df_pval_plot$chr <- as.numeric(df_pval_plot$chr)


#calculate the cumulative snp position function here for ggplot
df_ggplot <- df_pval_plot %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(df_pval_plot, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, BP) %>%
  mutate( BPcum=BP+tot)

#remove NW chromosomes
df_ggplot <- df_ggplot %>% filter (chr_name != "NW")

print('Finished prepping ggplot file and formatting snp names')

#save input for ggplot
write.csv(df_ggplot, "../outputs/figures/final_rda_plot.csv", row.names=FALSE)

#Then we need to prepare the X axis. Want to display the cumulative position of SNP in bp, but just show the chromosome number
axisdf = df_ggplot %>% 
  group_by(chr) %>% 
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#calculate minimum correlation that was above the threshold
cor_thresh <- df_ggplot %>%
  group_by(cand) %>%
  summarize(thresh=min(correlation))

#ggplot
man_plot <- ggplot(df_ggplot, aes(x=BPcum, y=correlation)) +
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("darkgrey", "dodgerblue4"), 32 )) +
  geom_hline(yintercept=as.numeric(cor_thresh[2,2]), color = "darkred") +
  scale_x_continuous(label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  ggtitle('RDA Correlation Value') +
  ylab("Maximum Environmental \n Correlation")+
  theme_classic(base_size = 14, base_family = "Times") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = "../outputs/figures/RDA_man.png", plot = man_plot, height = 6, width = 10, units = "in")

print(paste('Saved plot of RDA'))

#ggplot
man_plot <- 
  df_ggplot %>%
  filter(cand == "TRUE") %>%
  ggplot(aes(x=BPcum, y=correlation)) +
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("darkgrey", "dodgerblue4"), 32 )) +
  geom_hline(yintercept=as.numeric(cor_thresh[2,2]), color = "darkred") +
  scale_x_continuous(label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  ggtitle('RDA Correlation Value') +
  ylab("Maximum Environmental \n Correlation")+
  theme_classic(base_size = 14, base_family = "Times") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = "../outputs/figures/RDA_man_cands.png", plot = man_plot, height = 6, width = 10, units = "in")

print(paste('Saved plot of RDA correlations'))



# pvalue ------------------------------------------------------------------
#based on following scripts: https://www.biorxiv.org/content/biorxiv/suppl/2018/02/03/258988.DC1/258988-4.pdf
#see also https://www.biorxiv.org/content/10.1101/258988v2.full scripts

rbt_rda <- readRDS('../outputs/rbt_rda.RDS')

#function for p values
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

#look at inertia. potentially use k = 4... but only first axis significant so going with minimum needed (2)
ggplot() +
  geom_line(aes(x=c(1:length(rbt_rda$CCA$eig)), y=as.vector(rbt_rda$CCA$eig)), linetype="dotted",
            size = 1.5, color="darkgrey") +
  geom_point(aes(x=c(1:length(rbt_rda$CCA$eig)), y=as.vector(rbt_rda$CCA$eig)), size = 3,
             color="darkgrey") +
  scale_x_discrete(name = "Ordination axes", limits=c(1:9)) +
  ylab("Inertia") +
  theme_bw()

res_rdadapt<-rdadapt(rbt_rda, 2)
#saveRDS(res_rdadapt, '../outputs/res_rdadapt.RDS')
res_rdadapt <- readRDS('../outputs/res_rdadapt.RDS')

#need to add back in the names and then do separation as above

#Load allele frequency matrix
allele_freq <- fread('../outputs/allele_freq.csv')

allele_freq <- allele_freq[,-c(1)]

#add in the column names
res_rdadapt$snp <- colnames(allele_freq)

#next steps are same as above for splitting snp names for actual positions
sub2 <- res_rdadapt %>% filter(str_count(res_rdadapt$snp, "_") == 2)
sub3 <- res_rdadapt %>% filter(str_count(res_rdadapt$snp, "_") == 3)

newCols1 <- colsplit(sub2$snp, "_", c("chr_name","BP", "junk"))
newCols1$junk <- NULL
newCols2 <- colsplit(newCols1$chr_name, "y", c("junk", "chr"))
newCols2$junk <- NULL
df_pval_sub2 <- cbind(sub2, newCols1, newCols2)

newCols1 <- colsplit(sub3$snp, "_", c("chr_name","junk", "BP", "junk1"))
newCols1$junk <- NULL
newCols1$junk1 <- NULL
newCols2 <- colsplit(newCols1$chr_name, "y", c("junk", "chr"))
newCols2$junk <- NULL
df_pval_sub3 <- cbind(sub3, newCols1, newCols2)

df_pval_plot <-rbind(df_pval_sub2, df_pval_sub3)
df_pval_plot$BP <- as.numeric(df_pval_plot$BP)
df_pval_plot$chr <- as.numeric(df_pval_plot$chr)

df_ggplot <- df_pval_plot %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(df_pval_plot, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, BP) %>%
  mutate( BPcum=BP+tot)

df_ggplot <- df_ggplot %>% filter (chr_name != "NW")
print('Finished prepping ggplot file and formatting snp names')
write.csv(df_ggplot, "../outputs/figures/final_rda_plot_pvals.csv", row.names=FALSE)

axisdf = df_ggplot %>% 
  group_by(chr) %>% 
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#for pvalue lines, we can do the same calculations as lfmm
##calculate BH threshold + #Bonferroni correction
bonferroni <- .05/(nrow(df_ggplot))
print("Bonferroni criticla value at 0.05")

#bh test for each variable
criticalValue <- get_bh_threshold(df_ggplot$p.values, 0.05)
print("BH criticla value at 0.05")


##main manhattan plot
man_plot <- ggplot(df_ggplot, aes(x=BPcum, y=-log10(p.values))) +
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("darkgrey", "dodgerblue4"), 32 )) +
  geom_hline(yintercept=-log10(bonferroni), color = "darkred") + #line for Bonferroni correction
  geom_hline(yintercept=-log10(criticalValue), color = "blue") + #line for BH correction
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  ggtitle('RDA, Chi-squared p-values') +
  theme_classic(base_size = 14, base_family = "Times") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = "../outputs/figures/RDA_plot_p.png", plot = man_plot, height = 6, width = 10, units = "in")

#can do the same with q values. note: not sure on the threshold here that would be appropriate?
man_plot <- ggplot(df_ggplot, aes(x=BPcum, y=-log10(q.values))) +
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("darkgrey", "dodgerblue4"), 32 )) +
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  ggtitle('RDA, q-values') +
  theme_classic(base_size = 14, base_family = "Times") +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = "../outputs/figures/RDA_plot_q.png", plot = man_plot, height = 6, width = 10, units = "in")

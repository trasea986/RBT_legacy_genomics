#for viz
library(tidyverse)
library(reshape2)
library(vegan)
library(qqman)
library(IHW)

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

print('Finished prepping ggplot file')

#save input for ggplot
write.csv(df_ggplot, paste("../outputs/figures/final_rda_plot_",df[1,3],".csv", sep = ''), row.names=FALSE)

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


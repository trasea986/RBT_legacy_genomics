library(tidyverse)
library(reshape2)
library(vegan)
library(qqman)

#bring in outputs from 08 and row bind them together

####Clean and process data for Manhattan plots
#Load the allele frequency matrix with 20 simulated individuals per population
rbeta_data <- readRDS("../../outputs/no_koot/sim_data.rds")
rbeta_data$.id <- NULL #remove .id column that gets generated when reading the data


rbeta_data <- rbeta_data[, 1:50]


env_lfmm <- read.csv('../../outputs/no_koot/env_rda.csv')

#duplicate for each of the 20 individuals
# Create an index of the rows you want with duplications
env_lfmm <- env_lfmm[rep(seq_len(nrow(env_lfmm)), each = 20), ]

#remove columns that do not have the env variables
env_lfmm <- env_lfmm[,-c(1:3)]

print('loaded env and genomic data in')

#rename columns into "X" for chromsoomes, and "pval"
#colnames(df) <- c('X', 'pval')
df <- read.csv('../../outputs/no_koot/allraw_lfmm1_K9.csv')
df <- read.csv('../../outputs/no_koot/allraw_lfmm2_K9.csv')
df <- read.csv('../../outputs/no_koot/allraw_lfmm3_K9.csv')
df <- read.csv('../../outputs/no_koot/allraw_lfmm4_K9.csv')
df <- read.csv('../../outputs/no_koot/allraw_lfmm5_K9.csv')
df <- read.csv('../../outputs/no_koot/allraw_lfmm6_K9.csv')
df <- read.csv('../../outputs/no_koot/allraw_lfmm7_K9.csv')
df <- read.csv('../../outputs/no_koot/allraw_lfmm8_K9.csv')
df <- read.csv('../../outputs/no_koot/allraw_lfmm9_K9.csv')


#print("renames columns into X for chromsoomes, and pval")

##calculate BH threshold + #Bonferroni correction
bonferroni <- .05/(nrow(df))

#bh test for each variable
#create vector of ranks from the vector of p values
p_rank <- rank(df, na.last = TRUE, ties.method = c("average"))
#bring these two vectors together
data_ex <- data.frame(df, p_rank)
#calculate when p(i) <= rank(i) / k * 0.05. k calculated from number of rows function.
data_ex$pval_less <- ((data_ex$df.pval) <= (data_ex$p_rank/nrow(df$pval) * 0.05))
#summarize this new true/false column and find the maximum rank for a return of TRUE or False from previous step. this returns "10" for the TRUE group
data_ex %>%
  group_by(pval_less) %>%
  summarize(max_rank = max(p_rank))
#then call the row with the rank returned in the last line of code, in this case 10
bh <- data_ex[p_rank == 10,]
bh <- bh[,1]

#add back in SNP names
snps <- colnames(rbeta_data)
df <- cbind(snps, df)

names(df)[names(df) == 'X'] <- 'pval'

#separate the snp column (here, it is X) into "chr_number" and "snp_position" columns
newColNames2 <- c("chr","BP", "junk")
newCols2 <- colsplit(df$snps, "_", newColNames2)
newCols2$junk <- NULL
newCols2$chr <- as.numeric(strsplit(newCols2$chr, "\\D+")[[1]][-1])

df_pval_plot <- cbind(df, newCols2)

df_pval_plot$BP <- as.numeric(df_pval_plot$BP) #convert int to numeric
df_pval_plot$chr <- as.numeric(df_pval_plot$chr) #convert int to numeric

print("created chr column and BP column")

#write the fully processed data into a .csv, this format is ready for the Manhattan plot script
write.csv(df, "../../outputs/no_koot/final_lfmm_K9.csv") ###EDIT ME###

cand_SNP_K9 <- subset(df, pval <= bonferroni)
write.csv(cand_SNP_K9, "cand_SNP_lfmm1_K9_bon.csv", row.names=FALSE)
cand_SNP_K9 <- subset(df, pval <= bh)
write.csv(cand_SNP_K9, "cand_SNP_lfmm1_K9_bh.csv", row.names=FALSE)


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


#Then we need to prepare the X axis. Want to display the cumulative position of SNP in bp, but just show the chromosome name instead

axisdf = df_ggplot %>% 
  group_by(chr) %>% 
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )



#function to create plots and save each one as a .png
plotstuff <- function(dataset) {
  plot_title <- paste("Plot of LFMM, K=9")
  ggplot(dataset, aes(x=BPcum, y=-log10(pval))) + #change pval to score if desired
    # Show all points
    geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    #add corrections lines
    geom_hline(yintercept=-log10(bonferroni), color = "darkred") + #line for Bonferroni correction, calculated manually by dividing 0.05 by the total number of SNPs (0.05/3212126 = 1.187049e-08)
    #geom_hline(yintercept=-log10(criticalValue), color = "blue") + #line for BH correction as calculated in part 2
    # custom X axis:
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    # Custom the theme:
    ggtitle(plot_title) +
    theme_classic() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  ggsave(filename = paste(plot_title,".png", sep = ""), plot = last_plot())
}


#Create the Manhattan plot 
plotstuff(dataset = df_ggplot)
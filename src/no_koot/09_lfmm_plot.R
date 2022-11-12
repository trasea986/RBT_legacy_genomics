library(tidyverse)
library(reshape2)
library(vegan)
library(qqman)
library(IHW)

#bring in outputs from 08 and row bind them together

####Clean and process data for Manhattan plots
#Load the allele frequency matrix with 20 simulated individuals per population
rbeta_data <- readRDS("../../outputs/no_koot/sim_data.rds")
rbeta_data$.id <- NULL #remove .id column

###looks like one (the first one) output the columns are different. below code fixes column order.
###note to self to double check output if running different k values in the future
###only run this the first time
fix <- read.csv(files[1]) #noticed that the first file x and y columns were swapped
fix1 <- cbind(fix[,2],fix[,1])
colnames(fix1) <- c('', 'env')
write.csv(fix1, paste(files[1]), row.names = FALSE)

print('loaded env and genomic data in')
data_dir <- '../../outputs/no_koot/'
files <- list.files(data_dir, "allraw_*", full.names = TRUE)

print('created file list, running lapply')

lapply(files, FUN = function(i) {
  
  #load in output files from lfmm
  df <- read.csv(i)
  
  #rename columns into "X" for chromsoomes, and "pval"
  colnames(df) <- c('pval', 'env_var')
  
  #print("renames columns into X for chromsoomes, and pval")
  
  ##calculate BH threshold + #Bonferroni correction
  bonferroni <- .05/(nrow(df))
  print("Bonferroni criticla value at 0.05")
  
  #bh test for each variable
  criticalValue <- get_bh_threshold(df$pval, 0.05)
  print("BH criticla value at 0.05")
  
  
  #add back in SNP names
  snps <- colnames(rbeta_data)
  df <- cbind(snps, df)
  
  #write the fully processed data into a .csv, this format is ready for the Manhattan plot script
  write.csv(df, paste("../../outputs/no_koot/final_lfmm_K7_",df[1,3],".csv", sep = ''), row.names=FALSE)
  
  #write the candidate snp files using different p value thresholds
  cand_SNP_K9 <- subset(df, pval <= 0.05)
  write.csv(cand_SNP_K9, paste("../../outputs/no_koot/final_lfmm_K7_05_",df[1,3],".csv", sep = ''), row.names=FALSE)
  cand_SNP_K9 <- subset(df, pval <= bonferroni)
  write.csv(cand_SNP_K9, paste("../../outputs/no_koot/final_lfmm_K7_Bon_",df[1,3],".csv", sep = ''), row.names=FALSE)
  cand_SNP_K9 <- subset(df, pval <= criticalValue)
  write.csv(cand_SNP_K9, paste("../../outputs/no_koot/final_lfmm_K7_BH_",df[1,3],".csv", sep = ''), row.names=FALSE)
  
  print('output cand SNPs with multiple p-value thresholds')
  
  #for plotting, need to subset the strangely named chromosomes, format, and then recombine
  sub2 <- df %>% filter(str_count(df$snps, "_") == 2)
  sub3 <- df %>% filter(str_count(df$snps, "_") == 3)
  
  
  #separate the snp column to get chr and bp
  newCols1 <- colsplit(sub2$snps, "_", c("chr_name","BP", "junk"))
  newCols1$junk <- NULL
  newCols2 <- colsplit(newCols1$chr_name, "y", c("junk", "chr"))
  newCols2$junk <- NULL
  df_pval_sub2 <- cbind(sub2, newCols1, newCols2)
  
  #same as above but for those with three _ in the name
  newCols1 <- colsplit(sub3$snps, "_", c("chr_name","junk", "BP", "junk1"))
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
  write.csv(df_ggplot, paste("../../outputs/no_koot/figures/final_lfmm_K7_plot_",df[1,3],".csv", sep = ''), row.names=FALSE)
  
  #Then we need to prepare the X axis. Want to display the cumulative position of SNP in bp, but just show the chromosome number
  axisdf = df_ggplot %>% 
    group_by(chr) %>% 
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  #create plot title
  plot_title <- paste("LFMM, K=9,", df_ggplot[1,3])
  
  #ggplot
  man_plot <- ggplot(df_ggplot, aes(x=BPcum, y=-log10(as.numeric(pval)))) +
    geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("darkgrey", "dodgerblue4"), 32 )) +
    #geom_hline(yintercept=-log10(bonferroni), color = "darkred") + #line for Bonferroni correction
    geom_hline(yintercept=-log10(criticalValue), color = "blue") + #line for BH correction
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    ggtitle(plot_title) +
    ylab("-log10 P Value") +
    theme_classic(base_size = 14, base_family = "Times") +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  ggsave(filename = paste("../../outputs/no_koot/figures/", plot_title,".png", sep = ""), plot = man_plot, height = 6, width = 10, units = "in")
  
  print(paste('Saved plot of ', plot_title), sep = '')
  
})

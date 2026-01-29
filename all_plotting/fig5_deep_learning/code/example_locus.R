library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(ggExtra)
library(patchwork)
library(paletteer)
library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(scales)
library(UpSetR)
library(forcats)

data_dir <- "/Users/omarahmed/Downloads/nygc/autoimmune/workspaces/workspace23_perspective/data/"
output_dir <- "/Users/omarahmed/Downloads/nygc/autoimmune/workspaces/workspace23_perspective/plots/"

##############################################################
# plot 1: GWAS p-value plot
##############################################################

df <- read.csv(paste(data_dir, "example_locus_1_data.csv", sep=""))

# filter variants below a certain p-value
print(paste("Number of rows in dataframe:", nrow(df)))
df_filt <- df %>% filter(hg19_pos >= 56600000 & hg19_pos <= 57200000)
print(paste("Number of rows in dataframe:", nrow(df_filt)))

plot_1 <- ggplot(df_filt, aes(x=hg19_pos, 
                              y=neglog10p, 
                              fill=gws,
                              color=gws, 
                              shape=hit_variant, 
                              size=hit_variant,
                              alpha=hit_variant)) +
          geom_point() +
          scale_color_manual(values=c("Yes"="#1f78b4", "No"="#66a61e")) +
          scale_fill_manual(values=c("Yes"="#1f78b4", "No"="#66a61e")) +
          scale_shape_manual(values=c("No"=16, "Yes"=18)) +
          scale_size_manual(values=c("No"=2, "Yes"=6)) +
          scale_alpha_manual(values=c("No"=0.3, "Yes"=1.0)) +
          theme_classic() +
          scale_y_continuous(limits=c(0, 11), breaks=seq(0, 10, 5)) +
          labs(x="Position (hg19)", y="-log10p", title="") +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                legend.position = "none")
plot_1

ggsave(paste(output_dir, "examplelocus-plot1-pvalue.pdf", sep=""),
       device="pdf",
       plot=plot_1,
       height=4,
       width=8,
       units=c("in"),
       dpi=600)

##############################################################
# plot 2: PIP plot
##############################################################

df <- read.csv(paste(data_dir, "example_locus_1_data.csv", sep=""))

# filter variants below a certain p-value
print(paste("Number of rows in dataframe:", nrow(df)))
df_filt <- df %>% filter(hg19_pos >= 56600000 & hg19_pos <= 57200000 & pip != -1)
print(paste("Number of rows in dataframe:", nrow(df_filt)))

plot_2 <- ggplot(df_filt, aes(x=hg19_pos, 
                              y=pip, 
                              fill=gws,
                              color=gws, 
                              shape=hit_variant, 
                              size=hit_variant,
                              alpha=hit_variant)) +
          geom_point() +
          scale_color_manual(values=c("Yes"="#1f78b4", "No"="#66a61e")) +
          scale_fill_manual(values=c("Yes"="#1f78b4", "No"="#66a61e")) +
          scale_shape_manual(values=c("No"=16, "Yes"=18)) +
          scale_size_manual(values=c("No"=2, "Yes"=6)) +
          scale_alpha_manual(values=c("No"=0.3, "Yes"=1.0)) +
          theme_classic() +
          scale_y_continuous(limits=c(0, 0.09), breaks=seq(0, 0.08, 0.02)) +
          labs(x="Position (hg19)", y="PIP", title="") +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                legend.position = "none")
plot_2

ggsave(paste(output_dir, "examplelocus-plot2-pips.pdf", sep=""),
       device="pdf",
       plot=plot_2,
       height=4,
       width=8,
       units=c("in"),
       dpi=600)


##############################################################
# plot 3: alphagenome score percentiles
##############################################################

df <- read.csv(paste(data_dir, "example_locus_1_data.csv", sep=""))

# filter variants below a certain p-value
print(paste("Number of rows in dataframe:", nrow(df)))
df_filt <- df %>% filter(gws == "Yes" & ag_score != -1) %>%
                  mutate(log_ag_score=log10(1+ag_score),
                         ag_score_percentile=percent_rank(ag_score)*100)
print(paste("Number of rows in dataframe:", nrow(df_filt)))




plot_3 <- ggplot(df_filt, aes(x=log_ag_score, 
                              y=ag_score_percentile, 
                              fill=gws,
                              color=gws, 
                              shape=hit_variant, 
                              size=hit_variant,
                              alpha=hit_variant)) +
          geom_point() +
          scale_color_manual(values=c("Yes"="#1f78b4", "No"="#66a61e")) +
          scale_fill_manual(values=c("Yes"="#1f78b4", "No"="#66a61e")) +
          scale_shape_manual(values=c("No"=16, "Yes"=18)) +
          scale_size_manual(values=c("No"=2, "Yes"=6)) +
          scale_alpha_manual(values=c("No"=0.3, "Yes"=1.0)) +
          theme_classic() +
          labs(x="log10(1+score)", y="Percentile", title="AlphaGenome") +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                legend.position = "none")
plot_3

ggsave(paste(output_dir, "examplelocus-plot3-ag_score_percentile.pdf", sep=""),
       device="pdf",
       plot=plot_3,
       height=3,
       width=5,
       units=c("in"),
       dpi=600)

##############################################################
# plot 4: borzoi prime score percentiles
##############################################################

df <- read.csv(paste(data_dir, "example_locus_1_data.csv", sep=""))

# filter variants below a certain p-value
print(paste("Number of rows in dataframe:", nrow(df)))
df_filt <- df %>% filter(gws == "Yes" & bp_score != -1) %>%
                  mutate(log_bp_score=log10(1+bp_score),
                         bp_score_percentile=percent_rank(bp_score)*100)
print(paste("Number of rows in dataframe:", nrow(df_filt)))




plot_4 <- ggplot(df_filt, aes(x=log_bp_score, 
                              y=bp_score_percentile, 
                              fill=gws,
                              color=gws, 
                              shape=hit_variant, 
                              size=hit_variant,
                              alpha=hit_variant)) +
          geom_point() +
          scale_color_manual(values=c("Yes"="#1f78b4", "No"="#66a61e")) +
          scale_fill_manual(values=c("Yes"="#1f78b4", "No"="#66a61e")) +
          scale_shape_manual(values=c("No"=16, "Yes"=18)) +
          scale_size_manual(values=c("No"=2, "Yes"=6)) +
          scale_alpha_manual(values=c("No"=0.3, "Yes"=1.0)) +
          theme_classic() +
          labs(x="log10(1+score)", y="Percentile", title="Borzoi Prime") +
          scale_x_continuous(breaks=seq(0, 1.5, 0.5), limits=c(0,1.5)) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                legend.position = "none")
plot_4

ggsave(paste(output_dir, "examplelocus-plot4-bp_score_percentile.pdf", sep=""),
       device="pdf",
       plot=plot_4,
       height=3,
       width=5,
       units=c("in"),
       dpi=600)


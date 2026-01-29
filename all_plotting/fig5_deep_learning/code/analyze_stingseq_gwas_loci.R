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
library(purrr)
library(broom)

data_dir <- "/Users/omarahmed/Downloads/nygc/autoimmune/workspaces/workspace23_perspective/data/"
output_dir <- "/Users/omarahmed/Downloads/nygc/autoimmune/workspaces/workspace23_perspective/plots/"

##############################################################
# plot 1: number of variants in each hit locus when filtering
#         variants based on p-value threshold of 5e-8
##############################################################

df <- read.csv(paste(data_dir, "variant_counts_per_locus.hit_loci.prior_to_finalizing.csv", sep=""))
df <- df %>% filter(pval_threshold == 5e-08) %>%
             mutate(variant_bin = case_when(num_variants == 0 ~ "0",
                                            num_variants >= 1 & num_variants <= 10 ~ "1–10",
                                            num_variants > 10 & num_variants <= 100 ~ "11–100",
                                            num_variants > 100 & num_variants <= 1000  ~ "101–1000",
                                            num_variants > 1000 ~ "1000+")
                   )
df$variant_bin <- factor(df$variant_bin, levels = c("0", "1–10", "11–100", "101–1000", "1000+"))

df_counts <- df %>% count(variant_bin, consort) %>% complete(variant_bin, consort, fill = list(n = 0))

consort_colors <- c("ukbb" = "#1f78b4", "bcc" = "#66a61e")

plot_1 <- ggplot(df_counts, aes(x=variant_bin, y=n, fill=consort)) +
          geom_col(position=position_dodge(width=0.9), color="black", linewidth=0.5) +
          geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.3, size=5) +
          theme_classic() +
          labs(x="Number of variants", y="Loci", title="p-value threhsold = 5e-8") +
          scale_fill_manual(values=consort_colors, name="Consortium", labels=c("bcc"="Blood Cell",
                                                                               "ukbb"="UKBB")) +
          scale_y_continuous(expand=c(0,0), limits=c(0, 65), breaks=seq(0, 60, 20)) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"))
plot_1

ggsave(paste(output_dir, "gwasloci-plot1-num_variants_per_hit_locus_5e8_threshold.pdf", sep=""),
       device="pdf",
       plot=plot_1,
       height=4,
       width=8,
       units=c("in"),
       dpi=600)

##############################################################
# plot 2: number of variants in each hit locus when filtering
#         variants based on p-value threshold of 1e-6
##############################################################

df <- read.csv(paste(data_dir, "variant_counts_per_locus.hit_loci.prior_to_finalizing.csv", sep=""))
df <- df %>% filter(pval_threshold == 1e-06) %>%
  mutate(variant_bin = case_when(num_variants == 0 ~ "0",
                                 num_variants >= 1 & num_variants <= 10 ~ "1–10",
                                 num_variants > 10 & num_variants <= 100 ~ "11–100",
                                 num_variants > 100 & num_variants <= 1000  ~ "101–1000",
                                 num_variants > 1000 ~ "1000+"))
df$variant_bin <- factor(df$variant_bin, levels = c("0", "1–10", "11–100", "101–1000", "1000+"))

df_counts <- df %>% count(variant_bin, consort) %>% complete(variant_bin, consort, fill = list(n = 0))
consort_colors <- c("ukbb" = "#1f78b4", "bcc" = "#66a61e")

plot_2 <- ggplot(df_counts, aes(x=variant_bin, y=n, fill=consort)) +
          geom_col(position=position_dodge(width=0.9), color="black", linewidth=0.5) +
          geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.3, size=5) +
          theme_classic() +
          labs(x="Number of variants", y="Loci", title="p-value threhsold = 1e-6") +
          scale_fill_manual(values=consort_colors, name="Consortium", labels=c("bcc"="Blood Cell",
                                                                               "ukbb"="UKBB")) +
          scale_y_continuous(expand=c(0,0), limits=c(0, 85), breaks=seq(0, 80, 20)) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"))
plot_2

ggsave(paste(output_dir, "gwasloci-plot2-num_variants_per_hit_locus_1e6_threshold.pdf", sep=""),
       device="pdf",
       plot=plot_2,
       height=4,
       width=8,
       units=c("in"),
       dpi=600)


##############################################################
# plot 3: percentiles for hit/non-hit variants
#         - uses density curves
##############################################################

df <- read.csv(paste(data_dir, "percentile_scores.csv", sep=""))

tool_labels <- c(
  "alphagenome" = "Alphagenome",
  "enformer"    = "Enformer",
  "borzoi"  = "Borzoi",
  "borzoi_prime" = "Borzoi Prime"
)

plot_3 <- ggplot(df, aes(x=vep_percentile, fill=locus_type, color=locus_type)) +
          geom_density(alpha=0.35, linewidth=1) +
          facet_wrap(~tool, labeller = labeller(tool = tool_labels)) +
          theme_classic() +
          labs(x="Score percentile (within locus)", y="Density") +
          scale_y_continuous(expand=c(0,0), limits=c(0, 0.021)) +
          scale_fill_discrete(name="Locus type",
                              labels=c("hit_loci"="Hit loci",
                                       "nonhit_loci"="Non-hit loci")) +
          scale_color_discrete(name="Locus type",
                               labels=c("hit_loci"="Hit loci",
                                        "nonhit_loci"="Non-hit loci")) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                strip.text = element_text(size=18, color="black"),
                panel.border = element_blank(),
                strip.background = element_rect(color = NA),
                legend.title=element_text(size=16, color="black"),
                legend.text=element_text(size=16, color="black"),
                panel.spacing.x=unit(1, "lines"))
plot_3

ggsave(paste(output_dir, "gwasloci-plot3-library_vep_percentiles.pdf", sep=""),
       device="pdf",
       plot=plot_3,
       height=5,
       width=8,
       units=c("in"),
       dpi=600)

##############################################################
# plot 4: percentiles for hit/non-hit variants
#         - uses box-plots curves
##############################################################
df <- read.csv(paste(data_dir, "percentile_scores.csv", sep=""))

tool_labels <- c(
  "alphagenome" = "Alphagenome",
  "enformer"    = "Enformer",
  "borzoi"  = "Borzoi",
  "borzoi_prime" = "Borzoi Prime"
)

plot_4 <- ggplot(df, aes(x=vep_percentile, fill=locus_type)) +
          geom_boxplot(width=0.6, outlier.shape=NA, linewidth=0.8) +
          facet_wrap(~tool, labeller = labeller(tool = tool_labels)) +
          theme_classic() +
          labs(x="Score percentile (within locus)") +
          scale_fill_discrete(name="Locus type",
                              labels=c("hit_loci"="Hit loci",
                                       "nonhit_loci"="Non-hit loci")) +
          scale_color_discrete(name="Locus type",
                               labels=c("hit_loci"="Hit loci",
                                        "nonhit_loci"="Non-hit loci")) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),

                # remove y-axis entirely
                axis.title.y = element_blank(),
                axis.text.y  = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y  = element_blank(),
                
                strip.text = element_text(size=18, color="black"),
                panel.border = element_blank(),
                strip.background = element_rect(color = NA),
                legend.title=element_text(size=16, color="black"),
                legend.text=element_text(size=16, color="black"),
                panel.spacing.x=unit(1, "lines"))
plot_4

ggsave(paste(output_dir, "gwasloci-plot4-library_vep_percentiles.pdf", sep=""),
       device="pdf",
       plot=plot_4,
       height=5,
       width=8,
       units=c("in"),
       dpi=600)

#################################################
# analysis: see if difference in percentiles
# is significant or not
#################################################
df <- df %>% mutate(locus_type=factor(locus_type, levels=c("nonhit_loci", "hit_loci")))

wilcox_results <- df %>% group_by(tool) %>%
                         summarise(n_hit=sum(locus_type == "hit_loci"),
                                   n_nonhit=sum(locus_type == "nonhit_loci"),
                                   median_hit=median(vep_percentile[locus_type == "hit_loci"]),
                                   median_nonhit=median(vep_percentile[locus_type == "nonhit_loci"]),
                                   test = list(wilcox.test(vep_percentile[locus_type == "hit_loci"],
                                                           vep_percentile[locus_type == "nonhit_loci"],
                                                           alternative = "greater",
                                                           exact = FALSE)
                                              ),
                                   .groups = "drop"
                                   ) %>%
                          mutate(p_value=map_dbl(test, ~ .x$p.value)) %>%
                          select(-test)
wilcox_results


##############################################################
# plot 5: number of variants in loci used for deep-learning
#         scoring (both hits and non-hits)
##############################################################

df <- read.csv(paste(data_dir, "loci_size.csv", sep=""))

locus_labels <- c(
  "hit_loci" = "Hit loci",
  "nonhit_loci"    = "Non-hit loci"
)

plot_5 <- ggplot(df, aes(x=num_variants)) +
          geom_histogram(bins=40, fill="#1f78b4", color="black")+
          facet_wrap(~locus_type, labeller=labeller(locus_type=locus_labels)) +
          theme_classic() +
          labs(x="Number of variants", y="") +
          scale_y_continuous(expand=c(0,0), limits=c(0, 13), breaks=seq(0,12,4)) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                strip.text = element_text(size=18, color="black"),
                panel.border = element_blank(),
                strip.background = element_rect(color = NA),
                legend.title=element_text(size=16, color="black"),
                legend.text=element_text(size=16, color="black"),
                panel.spacing.x=unit(1, "lines"))
plot_5

ggsave(paste(output_dir, "gwasloci-plot5-num_variants_per_scored_locus.pdf", sep=""),
       device="pdf",
       plot=plot_5,
       height=3,
       width=7,
       units=c("in"),
       dpi=600)


##############################################################
# plot 6 - 11: looking at how the different tools correlate
#               - plot 6: alphagenome vs borzoi
#               - plot 7: alphagenome vs enformer
#               - plot 8: enformer vs borzoi
#               - plot 9: alphagenome vs borzoi_prime
#               - plot 10: borzoi vs borzoi_prime
#               - plot 11: enformer vs borzoi_prime 
##############################################################

################
# plot 6
################
df <- read.csv(paste(data_dir, "percentile_scores.csv", sep=""))

# group together same varian
df_wide <- df %>% filter(tool %in% c("alphagenome", "borzoi")) %>%
                             pivot_wider(id_cols = c(chrom, pos), names_from=tool, values_from=vep_percentile)

# compute spearman correlation
r <- cor(df_wide$alphagenome, df_wide$borzoi, method = "spearman")

# round for label
r_label <- paste0("r = ", round(r, 2))

plot_6 <- ggplot(df_wide, aes(x=alphagenome, y=borzoi)) +
          geom_point(size=3, alpha=0.7, color="#1f78b4") +
          geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
          theme_classic() +
          labs(x="Alphagenome percentile", y="Borzoi percentile") +
          coord_equal() +
          scale_x_continuous(limits=c(0, 105), expand=c(0, 0)) +
          scale_y_continuous(limits=c(0, 105), expand=c(0, 0)) +
          annotate("text", x=5, y=97, label=r_label, hjust=0, vjust=0, size=6) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                strip.text = element_text(size=18, color="black"),
                panel.border = element_blank(),
                strip.background = element_rect(color = NA),
                legend.title=element_text(size=16, color="black"),
                legend.text=element_text(size=16, color="black"),
                panel.spacing.x=unit(1, "lines")
  )
plot_6

ggsave(paste(output_dir, "gwasloci-plot6-alphagenome_vs_borzoi.pdf", sep=""),
       device="pdf",
       plot=plot_6,
       height=4,
       width=4,
       units=c("in"),
       dpi=600)

################
# plot 7
################
df <- read.csv(paste(data_dir, "percentile_scores.csv", sep=""))

# group together same varian
df_wide <- df %>% filter(tool %in% c("alphagenome", "enformer")) %>%
                             pivot_wider(id_cols = c(chrom, pos), names_from=tool, values_from=vep_percentile)

# compute spearman correlation
r <- cor(df_wide$alphagenome, df_wide$enformer, method = "spearman")

# round for label
r_label <- paste0("r = ", round(r, 2))

plot_7 <- ggplot(df_wide, aes(x=alphagenome, y=enformer)) +
          geom_point(size=3, alpha=0.7, color="#1f78b4") +
          geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
          theme_classic() +
          labs(x="Alphagenome percentile", y="Enformer percentile") +
          coord_equal() +
          scale_x_continuous(limits=c(0, 105), expand=c(0, 0)) +
          scale_y_continuous(limits=c(0, 105), expand=c(0, 0)) +
          annotate("text", x=5, y=97, label=r_label, hjust=0, vjust=0, size=6) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                strip.text = element_text(size=18, color="black"),
                panel.border = element_blank(),
                strip.background = element_rect(color = NA),
                legend.title=element_text(size=16, color="black"),
                legend.text=element_text(size=16, color="black"),
                panel.spacing.x=unit(1, "lines")
          )
plot_7

ggsave(paste(output_dir, "gwasloci-plot7-alphagenome_vs_enformer.pdf", sep=""),
       device="pdf",
       plot=plot_7,
       height=4,
       width=4,
       units=c("in"),
       dpi=600)

################
# plot 8
################
df <- read.csv(paste(data_dir, "percentile_scores.csv", sep=""))

# group together same varian
df_wide <- df %>% filter(tool %in% c("borzoi", "enformer")) %>%
                         pivot_wider(id_cols = c(chrom, pos), names_from=tool, values_from=vep_percentile)

# compute spearman correlation
r <- cor(df_wide$borzoi, df_wide$enformer, method = "spearman")

# round for label
r_label <- paste0("r = ", round(r, 2))

plot_8 <- ggplot(df_wide, aes(x=borzoi, y=enformer)) +
          geom_point(size=3, alpha=0.7, color="#1f78b4") +
          geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
          theme_classic() +
          labs(x="Borzoi percentile", y="Enformer percentile") +
          coord_equal() +
          scale_x_continuous(limits=c(0, 105), expand=c(0, 0)) +
          scale_y_continuous(limits=c(0, 105), expand=c(0, 0)) +
          annotate("text", x=5, y=97, label=r_label, hjust=0, vjust=0, size=6) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                strip.text = element_text(size=18, color="black"),
                panel.border = element_blank(),
                strip.background = element_rect(color = NA),
                legend.title=element_text(size=16, color="black"),
                legend.text=element_text(size=16, color="black"),
                panel.spacing.x=unit(1, "lines")
          )
plot_8

ggsave(paste(output_dir, "gwasloci-plot8-enformer_vs_borzoi.pdf", sep=""),
       device="pdf",
       plot=plot_8,
       height=4,
       width=4,
       units=c("in"),
       dpi=600)

################
# plot 9
################
df <- read.csv(paste(data_dir, "percentile_scores.csv", sep=""))

# group together same varian
df_wide <- df %>% filter(tool %in% c("alphagenome", "borzoi_prime")) %>%
  pivot_wider(id_cols = c(chrom, pos), names_from=tool, values_from=vep_percentile)

# compute spearman correlation
r <- cor(df_wide$alphagenome, df_wide$borzoi_prime, method = "spearman")

# round for label
r_label <- paste0("r = ", round(r, 2))

plot_9 <- ggplot(df_wide, aes(x=alphagenome, y=borzoi_prime)) +
          geom_point(size=3, alpha=0.7, color="#1f78b4") +
          geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
          theme_classic() +
          labs(x="Alphagenome percentile", y="Borzoi Prime percentile") +
          coord_equal() +
          scale_x_continuous(limits=c(0, 105), expand=c(0, 0)) +
          scale_y_continuous(limits=c(0, 105), expand=c(0, 0)) +
          annotate("text", x=5, y=97, label=r_label, hjust=0, vjust=0, size=6) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                strip.text = element_text(size=18, color="black"),
                panel.border = element_blank(),
                strip.background = element_rect(color = NA),
                legend.title=element_text(size=16, color="black"),
                legend.text=element_text(size=16, color="black"),
                panel.spacing.x=unit(1, "lines")
          )
        plot_9

ggsave(paste(output_dir, "gwasloci-plot9-alphagenome_vs_borzoi_prime.pdf", sep=""),
       device="pdf",
       plot=plot_9,
       height=4,
       width=4,
       units=c("in"),
       dpi=600)

################
# plot 10
################
df <- read.csv(paste(data_dir, "percentile_scores.csv", sep=""))

# group together same varian
df_wide <- df %>% filter(tool %in% c("borzoi", "borzoi_prime")) %>%
  pivot_wider(id_cols = c(chrom, pos), names_from=tool, values_from=vep_percentile)

# compute spearman correlation
r <- cor(df_wide$borzoi, df_wide$borzoi_prime, method = "spearman")

# round for label
r_label <- paste0("r = ", round(r, 2))

plot_10 <- ggplot(df_wide, aes(x=borzoi, y=borzoi_prime)) +
          geom_point(size=3, alpha=0.7, color="#1f78b4") +
          geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
          theme_classic() +
          labs(x="Borzoi percentile", y="Borzoi Prime percentile") +
          coord_equal() +
          scale_x_continuous(limits=c(0, 105), expand=c(0, 0)) +
          scale_y_continuous(limits=c(0, 105), expand=c(0, 0)) +
          annotate("text", x=5, y=97, label=r_label, hjust=0, vjust=0, size=6) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                strip.text = element_text(size=18, color="black"),
                panel.border = element_blank(),
                strip.background = element_rect(color = NA),
                legend.title=element_text(size=16, color="black"),
                legend.text=element_text(size=16, color="black"),
                panel.spacing.x=unit(1, "lines")
          )
plot_10

ggsave(paste(output_dir, "gwasloci-plot10-borzoi_vs_borzoi_prime.pdf", sep=""),
       device="pdf",
       plot=plot_10,
       height=4,
       width=4,
       units=c("in"),
       dpi=600)

################
# plot 11
################
df <- read.csv(paste(data_dir, "percentile_scores.csv", sep=""))

# group together same varian
df_wide <- df %>% filter(tool %in% c("enformer", "borzoi_prime")) %>%
  pivot_wider(id_cols = c(chrom, pos), names_from=tool, values_from=vep_percentile)

# compute spearman correlation
r <- cor(df_wide$enformer, df_wide$borzoi_prime, method = "spearman")

# round for label
r_label <- paste0("r = ", round(r, 2))

plot_11 <- ggplot(df_wide, aes(x=enformer, y=borzoi_prime)) +
          geom_point(size=3, alpha=0.7, color="#1f78b4") +
          geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
          theme_classic() +
          labs(x="Enformer percentile", y="Borzoi Prime percentile") +
          coord_equal() +
          scale_x_continuous(limits=c(0, 105), expand=c(0, 0)) +
          scale_y_continuous(limits=c(0, 105), expand=c(0, 0)) +
          annotate("text", x=5, y=97, label=r_label, hjust=0, vjust=0, size=6) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                strip.text = element_text(size=18, color="black"),
                panel.border = element_blank(),
                strip.background = element_rect(color = NA),
                legend.title=element_text(size=16, color="black"),
                legend.text=element_text(size=16, color="black"),
                panel.spacing.x=unit(1, "lines")
          )
plot_11

ggsave(paste(output_dir, "gwasloci-plot11-enformer_vs_borzoi_prime.pdf", sep=""),
       device="pdf",
       plot=plot_11,
       height=4,
       width=4,
       units=c("in"),
       dpi=600)

##############################################################
# plot 12: percentiles for hit/non-hit variants
#         - uses density curves
#         - subset tools to relevant tools
##############################################################

df <- read.csv(paste(data_dir, "percentile_scores.csv", sep=""))

tool_labels <- c(
  "alphagenome" = "Alphagenome",
  "enformer"    = "Enformer",
  "borzoi"  = "Borzoi",
  "borzoi_prime" = "Borzoi Prime"
)

df <- df %>% filter(tool %in% c("alphagenome", "borzoi_prime", "enformer"))

plot_12 <- ggplot(df, aes(x=vep_percentile, fill=locus_type, color=locus_type)) +
          geom_density(alpha=0.35, linewidth=1) +
          facet_wrap(~tool, labeller=labeller(tool=tool_labels), nrow=1) +
          theme_classic() +
          labs(x="Score percentile (within locus)", y="Density") +
          scale_y_continuous(expand=c(0,0), limits=c(0, 0.027)) +
          scale_x_continuous(limits=c(-25, 150), breaks=seq(0, 100, 25)) +
          scale_fill_discrete(name="Locus type",
                              labels=c("hit_loci"="Hit loci",
                                       "nonhit_loci"="Non-hit loci")) +
          scale_color_discrete(name="Locus type",
                               labels=c("hit_loci"="Hit loci",
                                        "nonhit_loci"="Non-hit loci")) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                strip.text = element_text(size=18, color="black"),
                panel.border = element_blank(),
                strip.background = element_rect(color = NA),
                legend.title=element_text(size=16, color="black"),
                legend.text=element_text(size=16, color="black"),
                panel.spacing.x=unit(1, "lines"))
plot_12

ggsave(paste(output_dir, "gwasloci-plot12-library_vep_percentiles.pdf", sep=""),
       device="pdf",
       plot=plot_12,
       height=3,
       width=14,
       units=c("in"),
       dpi=600)


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
# plot 1: auprc values across tools on STINGseq dataset
##############################################################

df <- read.csv(paste(data_dir, "auprc_on_stingseq_datasets.csv", sep=""))
stingseq_df <- df %>% filter(author == "homebrew", dataset == "stingseq")

stingseq_df <- stingseq_df %>%
  mutate(
    tool_label = case_when(
      tool == "chrombpnet"      ~ "ChromBPNet",
      tool == "enformer"        ~ "Enformer",
      tool == "borzoi"          ~ "Borzoi",
      tool == "borzoi_prime"    ~ "Borzoi Prime",
      tool == "gpn_msa"         ~ "GPN-MSA",
      tool == "alphagenome"     ~ "AlphaGenome",
      TRUE ~ tool
    )
  )

tool_levels <- c(
  "ChromBPNet",
  "Enformer",
  "Borzoi",
  "Borzoi Prime",
  "GPN-MSA",
  "AlphaGenome"
)

tool_colors <- c(
  "ChromBPNet"    = "#1b9e77",
  "Enformer"      = "#d95f02",
  "Borzoi"        = "#7570b3",
  "Borzoi Prime"  = "#1f78b4",
  "GPN-MSA"       = "#e7298a",
  "AlphaGenome"   = "#66a61e"
)

stingseq_df <- stingseq_df %>% mutate(tool_label = factor(tool_label, levels = tool_levels))

plot_1 <- ggplot(stingseq_df, aes(x=reorder(tool_label, auprc), y=auprc, fill=tool_label)) +
          geom_bar(stat="identity") +
          coord_flip() +
          theme_classic() +
          labs(x="Tool", y="AUPRC", title="STING-seq (Blood traits)") +
          scale_fill_manual(values=tool_colors) +
          geom_text(aes(label=sprintf("%.3f", auprc)), hjust = -0.15, size = 5) +
          scale_y_continuous(expand=c(0,0), limits=c(0,0.5), breaks=seq(0, 0.4, 0.2)) +
          geom_hline(yintercept=0.25, color="black", linetype="dashed", linewidth=1) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                legend.position = "none")
plot_1

ggsave(paste(output_dir, "stingseq-plot1-auprc.pdf", sep=""),
       device="pdf",
       plot=plot_1,
       height=4,
       width=6,
       units=c("in"),
       dpi=600)

##############################################################
# plot 2: normalized scores for STING-seq variants across
#         different tools
##############################################################

df <- read.csv(paste(data_dir, "stingseq_variant_scores_with_annotations.csv", sep=""))

tool_levels <- c(
  "borzoi_prime",
  "alphagenome",
  "borzoi",
  "enformer",
  "chrombpnet",
  "gpn_msa"
)

tool_labels <- c(
  "borzoi_prime" = "Borzoi Prime",
  "alphagenome"  = "AlphaGenome",
  "borzoi"       = "Borzoi",
  "enformer"     = "Enformer",
  "chrombpnet"   = "ChromBPNet",
  "gpn_msa"      = "GPN-MSA"
)


df <- df %>% mutate(tool_label=factor(tool, levels=tool_levels, 
                                            labels=tool_labels[tool_levels]))
df$class <- factor(df$class, levels = c("stingseq_hits", "stingseq_nonhits"),
                             labels = c("STING-seq hits (n=131)", "STING-seq non-hits (n=390)"))

plot_2 <- ggplot(df, aes(x=tool_label, y=normalized_score, fill=class)) +
          geom_boxplot(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = c("STING-seq hits (n=131)"     = "#66a61e",
                                       "STING-seq non-hits (n=390)" = "#1f78b4")) +
          theme_classic() +
          labs(x = "Tool", fill="Group", y="Normalized Score") +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                legend.text=element_text(size=12),
                legend.title=element_text(size=12))
plot_2

ggsave(paste(output_dir, "stingseq-plot2-variant_vep_scores_across_tools.pdf", sep=""),
       device="pdf",
       plot=plot_2,
       height=4,
       width=14,
       units=c("in"),
       dpi=600)

##############################################################
# plot 3: visualize the vep scores across tools when we 
#         stratify based on TSS distance
##############################################################

df <- read.csv(paste(data_dir, "stingseq_variant_scores_with_annotations.csv", sep=""))

tool_levels <- c(
  "borzoi_prime",
  "alphagenome",
  "borzoi",
  "enformer",
  "chrombpnet",
  "gpn_msa"
)

tool_labels <- c(
  "borzoi_prime" = "Borzoi Prime",
  "alphagenome"  = "AlphaGenome",
  "borzoi"       = "Borzoi",
  "enformer"     = "Enformer",
  "chrombpnet"   = "ChromBPNet",
  "gpn_msa"      = "GPN-MSA"
)

df <- df %>% filter(class=="stingseq_hits")

df <- df %>% mutate(tool_label=factor(tool, levels=tool_levels, 
                                      labels=tool_labels[tool_levels]))

df$dist_group <- factor(df$dist_group, levels=c(1, 2, 3),
                                       labels=c("[0, 10k) (n=69)", "[10k, 50k) (n=36)", "(50k, 500k] (n=26)"))

plot_3 <- ggplot(df, aes(x=tool_label, y=normalized_score, fill=dist_group)) +
          geom_boxplot(position=position_dodge(width=0.75)) +
          theme_classic() +
          labs(x = "Tool", fill="Dist. to TSS (bp)", y="Normalized Score") +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                legend.text=element_text(size=12),
                legend.title=element_text(size=12))
plot_3

ggsave(paste(output_dir, "stingseq-plot3-variant_vep_scores_across_distances.pdf", sep=""),
       device="pdf",
       plot=plot_3,
       height=4,
       width=14,
       units=c("in"),
       dpi=600)

##############################################################
# plot 4: visualize the vep scores across tools when we 
#         stratify based on functional marks
##############################################################

df <- read.csv(paste(data_dir, "stingseq_variant_scores_with_annotations.csv", sep=""))

tool_levels <- c(
  "borzoi_prime",
  "alphagenome",
  "borzoi",
  "enformer",
  "chrombpnet",
  "gpn_msa"
)

tool_labels <- c(
  "borzoi_prime" = "Borzoi Prime",
  "alphagenome"  = "AlphaGenome",
  "borzoi"       = "Borzoi",
  "enformer"     = "Enformer",
  "chrombpnet"   = "ChromBPNet",
  "gpn_msa"      = "GPN-MSA"
)

df <- df %>% filter(class=="stingseq_hits")

df <- df %>% mutate(tool_label=factor(tool, levels=tool_levels, 
                                      labels=tool_labels[tool_levels]))

df$func_group <- factor(df$func_group, levels=c(1, 2, 3),
                                       labels=c("3 (n=34)", "2 (n=43)", "1 (n=54)"))

plot_4 <- ggplot(df, aes(x=tool_label, y=normalized_score, fill=func_group)) +
          geom_boxplot(position=position_dodge(width=0.75)) +
          theme_classic() +
          labs(x = "Tool", fill="Functional Marks", y="Normalized Score") +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                legend.text=element_text(size=12),
                legend.title=element_text(size=12))
plot_4

ggsave(paste(output_dir, "stingseq-plot4-variant_vep_scores_across_functional_marks.pdf", sep=""),
       device="pdf",
       plot=plot_4,
       height=4,
       width=14,
       units=c("in"),
       dpi=600)

##############################################################
# plot 5: precision-recall curves for each tool trying to 
#         classify between hits and non-hits
##############################################################

df <- read.csv(paste(data_dir, "precision_recall_curve_on_stingseq_datasets.csv", sep=""))

tool_labels <- c("alphagenome"="Alphagenome",
                 "borzoi"="Borzoi",
                 "borzoi_prime"="Borzoi Prime",
                 "chrombpnet"="ChromBPNet",
                 "enformer"="Enformer",
                 "gpn_msa"="GPN-MSA")

plot_5 <- ggplot(df, aes(x=recall, y=precision, color=tool, group=tool)) +
          geom_line(size=1) +
          theme_classic() +
          labs(x="Recall", y="Precision", color="Tool") +
          scale_color_discrete(labels=tool_labels) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                legend.text=element_text(size=12),
                legend.title=element_text(size=12))
plot_5

ggsave(paste(output_dir, "stingseq-plot5-precision_recall_curves.pdf", sep=""),
       device="pdf",
       plot=plot_5,
       height=4,
       width=6,
       units=c("in"),
       dpi=600)

##############################################################
# plot 6: same as plot 4, but now I want to include a box
#         for non-hits as well
##############################################################

df <- read.csv(paste(data_dir, "stingseq_variant_scores_with_annotations.csv", sep=""))

tool_levels <- c(
  "borzoi_prime",
  "alphagenome",
  "borzoi",
  "enformer",
  "chrombpnet",
  "gpn_msa"
)

tool_labels <- c(
  "borzoi_prime" = "Borzoi Prime",
  "alphagenome"  = "AlphaGenome",
  "borzoi"       = "Borzoi",
  "enformer"     = "Enformer",
  "chrombpnet"   = "ChromBPNet",
  "gpn_msa"      = "GPN-MSA"
)


# add labels for tools
df <- df %>% mutate(tool_label=factor(tool, levels=tool_levels, 
                                      labels=tool_labels[tool_levels]))

# factorize and create levels for functional groups
df$func_group <- factor(df$func_group, levels=c(1, 2, 3),
                        labels=c("3 (n=34)", "2 (n=43)", "1 (n=54)"))

# add labels for hits vs non-hits
df$class <- factor(df$class,
                   levels = c("stingseq_hits", "stingseq_nonhits"),
                   labels = c("STING-seq hits", "STING-seq non-hits"))


# create new column in order to order the boxs in plot
df <- df %>% mutate(func_class=interaction(func_group, class, sep = " + "))

# create that ordering
df$func_class <- factor(df$func_class, levels=c("3 (n=34) + STING-seq hits",
                                                "3 (n=34) + STING-seq non-hits",
                                                "2 (n=43) + STING-seq hits",
                                                "2 (n=43) + STING-seq non-hits",
                                                "1 (n=54) + STING-seq hits",
                                                "1 (n=54) + STING-seq non-hits"))

plot_6 <- ggplot(df, aes(x=tool_label, 
                         y=normalized_score, 
                         fill=func_group, 
                         color=class, 
                         group=interaction(tool_label, func_class))) +
          geom_boxplot(position=position_dodge2(width = 0.8, preserve = "single"),
                       alpha=1.0,
                       size=0.75) +
          scale_color_manual(values = c("STING-seq hits"="black",
                                        "STING-seq non-hits" = "grey70")) +
          theme_classic() +
          labs(x="Tool", 
               y="Normalized Score",
               fill="Functional Marks", 
               color="Variant group") +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                legend.text=element_text(size=12),
                legend.title=element_text(size=12))
plot_6

ggsave(paste(output_dir, "stingseq-plot6-variant_vep_scores_across_functional_marks_with_both_classes.pdf", sep=""),
       device="pdf",
       plot=plot_6,
       height=4,
       width=15,
       units=c("in"),
       dpi=600)

##############################################################
# plot 7: Same as plot 3, but focus on specific tools
##############################################################

df <- read.csv(paste(data_dir, "stingseq_variant_scores_with_annotations.csv", sep=""))

tool_levels <- c(
  "borzoi_prime",
  "alphagenome",
  "borzoi",
  "enformer",
  "chrombpnet",
  "gpn_msa"
)

tool_labels <- c(
  "borzoi_prime" = "Borzoi Prime",
  "alphagenome"  = "AlphaGenome",
  "borzoi"       = "Borzoi",
  "enformer"     = "Enformer",
  "chrombpnet"   = "ChromBPNet",
  "gpn_msa"      = "GPN-MSA"
)

df <- df %>% filter(class=="stingseq_hits")

df <- df %>% filter(tool %in% c("borzoi_prime", "alphagenome", "enformer"))

df <- df %>% mutate(tool_label=factor(tool, levels=tool_levels, 
                                      labels=tool_labels[tool_levels]))

df$dist_group <- factor(df$dist_group, levels=c(1, 2, 3),
                        labels=c("[0, 10k) (n=69)", "[10k, 50k) (n=36)", "(50k, 500k] (n=26)"))

plot_7 <- ggplot(df, aes(x=tool_label, y=normalized_score, fill=dist_group)) +
          geom_boxplot(position=position_dodge(width=0.75)) +
          theme_classic() +
          labs(x = "Tool", fill="Dist. to TSS (bp)", y="Normalized Score") +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                legend.text=element_text(size=12),
                legend.title=element_text(size=12))
plot_7

ggsave(paste(output_dir, "stingseq-plot7-variant_vep_scores_across_distances_subsetted.pdf", sep=""),
       device="pdf",
       plot=plot_7,
       height=4,
       width=9,
       units=c("in"),
       dpi=600)

##############################################################
# plot 8: same as plot 6, but subsetting to certain tools
##############################################################

df <- read.csv(paste(data_dir, "stingseq_variant_scores_with_annotations.csv", sep=""))

tool_levels <- c(
  "borzoi_prime",
  "alphagenome",
  "borzoi",
  "enformer",
  "chrombpnet",
  "gpn_msa"
)

tool_labels <- c(
  "borzoi_prime" = "Borzoi Prime",
  "alphagenome"  = "AlphaGenome",
  "borzoi"       = "Borzoi",
  "enformer"     = "Enformer",
  "chrombpnet"   = "ChromBPNet",
  "gpn_msa"      = "GPN-MSA"
)

# subset to certain tools
df <- df %>% filter(tool %in% c("borzoi_prime", "alphagenome", "enformer"))

# add labels for tools
df <- df %>% mutate(tool_label=factor(tool, levels=tool_levels, 
                                      labels=tool_labels[tool_levels]))

# factorize and create levels for functional groups
df$func_group <- factor(df$func_group, levels=c(1, 2, 3),
                        labels=c("3 (n=34)", "2 (n=43)", "1 (n=54)"))

# add labels for hits vs non-hits
df$class <- factor(df$class,
                   levels = c("stingseq_hits", "stingseq_nonhits"),
                   labels = c("STING-seq hits", "STING-seq non-hits"))


# create new column in order to order the boxs in plot
df <- df %>% mutate(func_class=interaction(func_group, class, sep = " + "))

# create that ordering
df$func_class <- factor(df$func_class, levels=c("3 (n=34) + STING-seq hits",
                                                "3 (n=34) + STING-seq non-hits",
                                                "2 (n=43) + STING-seq hits",
                                                "2 (n=43) + STING-seq non-hits",
                                                "1 (n=54) + STING-seq hits",
                                                "1 (n=54) + STING-seq non-hits"))

plot_8 <- ggplot(df, aes(x=tool_label, 
                                 y=normalized_score, 
                                 fill=func_group, 
                                 color=class, 
                                 group=interaction(tool_label, func_class))) +
          geom_boxplot(position=position_dodge2(width = 0.8, preserve = "single"),
                       alpha=1.0,
                       size=0.75) +
          scale_color_manual(values = c("STING-seq hits"="black",
                                        "STING-seq non-hits" = "grey70")) +
          theme_classic() +
          labs(x="Tool", 
               y="Normalized Score",
               fill="Functional Marks", 
               color="Variant group") +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                legend.text=element_text(size=12),
                legend.title=element_text(size=12))
plot_8

ggsave(paste(output_dir, "stingseq-plot8-variant_vep_scores_across_functional_marks_with_both_classes.pdf", sep=""),
       device="pdf",
       plot=plot_8,
       height=4,
       width=10,
       units=c("in"),
       dpi=600)

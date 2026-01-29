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


#########################################################
# plot 1, 2, 3: auprc values across tools on OMIM dataset
#               - plot 1: homebrew
#               - plot 2: traitgym
#               - plot 3: combined plot
#########################################################

#########
# plot 1
#########
df <- read.csv(paste(data_dir, "auprc_on_traitgym_datasets.csv", sep=""))
omim_df <- df %>% filter(author == "homebrew", dataset == "omim")

omim_df <- omim_df %>%
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

omim_df <- omim_df %>% mutate(tool_label = factor(tool_label, levels = tool_levels))

plot_1 <- ggplot(omim_df, aes(x=reorder(tool_label, auprc), y=auprc, fill=tool_label)) +
          geom_bar(stat="identity") +
          coord_flip() +
          theme_classic() +
          labs(x="Tool", y="AUPRC", title="OMIM (in-house)") +
          scale_fill_manual(values=tool_colors) +
          geom_text(aes(label=sprintf("%.3f", auprc)), hjust = -0.15, size = 5) +
          scale_y_continuous(expand=c(0,0), limits=c(0,0.85), breaks=seq(0, 0.8, 0.2)) +
          # geom_text(aes(label=count), hjust = -0.3, size=6, color = "black") +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                legend.position = "none")
plot_1

##########
# plot 2
##########
df <- read.csv(paste(data_dir, "auprc_on_traitgym_datasets.csv", sep=""))
omim_df <- df %>% filter(author == "traitgym", dataset == "omim")

omim_df <- omim_df %>%
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

omim_df <- omim_df %>% mutate(tool_label = factor(tool_label, levels = tool_levels))

plot_2 <- ggplot(omim_df, aes(x=reorder(tool_label, auprc), y=auprc, fill=tool_label)) +
          geom_bar(stat="identity") +
          coord_flip() +
          theme_classic() +
          labs(x="Tool", y="AUPRC", title="OMIM (TraitGym)") +
          scale_fill_manual(values=tool_colors) +
          geom_text(aes(label=sprintf("%.3f", auprc)), hjust = -0.15, size = 5) +
          scale_y_continuous(expand=c(0,0), limits=c(0,0.85), breaks=seq(0, 0.8, 0.2)) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                legend.position = "none")
plot_2

##########
# plot 3
##########
plot_3 <- plot_1 + plot_2
plot_3

ggsave(paste(output_dir, "traitgym-plot1-omim_auprc_comparison.pdf", sep=""),
       device="pdf",
       plot=plot_3,
       height=4,
       width=10,
       units=c("in"),
       dpi=600)

#########################################################
# plot 4, 5, 6: auprc values across tools on GWAS dataset
#               - plot 4: homebrew
#               - plot 5: traitgym
#               - plot 6: combined plot
#########################################################

#########
# plot 4
#########
df <- read.csv(paste(data_dir, "auprc_on_traitgym_datasets.csv", sep=""))
gwas_df <- df %>% filter(author == "homebrew", dataset == "gwas")

gwas_df <- gwas_df %>%
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

gwas_df <- gwas_df %>% mutate(tool_label = factor(tool_label, levels = tool_levels))

plot_4 <- ggplot(gwas_df, aes(x=reorder(tool_label, auprc), y=auprc, fill=tool_label)) +
          geom_bar(stat="identity") +
          coord_flip() +
          theme_classic() +
          labs(x="Tool", y="AUPRC", title="GWAS (in-house)") +
          scale_fill_manual(values=tool_colors) +
          geom_text(aes(label=sprintf("%.3f", auprc)), hjust = -0.15, size = 5) +
          scale_y_continuous(expand=c(0,0), limits=c(0,0.35), breaks=seq(0, 0.3, 0.1)) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                legend.position = "none")
plot_4

##########
# plot 5
##########
df <- read.csv(paste(data_dir, "auprc_on_traitgym_datasets.csv", sep=""))
gwas_df <- df %>% filter(author == "traitgym", dataset == "gwas")

gwas_df <- gwas_df %>%
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

gwas_df <- gwas_df %>% mutate(tool_label = factor(tool_label, levels = tool_levels))

plot_5 <- ggplot(gwas_df, aes(x=reorder(tool_label, auprc), y=auprc, fill=tool_label)) +
          geom_bar(stat="identity") +
          coord_flip() +
          theme_classic() +
          labs(x="Tool", y="AUPRC", title="GWAS (TraitGym)") +
          scale_fill_manual(values=tool_colors) +
          geom_text(aes(label=sprintf("%.3f", auprc)), hjust = -0.15, size = 5) +
          scale_y_continuous(expand=c(0,0), limits=c(0,0.35), breaks=seq(0, 0.3, 0.1)) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"),
                legend.position = "none")
plot_5

##########
# plot 6
##########
plot_6 <- plot_4 + plot_5
plot_6

ggsave(paste(output_dir, "traitgym-plot2-gwas_auprc_comparison.pdf", sep=""),
       device="pdf",
       plot=plot_6,
       height=4,
       width=10,
       units=c("in"),
       dpi=600)


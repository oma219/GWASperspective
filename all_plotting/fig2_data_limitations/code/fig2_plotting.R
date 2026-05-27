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
library(readxl)
options(bitmapType='cairo')

data_dir <- "/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/all_plotting/fig2_data_limitations/data/"
output_dir <- "/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/all_plotting/fig2_data_limitations/plots/"

######################################################
# plot 1: number of variants per GWAS (histogram)
######################################################

df <- read.csv(paste(data_dir, "snp_density_per_study.csv", sep=""), header=FALSE)
colnames(df) <- c("study", "snp_density")
df$snp_density_millions <- df$snp_density/1e6

plot_1 <- ggplot(df, aes(x=snp_density_millions)) +
          geom_histogram(fill="#348feb", color="black") +
          theme_classic() +
          labs(x="SNPs in GWAS summary statistic (millions)", 
               y="Count",
               title="Varying SNP density") +
          scale_y_continuous(expand=c(0,0), breaks=seq(0, 30, 5)) +
          scale_x_continuous(breaks=seq(0, 60, 10)) +
          theme(axis.line=element_line(linewidth=1),
                axis.ticks=element_line(linewidth=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=16),
                axis.text=element_text(size=16, color="black"),
                plot.title=element_text(hjust=0.5, size=16, face="bold"))
plot_1

ggsave(paste(output_dir, "fig2-snp_density_histogram.pdf", sep=""),
       plot=plot_1,
       height=4,
       width=6 ,
       units=c("in"),
       dpi=600)

######################################################
# plot 2: number of discrepancies in GWAS loci
######################################################

df <- read.csv(paste(data_dir, "gwas_data_discrepancies.csv", sep=""), header=FALSE)
colnames(df) <- c("diseaes","num_loci", "percent_loci")

df$group <- ifelse(df$num_loci == 0, 1,
                   ifelse(df$percent_loci <= 20, 2,
                          ifelse(df$percent_loci <= 40, 3,
                                 ifelse(df$percent_loci <= 60, 4, 
                                        ifelse(df$percent_loci <= 80, 5, 6)))))
group_counts <- df %>% group_by(group) %>% summarise(n = n()) %>% mutate(fraction = n / sum(n)) %>% arrange(desc(group)) %>%   mutate(ypos = cumsum(fraction) - 0.5 * fraction)
 

# Define colors for groups (customize here)
custom_colors <- c(
  "1" = "#348feb",
  "2" = "#ffb6c1",
  "3" = "#F68EA0",
  "4" = "#EE657F",
  "5" = "#E53D5D",
  "6" = "#dc143c"
)

plot_2 <-ggplot(group_counts, aes(x = 2, y = fraction, fill = as.factor(group))) +
          geom_bar(stat = "identity", width = 1, color = "black", size = 0.5) +
          coord_polar(theta = "y") +
          geom_text(aes(y = ypos, label = n), color = "white", size = 5) +
          xlim(0.5, 2.5) +
          theme_void() +
          theme(
            legend.position = "right",
            legend.title = element_blank(),
            legend.text = element_text(size = 12),
            legend.key = element_rect(fill = "transparent", color = NA),
            legend.key.size = unit(1.2, "lines")
          ) +
          scale_fill_manual(values = custom_colors) +
          guides(fill = guide_legend(
            override.aes = list(
              shape = 21,          # circle with fill and border
              size = 5,
              color = "black",     # border color (outline)
              stroke = 1.5         # border thickness
            )
          )) 
plot_2

ggsave(paste(output_dir, "fig2-discrepancies_doughnut.pdf", sep=""),
       plot=plot_2,
       height=4,
       width=5 ,
       units=c("in"),
       dpi=600)

##########################################################
# plot 3: in-sample data availability (bar graph)
##########################################################
all_categories <- c("LD", "FM", "No")

df <- read.csv(paste(data_dir, "in_sample_ld_available_v2.csv", sep=""), header=TRUE)
df$finemapping_available <- factor(df$finemapping_available,
                                   levels=all_categories)

counts <- as.data.frame(table(df$finemapping_available))
colnames(counts) <- c("finemapping_available", "count")


custom_colors <- c("LD"="black",
                   "FM"="#dc143c", 
                   "No"="#348feb")

plot_3 <- ggplot(counts, aes(x=finemapping_available, y=count, fill=finemapping_available)) +
          geom_bar(stat="identity") +
          geom_text(aes(label=count), hjust=-0.5, size=7) +
          scale_fill_manual(values=custom_colors) +
          coord_flip() +
          theme_classic() +
          labs(x="Data availability",
               y="Number of GWAS") +
          scale_y_continuous(expand=c(0,0), breaks=seq(0, 35, 5), limits=c(0, 36)) +
          theme(axis.line=element_line(linewidth=1),
                axis.ticks=element_line(linewidth=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                plot.title=element_text(hjust=0.5, size=16, face="bold"),
                plot.margin=margin(0.5,1,0.5,0.5,"cm"),
                legend.position="none")

plot_3


ggsave(paste(output_dir, "fig2-data_availability.pdf", sep=""),
       plot=plot_3,
       height=5,
       width=6 ,
       units=c("in"),
       dpi=600)


##########################################################
# plot 4: number of variants per GWAS vs year
#         - colored by genome build
##########################################################

# merge the dataframes with genome build and snp density
df <- read_xlsx(paste(data_dir, "builds_and_tech.xlsx", sep=""), col_names = FALSE)
colnames(df)<- c("name", "build", "tech", "numtech", "study")

snp_density <- read.csv(paste(data_dir, "snp_density_per_study.csv", sep=""), header=FALSE)
colnames(snp_density) <- c("study", "snp_density")

match_idx <- sapply(df[[5]], function(x) {
  hits <- which(str_detect(snp_density[[1]], fixed(x)) | str_detect(x, fixed(snp_density[[1]])))
  if (length(hits) == 0) return(NA_integer_)
  return(hits[1]) 
})

df$snp_density_col2 <- snp_density[[2]][match_idx]

df <- df %>% mutate(year=str_extract(as.character(study), "^[0-9]+"),
                    year=as.numeric(year),
                    snp_density_numeric=suppressWarnings(as.numeric(snp_density_col2)),
                    snp_density_millions = snp_density_numeric / 1e6)

year_min <- 2008
year_max <- 2024

year_breaks <- seq(year_min, year_max, by=2)

# plot the figure
plot_4 <- ggplot(df, aes(x=snp_density_millions, y=year, color=as.factor(build))) +
          geom_point(size = 3.2, alpha = 0.6) +   
          theme_classic() +
          labs(x="SNPs in GWAS summary statistic (millions)",
               y="Year",
               color="Build") +
          scale_x_continuous(expand=expansion(mult=c(0.02, 0.08)), 
                             breaks = scales::pretty_breaks(n = 6)) +
          scale_y_continuous(limits=c(year_min, year_max),
                             breaks=year_breaks,
                             labels=year_breaks) +
          scale_color_brewer(palette = "Set1") +
          theme(axis.line = element_line(size = 1),
                axis.ticks = element_line(size = 1),
                axis.ticks.length = unit(0.25, "cm"),
                axis.title = element_text(size = 18),
                axis.text = element_text(size = 18, color = "black"),
                legend.title = element_text(size = 18),
                legend.text = element_text(size = 18))

plot_4


ggsave(filename=paste(output_dir, "fig2-snp_density_histogram_vs_year_with_genome_build.pdf", sep=""),,
       plot=plot_4,
       width=8,
       height=6,
       units="in",
       dpi = 600)

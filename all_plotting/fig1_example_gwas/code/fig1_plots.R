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

output_dir <- "/Users/omarahmed/Downloads/"

##############################################################
# plot 1: PIP threshold vs Optimistic Signal-to-Noise
#         - bar plot
##############################################################

# version 1: all data
df <- data.frame(
        min_pip <- c(99, 90, 80, 50, 10, 1),
        snr <- c(0.995, 0.9815, 0.9511538462, 0.7762903226, 0.3480944625, 0.0846123298)
        )
# version 2: excluding 1 percent
df <- data.frame(
  min_pip <- c(99, 90, 80, 50, 10),
  snr <- c(0.995, 0.9815, 0.9511538462, 0.7762903226, 0.3480944625)
)
colnames(df) <- c("min_pip", "snr")

df$min_pip <- factor(df$min_pip, level=c(1, 10, 50, 80, 90, 99))

plot_1 <- ggplot(df, aes(x=min_pip, y=snr)) +
          geom_bar(stat="identity", fill="skyblue") +
          theme_classic() +
          labs(x="Minimum PIP threshold (%)", y="Signal-to-noise ratio") +
          scale_y_continuous(expand=c(0,0), limits=c(0, 1.05), breaks=seq(0, 1, 0.25)) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"))
plot_1         

ggsave(paste(output_dir, "fig1-plot1_signal_to_noise.pdf", sep=""),
       device="pdf",
       plot=plot_1,
       height=4,
       width=7,
       units=c("in"),
       dpi=600)


##############################################################
# plot 2: Number of CRISPR gRNAs across thresholds
##############################################################

# version 1: all data
df <- data.frame(
  min_pip <- c(99, 90, 80, 50, 10, 1),
  num_grnas <- c(28, 40, 52, 124, 1228, 10576)
)

# version 2: excluding 1 percent
df <- data.frame(
  min_pip <- c(99, 90, 80, 50, 10),
  num_grnas <- c(484, 608, 732, 1024, 4156)
)

colnames(df) <- c("min_pip", "num_grnas")

df$min_pip <- factor(df$min_pip, level=c(1, 10, 50, 80, 90, 99))

plot_2 <- ggplot(df, aes(x=min_pip, y=num_grnas)) +
          geom_bar(stat="identity", fill="skyblue") +
          geom_text(aes(label=num_grnas), vjust=-0.5, size=5) +
          theme_classic() +
          labs(x="Minimum PIP threshold (%)", 
               y="Library size",
               title="CRISPR (4 gRNAS/variant)") +
          scale_y_continuous(expand=c(0,0), labels=scales::comma, limits=c(0, 4250)) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=18, face="bold"))
plot_2         

ggsave(paste(output_dir, "fig1-plot2_crispr_grna_library.pdf", sep=""),
       device="pdf",
       plot=plot_2,
       height=4,
       width=8,
       units=c("in"),
       dpi=600)


##############################################################
# plot 3: Number of MPRA oligos across thresholds
##############################################################
# version 1: all data
df <- data.frame(
  min_pip <- c(99, 90, 80, 50, 10, 1),
  num_grnas <- c(700, 1000, 1300, 3100, 30700, 264400)
)

# version 2: excluding 1 percent
df <- data.frame(
  min_pip <- c(99, 90, 80, 50, 10),
  num_grnas <- c(6050, 7600, 9150, 12800, 51950)
)

colnames(df) <- c("min_pip", "num_grnas")

df$min_pip <- factor(df$min_pip, level=c(1, 10, 50, 80, 90, 99))
#df$num_grnas <- df$num_grnas/1000

plot_3 <- ggplot(df, aes(x=min_pip, y=num_grnas)) +
          geom_bar(stat="identity", fill="skyblue") +
          geom_text(aes(label=num_grnas), vjust=-0.5, size=5) +
          theme_classic() +
          labs(x="Minimum PIP threshold (%)", 
               y="Library size",
               title="MPRA (50 barcodes/allele)") +
          scale_y_continuous(expand=c(0,0), labels=scales::comma, limits=c(0, 55000)) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=18, face="bold"))
plot_3         

ggsave(paste(output_dir, "fig1-plot3_mpra_library.pdf", sep=""),
       device="pdf",
       plot=plot_3,
       height=4,
       width=8,
       units=c("in"),
       dpi=600)

##############################################################
# plot 4: Library size vs Number of cells needed
##############################################################
# version 1: all data
df <- data.frame(
  library_size <- c(28, 40, 52, 124, 1228, 10576, 28, 40, 52, 124, 1228, 10576, 28, 40, 52, 124, 1228, 10576),
  moi <- c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5),
  num_cells <- c(14000, 20000, 26000, 62000, 614000, 5288000, 2800, 4000, 5200, 12400, 122800, 1057600, 560, 800, 1040, 2480, 24560, 211520)
)

# version 2: exclusing 1 percent
df <- data.frame(
  library_size <- c(484, 608, 732, 1024, 4156, 484, 608, 732, 1024, 4156, 484, 608, 732, 1024, 4156),
  moi <- c(0.2, 0.2, 0.2, 0.2, 0.2,  1, 1, 1, 1, 1, 5, 5, 5, 5, 5),
  num_cells <- c(242000,304000,366000,512000,2078000, 48400,60800, 73200,102400,415600, 9680,12160,14640,20480,83120)
)

colnames(df) <- c("library_size", "moi", "num_cells")
df$moi <- factor(df$moi)

plot_4 <- ggplot(df, aes(x=library_size, y=num_cells, group=moi)) +
          geom_point(aes(color=moi), size=2) +
          geom_line(aes(color=moi))+
          theme_classic() +
          labs(x="Library size", 
               y="Number of cells needed") +
          #scale_y_continuous(expand=c(0,0), limits=c(0,325), breaks=seq(0, 300, 100)) +
          scale_color_manual(values=c("0.2"="red",
                                        "1"="#f2fa11",
                                        "5"="#1db512"),
                             name="MOI") +
          #scale_y_log10(labels = scales::comma) +
          scale_x_log10(breaks = c(500, 1000, 2500, 5000), 
                labels = scales::comma,
                limits = c(500, 6000)) +
          scale_y_log10(breaks = c(10000, 100000, 500000, 2000000), 
                labels = scales::comma,
                limits = c(5000, 2100000)) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=18, face="bold"),
                legend.title=element_text(size=18, color="black"),
                legend.text=element_text(size=18, color="black"))
plot_4   

ggsave(paste(output_dir, "fig1-plot4_num_cells_needed.pdf", sep=""),
       device="pdf",
       plot=plot_4,
       height=4,
       width=6,
       units=c("in"),
       dpi=600)

##############################################################
# plot 5: PIP threshold vs Optimistic Signal-to-Noise
#         - line plot
##############################################################

df <- data.frame(
  min_pip <- c(1, 10, 50, 80, 90, 100),
  num_variants <- c(12459, 1039, 256, 193, 152, 121)
)
colnames(df) <- c("min_pip", "num_variants")
df$min_pip <- as.integer(df$min_pip)

plot_5 <- ggplot(df, aes(x=min_pip, y=num_variants)) +
          geom_point(size=3) +
          geom_line(linewidth=0.5) +
          theme_classic() +
          labs(x="Minimum PIP threshold (%)", y="Number of variants") +
          #scale_y_continuous(expand=c(0,0), limits=c(0, 1.05), breaks=seq(0, 1, 0.25)) +
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1, color="black"),
                axis.ticks.length=unit(0.25, "cm"),
                axis.title=element_text(size=18),
                axis.text=element_text(size=18, color="black"),
                axis.text.x=element_text(size=18),
                plot.title=element_text(hjust=0.5, size=20, face="bold"))
plot_5         

ggsave(paste(output_dir, "fig1-plot5_num_variants_above_thresholds.pdf", sep=""),
       device="pdf",
       plot=plot_5,
       height=5,
       width=7,
       units=c("in"),
       dpi=600)

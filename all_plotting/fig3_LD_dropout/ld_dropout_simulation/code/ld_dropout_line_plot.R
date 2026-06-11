library(data.table)
library(stringr)
library(ggplot2)
library(grid)

in_fn <- "/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/all_plotting/fig3_LD_dropout/LD_dropout_simulation/data/merged_0.03_5e+05_2.csv"
dt <- fread(in_fn)
dt[, `1` := as.character(`1`)]
dt[, replicate := as.character(replicate)]
dt[, dropout_replicate := as.character(dropout_replicate)]
num_pat <- "^[-+]?(?:\\d*\\.?\\d+)(?:[eE][-+]?\\d+)?$"
dt <- dt[
  grepl(num_pat, str_trim(`1`)) &
    grepl(num_pat, str_trim(replicate)) &
    grepl(num_pat, str_trim(dropout_replicate))
]
dt[, `1` := as.numeric(`1`)]
dt[, replicate := as.integer(as.numeric(replicate))]
dt[, dropout_replicate := as.integer(as.numeric(dropout_replicate))]
dt[, causal := as.character(causal)]
dt[, PIP := as.character(PIP)]
dt[, PIP := str_trim(PIP)]
dt[, PIP := str_replace_all(PIP, ",", "")]
dt[PIP == "" | PIP == "NA" | PIP == "NaN", PIP := NA_character_]
dt[, PIP := as.numeric(PIP)]   # preserves decimals like 0.134, 1e-3, etc.
dropout_pct_col <- NULL
cands <- c("dropout_pct","dropout_percent","dropout_percentage","dropout","dropout_rate")

for (cc in cands) {
  if (cc %in% names(dt)) {
    dropout_pct_col <- cc
    break
  }
}
if (is.null(dropout_pct_col)) {
  dt[, dropout_pct := dropout_replicate]
  dropout_pct_col <- "dropout_pct"
} else {
  dt[, (dropout_pct_col) := as.character(get(dropout_pct_col))]
}
dt[, `1` := as.character(`1`)]
dt <- dt[grepl("^[-+]?[0-9]*\\.?[0-9]+$", `1`)]

dt[, `1` := as.numeric(`1`)]

cat("Rows remaining after numeric filter on column `1`:", nrow(dt), "\n")

tp_fp_df <- dt[, .(
  TP = sum(causal == "Y" & !is.na(PIP) & PIP > 0.1),
  FP = sum(causal == "N" & !is.na(PIP) & PIP > 0.1),
  n_rows = .N
), by = .(
  dropout_pct = get(dropout_pct_col),
  `1`,
  dropout_replicate,
  replicate
)]

tp_fp_df[order(dropout_pct, `1`, dropout_replicate, replicate)][1:25]

tp_fp_long <- melt(
  tp_fp_df,
  id.vars = c("dropout_pct", "1", "dropout_replicate", "replicate"),
  measure.vars = c("TP", "FP"),
  variable.name = "metric",
  value.name = "value"
)
tp_fp_long[, metric := factor(metric, levels = c("TP", "FP"))]
tp_fp_long[, dropout_pct := as.character(dropout_pct)]
tp_fp_long[, dropout_pct := factor(dropout_pct, levels = sort(unique(dropout_pct)))]
print(tp_fp_long[order(dropout_pct, `1`, dropout_replicate, replicate, metric)][1:30])
tp_fp_long[, dropout_pct := as.character(dropout_pct)]
tp_fp_long[, dropout_pct := factor(dropout_pct, levels = sort(unique(dropout_pct)))]

out_dir <- "/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/all_plotting/fig3_LD_dropout/LD_dropout_simulation/plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_pdf <- file.path(out_dir, "tp_fp_bar_sum_by_dropout_1.pdf")

tp_fp_sum <- as.data.table(tp_fp_long)[
  , .(value = sum(as.numeric(value), na.rm = TRUE)),
  by = .(dropout_pct, metric)
]
tp_fp_sum[, dropout_pct_chr := as.character(dropout_pct)]
tp_fp_sum[, dropout_pct_num := suppressWarnings(as.numeric(dropout_pct_chr))]
tp_fp_sum <- tp_fp_sum[!is.na(dropout_pct_num)]

setorder(tp_fp_sum, dropout_pct_num)
tp_fp_sum[, dropout_pct := factor(dropout_pct_chr, levels = unique(dropout_pct_chr))]

# colors
cols <- c("TP" = "#FBB040", "FP" = "#0073B2")

p <- ggplot(
  tp_fp_sum,
  aes(
    x = dropout_pct,
    y = value,
    color = metric,
    group = metric
  )
) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = cols, name = NULL) +
  labs(
    x = "Dropout (%)",
    y = "Total count",
    title = "Total TP vs FP by dropout percentage"
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(0, 2000, 500),
    limits = c(0, 2100)
  ) +
  geom_hline(
    yintercept = 600,
    color = "#FBB040",
    linetype = "dashed",
    size = 1
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 1, color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 18, color = "black"),
    plot.title   = element_blank(),
    legend.text  = element_text(size = 18),
    legend.title = element_text(size = 18)
  )
p


ggsave(out_pdf, 
       p, 
       width=12, 
       height=6, 
       units="in", 
       dpi=300)
cat("Saved plot to:", out_pdf, "\n")


#########################################################################
# Update by Omar on 2/3/25
# - Make the axis linear, it doesn't scale correctly
#########################################################################

data_dir <- "/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/all_plotting/fig3_LD_dropout/LD_dropout_simulation/data/"
output_dir <- "/gpfs/commons/groups/sanjana_lab/oahmed/gwas_perspective/all_plotting/fig3_LD_dropout/LD_dropout_simulation/plots/"

# Omar: just ran this once to save the processed data to a file
write.csv(tp_fp_sum, paste(data_dir, "tp_fp_sum.csv", sep=""), row.names = FALSE)


# Start here and load data ...
df <- read.csv(paste(data_dir, "tp_fp_sum.csv", sep=""))
df <- df %>% mutate(dropout_pct=as.numeric(dropout_pct))

color_pallete <- c("TP" = "#FBB040", "FP" = "#0073B2")

line_plot <- ggplot(df, aes(x=dropout_pct, y=value, color=metric, group=metric), name="Metri") +
             geom_line(size=1.2) +
             geom_point(size = 3) +
             scale_color_manual(values=color_pallete, name="Metric", labels=c("TP"="True positives", "FP"="False positives")) +
             labs(x="Dropout (%)",
                  y="Total predicted positives") +
             scale_y_continuous(expand=c(0, 0), breaks=seq(0, 2000, 500), limits = c(0, 2100)) +
             scale_x_continuous(limits=c(0, 85), breaks=seq(0, 80, 10)) +
             geom_hline(yintercept = 600, color="#FBB040", linetype="dashed", size = 1) +
             theme_classic() +
             theme(axis.line = element_line(size = 1),
                   axis.ticks = element_line(size = 1, color = "black"),
                   axis.ticks.length = unit(0.25, "cm"),
                   axis.title.x = element_text(size = 18),
                   axis.text.x  = element_text(size = 18, color = "black"),
                   axis.title.y = element_text(size = 18),
                   axis.text.y  = element_text(size = 18, color = "black"),
                   plot.title   = element_blank(),
                   legend.text  = element_text(size = 18),
                   legend.title = element_text(size = 18))
line_plot

ggsave(paste(output_dir, "fig3-ld_simulation_line_plot.pdf", sep=""),
       plot=line_plot,
       height=5,
       width=8 ,
       units=c("in"),
       dpi=600)

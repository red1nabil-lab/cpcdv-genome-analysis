#!/usr/bin/env Rscript
# Step 17: plot the CpCDV phylogeny. Uses the ML tree if it exists, otherwise
# builds a neighbour-joining tree (TN93) as a stand-in so the figure pipeline
# can be checked before the ML run finishes.
suppressMessages({library(ape); library(phangorn); library(ggtree)
                  library(ggtreeExtra); library(ggnewscale)
                  library(ggplot2); library(jsonlite); library(dplyr)})
setwd(Sys.getenv("CPCDV_DATA",
                 unset = file.path(dirname(getwd()), "data")))
options(ignore.negative.edge = TRUE)

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) > 0) args[1] else "auto"

keep <- fromJSON("10_analysis_set.json")$keep
meta <- read.csv("03_curated_metadata.csv", fileEncoding = "UTF-8-BOM",
                 stringsAsFactors = FALSE)
sass <- read.csv("06_strain_assignment.csv", fileEncoding = "UTF-8-BOM",
                 stringsAsFactors = FALSE)

mlfile <- "12_ml_tree.newick"
if (mode != "nj" && file.exists(mlfile)) {
  tr <- read.tree(mlfile); src <- "maximum likelihood (GTR+G+I), ultrafast bootstrap"
  tag <- "ML"
} else {
  aln <- read.phyDat("05_full_aln.fasta", format = "fasta", type = "DNA")
  aln <- subset(aln, keep)
  d <- dist.ml(aln, model = "F81")
  tr <- NJ(d); src <- "neighbour-joining (F81 distance)"; tag <- "NJ"
}
tr <- phangorn::midpoint(tr)
cat("tree:", length(tr$tip.label), "tips -", src, "\n")

# ---- tip annotation -----------------------------------------------------
m <- meta[match(tr$tip.label, meta$accession), ]
s <- sass[match(tr$tip.label, sass$accession), ]
dat <- data.frame(label = tr$tip.label,
                  strain = s$strain_final,
                  region = m$subregion,
                  country = m$country,
                  host = m$host_family,
                  stringsAsFactors = FALSE)
dat$strain[is.na(dat$strain)] <- "unassigned"

# fixed hue order, never cycled; rarer strains fold into "other"
topS <- names(sort(table(dat$strain[dat$strain != "unassigned"]),
                   decreasing = TRUE))[1:7]
dat$strain_grp <- ifelse(dat$strain %in% topS, paste("strain", dat$strain),
                  ifelse(dat$strain == "unassigned", "unassigned", "other strain"))
lev <- c(paste("strain", topS), "other strain", "unassigned")
dat$strain_grp <- factor(dat$strain_grp, levels = lev)
CAT <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
         "#008300", "#6250d6")
cols <- setNames(c(CAT[seq_along(topS)], "#b8b7b1", "#6b6a66"), lev)

REG <- c("East Africa" = "#2a78d6", "South Asia" = "#eb6834",
         "West Asia" = "#1baf7a", "North Africa" = "#eda100",
         "Arabian Peninsula" = "#e87ba4", "West Africa" = "#008300",
         "Southern Africa" = "#6250d6")

# Two panels of the same tree, rather than one tree carrying a ring as well:
# a ring plus coloured tips would give one hue two meanings (blue would be
# both strain H and East Africa).
base_theme <- theme(
  legend.position = "bottom",
  legend.key.size = unit(7, "pt"),
  legend.text = element_text(size = 6.2, colour = "#0b0b0b"),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.margin = margin(-8, 0, 0, 0),
  plot.title = element_text(size = 9, colour = "#0b0b0b", hjust = 0),
  plot.background = element_rect(fill = "#fcfcfb", colour = NA),
  panel.background = element_rect(fill = "#fcfcfb", colour = NA),
  panel.spacing = unit(0, "pt"),
  plot.margin = margin(1, 1, 1, 1))

mk <- function(var, values, title, nrow, support = FALSE) {
  g <- ggtree(tr, layout = "circular", size = 0.2, colour = "#a8a7a2")
  if (support) {
    g <- g %<+% supdf +
      geom_nodepoint(aes(subset = !is.na(bs_num) & bs_num >= 0.95),
                     shape = 16, size = 0.45, colour = "#52514e")
  }
  g <- g %<+% dat
  g + geom_tippoint(aes(colour = .data[[var]]), size = 1.1, stroke = 0) +
    scale_colour_manual(values = values, name = NULL, na.translate = FALSE) +
    guides(colour = guide_legend(nrow = nrow, byrow = TRUE,
                                 override.aes = list(size = 1.9))) +
    ggtitle(title) + base_theme
}

# bootstrap support, keyed by node number for ggtree's %<+% join
bs <- suppressWarnings(as.numeric(tr$node.label))
if (length(bs) && max(bs, na.rm = TRUE) > 1.5) bs <- bs / 100   # accept 0-100
supdf <- data.frame(node = seq_len(length(tr$tip.label) + tr$Nnode),
                    bs_num = c(rep(NA_real_, length(tr$tip.label)), bs))
cat("nodes with bootstrap >= 0.95:", sum(bs >= 0.95, na.rm = TRUE),
    "of", length(bs), "\n")
dat$region_f <- factor(dat$region, levels = names(REG))
pa <- mk("strain_grp", cols, "a   tips coloured by strain", 3, support = TRUE)
pb <- mk("region_f", REG, "b   tips coloured by geographic subregion", 3, support = TRUE)

comb <- try({
  library(patchwork)
  pa + pb + plot_layout(ncol = 2)
}, silent = TRUE)
if (inherits(comb, "try-error")) comb <- pa

ggsave(sprintf("figures/Fig_S5_phylogeny_%s.pdf", tag), comb,
       width = 7.4, height = 3.0)
ggsave(sprintf("figures/Fig_S5_phylogeny_%s.png", tag), comb,
       width = 7.4, height = 3.0, dpi = 300)
cat("wrote figures/Fig_S5_phylogeny_", tag, "\n", sep = "")

# ---- how well do strains form clades? -----------------------------------
cat("\nmonophyly of assigned strains:\n")
res <- data.frame()
for (st in sort(unique(dat$strain))) {
  tips <- dat$label[dat$strain == st]
  if (length(tips) < 2) next
  mono <- is.monophyletic(tr, tips)
  cat(sprintf("  strain %-12s n=%-4d monophyletic: %s\n", st, length(tips), mono))
  res <- rbind(res, data.frame(strain = st, n = length(tips), monophyletic = mono))
}
write.csv(res, "17_strain_monophyly.csv", row.names = FALSE)
cat("wrote 17_strain_monophyly.csv\n")

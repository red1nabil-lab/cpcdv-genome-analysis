#!/usr/bin/env Rscript
# Step 12: maximum-likelihood phylogeny of the CpCDV full-genome analysis set,
# with ultrafast bootstrap support, plus a model-selection table.
suppressMessages({library(phangorn); library(ape); library(jsonlite)})
setwd(Sys.getenv("CPCDV_DATA",
                 unset = file.path(dirname(getwd()), "data")))

keep <- fromJSON("10_analysis_set.json")$keep
aln  <- read.phyDat("05_full_aln.fasta", format = "fasta", type = "DNA")
aln  <- subset(aln, keep)
cat("taxa:", length(aln), " sites:", sum(attr(aln, "weight")), "\n")

# ---- model selection ----------------------------------------------------
cat("\n== model test ==\n")
mt <- modelTest(aln, model = c("JC", "K80", "HKY", "TrN", "GTR"),
                G = TRUE, I = TRUE, k = 4, multicore = TRUE, mc.cores = 8)
mt <- mt[order(mt$BIC), ]
print(head(mt[, c("Model", "df", "logLik", "AIC", "BIC")], 8))
write.csv(mt, "12_model_selection.csv", row.names = FALSE)
best <- mt$Model[1]
cat("\nbest model by BIC:", best, "\n")

# ---- ML tree with ultrafast bootstrap ----------------------------------
set.seed(42)
cat("\n== ML inference with ultrafast bootstrap (1000 replicates) ==\n")
fit <- pml_bb(aln, model = best, rearrangement = "stochastic",
              control = pml.control(trace = 0))
cat("final logLik:", logLik(fit), "\n")

tr <- fit$tree
write.tree(tr, "12_ml_tree.newick")
cat("wrote 12_ml_tree.newick\n")

# midpoint-rooted copy for display
tr2 <- try(phangorn::midpoint(tr), silent = TRUE)
if (!inherits(tr2, "try-error")) {
  write.tree(tr2, "12_ml_tree_midpoint.newick")
  cat("wrote 12_ml_tree_midpoint.newick\n")
}

# ---- patristic distance summary ----------------------------------------
pd <- cophenetic(tr)
cat("\npatristic distance: mean", round(mean(pd[upper.tri(pd)]), 4),
    " max", round(max(pd), 4), "\n")
saveRDS(fit, "12_pml_fit.rds")
cat("done\n")

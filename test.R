rm(list = ls())

library(Informeasure)
library(sycomore)

load(system.file("extdata/tcga.brca.testdata.Rdata", package="sycomore"))

miRNAexpression  <- log2(miRNAexpression  + 1)
mRNAexpression   <- log2(mRNAexpression   + 1)

results <- sycomore(miRtarget, miRNAexpression, mRNAexpression,
                    permutation = TRUE,
                    B = 10000,
                    permute_target = "mRNA",
                    one_sided = TRUE,
                    seed = 1)

results_synergy <- results %>% filter(interaction == "neighboring" & SRscore > 0 & q.SRscore < 0.05)
results_redundancy <- results %>% filter(interaction == "overlapping" & SRscore < 0 & q.SRscore < 0.05)

library(GenomicRanges)
library(dplyr)
library(data.table)
library(R.utils)

source("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/code/benchmark_EG_finemapUKBB_gene_NUMEG.R")

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
EGpred_file = toString(args[1])
tissue = toString(args[2])
output_file = toString(args[3])


oo1 = benchmark_EG_finemapUKBB_gene(EGpred_file, tissue = tissue, use_pops=TRUE, num_pops=1, num_EGgenes=1)
oo2 = benchmark_EG_finemapUKBB_gene(EGpred_file, tissue = tissue, use_pops=TRUE, num_pops=2, num_EGgenes=2)
oo3 = benchmark_EG_finemapUKBB_gene(EGpred_file, tissue = tissue, use_pops=TRUE, num_pops=2, num_EGgenes=1)
oo4 = benchmark_EG_finemapUKBB_gene(EGpred_file, tissue = tissue, use_pops=TRUE, num_pops=1, num_EGgenes = 500)
oo5 = benchmark_EG_finemapUKBB_gene(EGpred_file, tissue = tissue, use_pops=TRUE, num_pops=2, num_EGgenes = 500)
oo6 = benchmark_EG_finemapUKBB_gene(EGpred_file, tissue = tissue, use_pops=FALSE, num_EGgenes=1)
oo7 = benchmark_EG_finemapUKBB_gene(EGpred_file, tissue = tissue, use_pops=FALSE, num_EGgenes=2)
oo8 = benchmark_EG_finemapUKBB_gene(EGpred_file, tissue = tissue, use_pops=FALSE, num_EGgenes=500)
oo9 = benchmark_EG_finemapUKBB_gene(EGpred_file, tissue = tissue, use_pops=TRUE, num_pops=1, num_EGgenes=2)

outdff = rbind(c(oo1$precision, oo1$sd.precision,  oo1$recall, oo1$sd.recall,  oo1$precision2, oo1$sd.precision2,  oo1$recall2, oo1$sd.recall2),
               c(oo2$precision, oo2$sd.precision,  oo2$recall, oo2$sd.recall,  oo2$precision2, oo2$sd.precision2,  oo2$recall2, oo2$sd.recall2),
               c(oo3$precision, oo3$sd.precision,  oo3$recall, oo3$sd.recall,  oo3$precision2, oo3$sd.precision2,  oo3$recall2, oo3$sd.recall2),
               c(oo4$precision, oo4$sd.precision,  oo4$recall, oo4$sd.recall,  oo4$precision2, oo4$sd.precision2,  oo4$recall2, oo4$sd.recall2),
               c(oo5$precision, oo5$sd.precision,  oo5$recall, oo5$sd.recall,  oo5$precision2, oo5$sd.precision2,  oo5$recall2, oo5$sd.recall2),
               c(oo6$precision, oo6$sd.precision,  oo6$recall, oo6$sd.recall,  oo6$precision2, oo6$sd.precision2,  oo6$recall2, oo6$sd.recall2),
               c(oo7$precision, oo7$sd.precision,  oo7$recall, oo7$sd.recall,  oo7$precision2, oo7$sd.precision2,  oo7$recall2, oo7$sd.recall2),
               c(oo8$precision, oo8$sd.precision,  oo8$recall, oo8$sd.recall,  oo8$precision2, oo8$sd.precision2,  oo8$recall2, oo8$sd.recall2),
               c(oo9$precision, oo9$sd.precision,  oo9$recall, oo9$sd.recall,  oo9$precision2, oo9$sd.precision2,  oo9$recall2, oo9$sd.recall2))

colnames(outdff) = c("Precision","sd.Precision", "Recall", "sd.Recall", "Precision2", "sd.Precision2", "Recall2", "sd.Recall2")
rownames(outdff) = c("PoPS:1; EG:1", "PoPS:2; EG:2", "PoPS:2; EG:1", "PoPS:1; EG:50", "PoPS:2; EG:50", "EG:1", "EG:2", "EG:50",
                     "PoPS:1; EG:2")
write.table(outdff, file = paste0(output_file), col.names = T, row.names = T, sep = "\t", quote=F)



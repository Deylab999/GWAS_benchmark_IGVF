library(data.table)
library(R.utils)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
genescore_dir <- toString(args[1])
bed_dir <- toString(args[2])
annot_name <- toString(args[3])

score_file = paste0(genescore_dir, "/", annot_name, ".txt")
gene_scores = read.delim(score_file, header=F)

if(!dir.exists(paste0(bed_dir, "/", annot_name))){
  dir.create(paste0(bed_dir, "/", annot_name))
}

source("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/code/Primary_CRE_gene_linkfiles_ENCODE_E2G_2024.R")

gene_scores = read.delim(score_file, header=F)

scores = gene_scores[,2]
names(scores) = gene_scores[,1]

#out1 = ENCODE_E2G_2024_DNase_calc(scores,
#                                  pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_Only/thresholded_predictions",
#                                  output_cell = paste0(bed_dir, "/", annot_name),
#                                  output_bed = paste0("ENCODE_E2G_DNase_only"))

#out1 = ENCODE_E2G_2024_calc(scores,
#                            pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_HiC/thresholded_predictions",
#                            output_cell = paste0(bed_dir, "/", annot_name),
#                            output_bed = paste0("ENCODE_E2G_DNase_HiC"))

#out1 = ABC_2024_DNase_calc(scores,
#                            pred_cell="/data/deyk/ENCODE/ABC_2024/thresholded_predictions",
#                            output_cell = paste0(bed_dir, "/", annot_name),
#                            output_bed = paste0("ABC_2024_DNase_only"))


out1 = ENCODE_E2G_2024_ext_calc(scores,
                           pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/predictions_ext",
                           output_cell = paste0(bed_dir, "/", annot_name),
                           output_bed = paste0("ENCODE_E2G_2024_ext"))

out1 = ABC_2024_ext_calc(scores,
                        pred_cell="/data/deyk/ENCODE/ABC_2024/predictions_ext",
                        output_cell = paste0(bed_dir, "/", annot_name),
                        output_bed = paste0("ABC_2024_ext"))

out1 = scE2G_calc(scores,
                 pred_cell="/data/deyk/IGVF/IGVF_scE2G",
                 output_cell = paste0(bed_dir, "/", annot_name),
                 output_bed = paste0("scE2G"))


ENCODE_E2G_2024_DNase_BLD_EG_preds <- function(scores,
                                       pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_Only/thresholded_predictions",
                                       output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/ENCODE_E2G_Blood_EG_predictions.txt"){
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(logreg_dnase_biosamples[c(grep("T-", logreg_dnase_biosamples),
                                               grep("T_cell", logreg_dnase_biosamples),
                                               grep("T_follicular_helper_cell", logreg_dnase_biosamples),
                                               grep("hemato", logreg_dnase_biosamples),
                                               grep("B_cell", logreg_dnase_biosamples),
                                               grep("_CD", logreg_dnase_biosamples),
                                               grep("B-", logreg_dnase_biosamples),
                                               grep("monocyte", logreg_dnase_biosamples),
                                               grep("killer", logreg_dnase_biosamples),
                                               grep("GM1", logreg_dnase_biosamples),
                                               grep("K562", logreg_dnase_biosamples),
                                               grep("myeloid", logreg_dnase_biosamples),
                                               grep("lymphoid", logreg_dnase_biosamples))])
  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
    tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numl]] = tabb2
    cat("We are at blood sample:", numl, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out1 = ENCODE_E2G_2024_DNase_BLD_EG_preds()

dist_tss_BLD_EG_preds <- function(scores,
                                  pred_cell="/data/deyk/ENCODE/Baseline_Preds_2024/",
                                  output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/dist_TSS_Blood_EG_predictions.txt"){
  library(data.table)
  logreg_dnase_biosamples = read.table("/data/deyk/ENCODE/Baseline_Preds_2024/preds.txt", header = F)[,1]
  keep_eids = unique(logreg_dnase_biosamples[c(grep("T-", logreg_dnase_biosamples),
                                               grep("T_cell", logreg_dnase_biosamples),
                                               grep("T_follicular_helper_cell", logreg_dnase_biosamples),
                                               grep("hemato", logreg_dnase_biosamples),
                                               grep("B_cell", logreg_dnase_biosamples),
                                               grep("_CD", logreg_dnase_biosamples),
                                               grep("B-", logreg_dnase_biosamples),
                                               grep("monocyte", logreg_dnase_biosamples),
                                               grep("killer", logreg_dnase_biosamples),
                                               grep("GM1", logreg_dnase_biosamples),
                                               grep("K562", logreg_dnase_biosamples),
                                               grep("myeloid", logreg_dnase_biosamples),
                                               grep("lymphoid", logreg_dnase_biosamples))])
  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0("/data/deyk/ENCODE/Baseline_Preds_2024/",
                                         keep_eids[numl], "/", "dist_to_tss", ".tsv.gz")))
    preds_tabb = preds_tabb[which(preds_tabb$Score < 53976), ]
    tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, 1)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numl]] = tabb2
    cat("We are at blood sample:", numl, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out1 = dist_tss_BLD_EG_preds()

nearest_gene_BLD_EG_preds <- function(scores,
                                  pred_cell="/data/deyk/ENCODE/Baseline_Preds_2024/",
                                  output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/nearest_gene_Blood_EG_predictions.txt"){
  library(data.table)
  logreg_dnase_biosamples = read.table("/data/deyk/ENCODE/Baseline_Preds_2024/preds.txt", header = F)[,1]
  keep_eids = unique(logreg_dnase_biosamples[c(grep("T-", logreg_dnase_biosamples),
                                               grep("T_cell", logreg_dnase_biosamples),
                                               grep("T_follicular_helper_cell", logreg_dnase_biosamples),
                                               grep("hemato", logreg_dnase_biosamples),
                                               grep("B_cell", logreg_dnase_biosamples),
                                               grep("_CD", logreg_dnase_biosamples),
                                               grep("B-", logreg_dnase_biosamples),
                                               grep("monocyte", logreg_dnase_biosamples),
                                               grep("killer", logreg_dnase_biosamples),
                                               grep("GM1", logreg_dnase_biosamples),
                                               grep("K562", logreg_dnase_biosamples),
                                               grep("myeloid", logreg_dnase_biosamples),
                                               grep("lymphoid", logreg_dnase_biosamples))])
  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0("/data/deyk/ENCODE/Baseline_Preds_2024/",
                                         keep_eids[numl], "/", "nearest_expressed_gene", ".tsv.gz")))
    preds_tabb = preds_tabb[which(preds_tabb$Score > 0), ]
    tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, 1)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numl]] = tabb2
    cat("We are at blood sample:", numl, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out1 = nearest_gene_BLD_EG_preds()


ABC_2024_DNase_BLD_EG_preds <- function(scores,
                                        pred_cell="/data/deyk/ENCODE/ABC_2024/thresholded_predictions",
                                        output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/ABC_Blood_EG_predictions.txt"){
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(logreg_dnase_biosamples[c(grep("T-", logreg_dnase_biosamples),
                                               grep("T_cell", logreg_dnase_biosamples),
                                               grep("T_follicular_helper_cell", logreg_dnase_biosamples),
                                               grep("hemato", logreg_dnase_biosamples),
                                               grep("B_cell", logreg_dnase_biosamples),
                                               grep("_CD", logreg_dnase_biosamples),
                                               grep("B-", logreg_dnase_biosamples),
                                               grep("monocyte", logreg_dnase_biosamples),
                                               grep("killer", logreg_dnase_biosamples),
                                               grep("GM1", logreg_dnase_biosamples),
                                               grep("K562", logreg_dnase_biosamples),
                                               grep("myeloid", logreg_dnase_biosamples),
                                               grep("lymphoid", logreg_dnase_biosamples))])
  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
    tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numl]] = tabb2
    cat("We are at blood sample:", numl, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out1 = ABC_2024_DNase_BLD_EG_preds()


ENCODE_E2G_2024_ext_BLD_EG_preds <- function(scores,
                                               pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_HiC/thresholded_predictions",
                                               output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/ENCODE_E2G_ext_Blood_EG_predictions.txt"){
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(logreg_dnase_biosamples[c(grep("T-", logreg_dnase_biosamples),
                                               grep("T_cell", logreg_dnase_biosamples),
                                               grep("T_follicular_helper_cell", logreg_dnase_biosamples),
                                               grep("hemato", logreg_dnase_biosamples),
                                               grep("B_cell", logreg_dnase_biosamples),
                                               grep("_CD", logreg_dnase_biosamples),
                                               grep("B-", logreg_dnase_biosamples),
                                               grep("monocyte", logreg_dnase_biosamples),
                                               grep("killer", logreg_dnase_biosamples),
                                               grep("GM1", logreg_dnase_biosamples),
                                               grep("K562", logreg_dnase_biosamples),
                                               grep("myeloid", logreg_dnase_biosamples),
                                               grep("lymphoid", logreg_dnase_biosamples))])
  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
    tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numl]] = tabb2
    cat("We are at blood sample:", numl, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out1 = ENCODE_E2G_2024_ext_BLD_EG_preds()

ENCODE_E2G_2024_DNase_All_EG_preds <- function(scores,
                                               pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_Only/thresholded_predictions",
                                               output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/ENCODE_E2G_All_EG_predictions.txt"){
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(logreg_dnase_biosamples)
  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
    tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numl]] = tabb2
    cat("We are at blood sample:", numl, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out2 = ENCODE_E2G_2024_DNase_All_EG_preds()


ENCODE_E2G_2024_ext_All_EG_preds <- function(scores,
                                             pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_HiC/thresholded_predictions",
                                             output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/ENCODE_E2G_ext_All_EG_predictions.txt"){
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(logreg_dnase_biosamples)
  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
    tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numl]] = tabb2
    cat("We are at blood sample:", numl, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out2 = ENCODE_E2G_2024_ext_All_EG_preds()

ABC_2024_All_EG_preds <- function(scores,
                                      pred_cell="/data/deyk/ENCODE/ABC_2024/thresholded_predictions",
                                      output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/ABC_All_EG_predictions.txt"){
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  keep_eids = unique(logreg_dnase_biosamples)
  pooled_tabb_list = list()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
    tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numl]] = tabb2
    cat("We are at sample:", numl, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out2 = ABC_2024_All_EG_preds()

ll = list.files("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions")
oo = cbind.data.frame(ll, c("ALL", "BLD", "ALL", "BLD", "ALL", "BLD"),  paste0(ll, ".prec_recall"))
write.table(oo, file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/finemap_precision_recall.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)



Ramil_BLD_EG_preds <- function(
                       pred_file="/n/groups/price/kushal/ENCODE/data/EG_predictions/Ramil_predictions.complete.thresh0_009.tsv.gz",
                       key_file="/n/groups/price/kushal/ENCODE/data/EG_predictions/Ramil_biosamples_key.csv",
                       output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/Ramil_Blood_EG_predictions.txt")
{
    library(data.table)
    ramil_biosamples = read.csv(paste0(key_file))
    tabb_pre = data.frame(fread(paste0(pred_file)))
    keep_eids = ramil_biosamples[grep("BLD", ramil_biosamples[,2]), 1]

    tabb = tabb_pre[which(tabb_pre$CellType %in% keep_eids == T), ]
    tabb2 = tabb[, c("chr", "start", "end", "TargetGene", "Score")]
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    fwrite(tabb2, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


BingABC_BLD_EG_preds <- function(
                         pred_file="/n/groups/price/kushal/ENCODE/data/EG_predictions/2022_BingRen_scATAC_EnhancerPredictionsFull.txt.gz",
                         key_file="/n/groups/price/kushal/ENCODE/data/EG_predictions/BingABC_biosamples.csv",
                         output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/BingABC_Blood_EG_predictions.txt")
{
  library(data.table)
  bingabc_biosamples = read.csv(paste0(key_file))
  tabb_pre = data.frame(fread(paste0(pred_file)))
  tabb_pre = tabb_pre[which(tabb_pre$Score > 0.025), ]
  keep_eids = bingabc_biosamples[grep("BLD", bingabc_biosamples[,2]), 1]

  tabb = tabb_pre[which(tabb_pre$CellType %in% keep_eids == T), ]
  tabb2 = cbind.data.frame(tabb$chr, tabb$start, tabb$end, tabb$TargetGene, tabb$Score)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(tabb2, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

EpiMap_BLD_EG_preds <- function(scores,
                        pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/EpiMap",
                        key_file="/n/groups/price/kushal/ENCODE/data/EG_predictions/EpiMap/EpiMap_biosample_key.tsv",
                        output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/EpiMap_Blood_EG_predictions.txt"){
    library(data.table)
    biosamples_listed = as.character(sapply(as.character(sapply(list.files(pred_cell, pattern = "linking"),
                                                                function(x) return(strsplit(x, "_")[[1]][3]))),
                                            function(xx) return(strsplit(xx, "[.]")[[1]][1])))
    epimap_biosamples = read.delim(paste0(key_file))
    keep_eids = epimap_biosamples$id[c(grep("Blood & T-cell", epimap_biosamples$GROUP),
                                       grep("HSC & B-cell", epimap_biosamples$GROUP),
                                       grep("Lymphoblastoid", epimap_biosamples$GROUP))]

    keep_eids = intersect(keep_eids, biosamples_listed)

    pooled_list = list()
    for(numl in 1:length(keep_eids)){
      preds_tabb = data.frame(fread(paste0("/n/groups/price/kushal/ENCODE/data/EG_predictions/EpiMap/",
                                           "linking_collated_", keep_eids[numl],".bed.gz")))
      preds_tabb = preds_tabb[which(preds_tabb$Score > 0.0045), ]
      tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
      colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
      pooled_list[[numl]] = tabb2
      cat("Read file", numl, "out of ", length(keep_eids), " files \n")
    }
    pooled_tabb = do.call(rbind, pooled_list)
    fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


Baseline_BLD_EG_preds <- function(
                             pred_type = "dist_to_tss",
                             preds_file="/n/groups/price/kushal/ENCODE/data/EG_predictions/Baseline_EG/preds.txt",
                             output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/Baseline_dist_to_tss_Blood_EG_predictions.txt")
{
  library(data.table)
  preds_biosample = read.table(preds_file, header = F)[,1]
  keep_eids = c("CD14-positive_monocyte_ID_370", "CD4-positive_alpha-beta_T_cell_ID_120", "CD8-positive_alpha-beta_T_cell_ID_148",
                "GM12878_ID_2724", "GM23338_ID_2088", "K562_ID_2644", "natural_killer_cell_ID_2608", "T-cell_ID_61")

  pooled_tabb = c()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0("/n/groups/price/kushal/ENCODE/data/EG_predictions/Baseline_EG/",
                                         keep_eids[numl], "/", pred_type, ".tsv.gz")))
    if(pred_type == "dist_to_tss"){
      preds_tabb = preds_tabb[which(preds_tabb$Score < 58111), ]
    }
    if(pred_type == "dist_to_gene"){
      preds_tabb = preds_tabb[which(preds_tabb$Score < 44686), ]
    }
    if(pred_type == "reads_by_dist_to_tss_norm"){
      preds_tabb = preds_tabb[which(preds_tabb$Score > 0.002), ]
    }
    if(pred_type == "H3K27ac_reads_by_dist_to_tss_norm"){
      preds_tabb = preds_tabb[which(preds_tabb$Score > 0.002), ]
    }
    tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb = rbind(pooled_tabb, tabb2)
  }
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


LogReg_DNase_EG_preds = function(pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/LogRegClassifiers_DNase",
                                 output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/LogReg_DNase_Blood_EG_predictions.txt")
{
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "features.predictions")
  keep_eids = logreg_dnase_biosamples[c(grep("T-", logreg_dnase_biosamples),
                                        grep("T_cell", logreg_dnase_biosamples),
                                        grep("T_follicular_helper_cell", logreg_dnase_biosamples),
                                        grep("B_cell", logreg_dnase_biosamples),
                                        grep("B-", logreg_dnase_biosamples),
                                        grep("monocyte", logreg_dnase_biosamples),
                                        grep("GM1", logreg_dnase_biosamples),
                                        grep("K562", logreg_dnase_biosamples),
                                        grep("myeloid", logreg_dnase_biosamples),
                                        grep("lymphoid", logreg_dnase_biosamples))]
  pooled_tabb = c()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
    preds_tabb = preds_tabb[which(preds_tabb$Full.Score > 0.20), ]
    tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Full.Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb = rbind(pooled_tabb, tabb2)
    cat("We are at celltype:", numl, "\n")
  }

  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


LogReg_DNase_All_EG_preds = function(pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/LogRegClassifiers_DNase",
                                 output_file = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/LogReg_DNase_All_EG_predictions.txt")
{
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "features.predictions")
  keep_eids = logreg_dnase_biosamples
  pooled_tabb = c()
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
    preds_tabb = preds_tabb[which(preds_tabb$Full.Score > 0.20), ]
    tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Full.Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb = rbind(pooled_tabb, tabb2)
    cat("We are at celltype:", numl, "\n")
  }

  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


LogReg_DNase_Biosample_EG_preds = function(pred_cell="/n/groups/price/kushal/ENCODE/data/EG_predictions/LogRegClassifiers_DNase",
                                     output_cell = "/n/groups/price/kushal/ENCODE/data/EG_predictions/Biosample_EG_predictions")
{
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "features.predictions")
  logreg_dnase_biosamples2 = as.character(sapply(logreg_dnase_biosamples, function(x) return(strsplit(x, "_DNaseOnly")[[1]][1])))

  keep_eids = logreg_dnase_biosamples
  for(numl in 1:length(keep_eids)){
    preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
    preds_tabb = preds_tabb[which(preds_tabb$Full.Score > 0.20), ]
    tabb2 = cbind.data.frame(preds_tabb$chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene, preds_tabb$Full.Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    fwrite(tabb2, paste0(output_cell, "/", "LogReg_DNase_", logreg_dnase_biosamples2[numl], "_EG_predictions.txt"),
           sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    cat("We are at celltype:", numl, "\n")
  }
}



###########################################  Output data files   ####################################################


features = c("H3K27ac_reads_by_dist_to_tss_norm",
             "dist_to_gene",
             "dist_to_tss",
             "nearest_gene",
             "nearest_tss",
             "reads_by_dist_to_tss_norm",
             "within_100kb_of_tss")

for(numf in 1:length(features)){
  out1 = Baseline_BLD_EG_preds (pred_type = features[numf],
                                output_file = paste0("/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/Baseline_",
                                features[numf], "_Blood_EG_predictions.txt"))
  cat("We are at feature:", features[numf])
}

out2 = Ramil_BLD_EG_preds()
out3 = BingABC_BLD_EG_preds()
out4 = EpiMap_BLD_EG_preds()
out5 = LogReg_DNase_EG_preds()

out6 = LogReg_DNase_Biosample_EG_preds()



ll = list.files("/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions")
oo = cbind.data.frame(ll, "BLD", paste0(ll, ".prec_recall"))
write.table(oo, file = "/n/groups/price/kushal/ENCODE/code/Figure4/Primary_onlyEG/finemap_precision_recall.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)


dff = data.frame(fread(paste0("/n/groups/price/kushal/ENCODE/data/EG_predictions/Blood_EG_predictions/Baseline_",
                              "nearest_tss", "_Blood_EG_predictions.txt")))


ENCODE_E2G_2024_ext_EG_preds <- function(
                        pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/predictions_ext",
                        output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/ENCODE_E2G_extended_Blood_EG_predictions.txt"){
  library(data.table)
  tissuenames = c("K562", "GM12878")
  pooled_tabb_list = list()
  for(numtt in 1:length(tissuenames)){
    if(tissuenames[numtt] == "K562"){
      tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562", "/", "encode_e2g_predictions_threshold0.336.tsv.gz")))
    }else if(tissuenames[numtt] == "GM12878"){
      tabb_pre = data.frame(fread(paste0(pred_cell, "/", "GM12878", "/", "encode_e2g_predictions_threshold0.336.tsv.gz")))
    }
    tabb2 = cbind.data.frame(tabb_pre$chr, tabb_pre$start, tabb_pre$end, tabb_pre$TargetGene, tabb_pre$ENCODE.rE2G.Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numtt]] = tabb2
    cat("We are at blood sample:", numtt, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out2 = ENCODE_E2G_2024_ext_EG_preds()


ABC_2024_ext_EG_preds <- function(
                    pred_cell="/data/deyk/ENCODE/ABC_2024/predictions_ext",
                    output_file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/ABC_extended_Blood_EG_predictions.txt"){
  library(data.table)
  tissuenames = c("K562", "GM12878")
  pooled_tabb_list = list()
  for(numtt in 1:length(tissuenames)){
    if(tissuenames[numtt] == "K562"){
      tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562/Predictions", "/", "ABC_pred_filtered_0.027_selfPromoters.tsv.gz")))
    }else if(tissuenames[numtt] == "GM12878"){
      tabb_pre = data.frame(fread(paste0(pred_cell, "/", "GM12878/Predictions/", "/", "ABC_pred_filtered_0.027_selfPromoters.tsv.gz")))
    }
    tabb2 = cbind.data.frame(tabb_pre$chr, tabb_pre$start, tabb_pre$end, tabb_pre$TargetGene, tabb_pre$ABC.Score)
    colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
    pooled_tabb_list[[numtt]] = tabb2
    cat("We are at blood sample:", numtt, "\n")
  }
  pooled_tabb = do.call(rbind, pooled_tabb_list)
  fwrite(pooled_tabb, paste0(output_file), sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

out2 = ABC_2024_ext_EG_preds()


scE2G_2024_ext_EG_preds <- function(
                       pred_cell="/data/deyk/IGVF/IGVF_scE2G",
                       output_cell,
                       output_bed = "temp.bed"){

  ##############################. multiome-7features-ENCODE-E2G  ##################################################

  tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562_Xu_multiome_7features_genomewide_predictions.tsv.gz")))
  tabb_pre = tabb_pre[which(tabb_pre$ENCODE.rE2G.Score > 0.171), ]
  tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,6], tabb_pre$ENCODE.rE2G.Score)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(tabb2, paste0("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/K562_Xu_multiome_7features_genomewide_predictions.txt"),
         sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

  ##############################. multiome-7features-intx-ENCODE-E2G  ##################################################

  tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562_Xu_multiome_7features_intx_genomewide_predictions.tsv.gz")))
  tabb_pre = tabb_pre[which(tabb_pre$ENCODE.rE2G.Score > 0.176), ]
  tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,6], tabb_pre$ENCODE.rE2G.Score)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(tabb2, paste0("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/K562_Xu_multiome_7features_intx_genomewide_predictions.txt"),
         sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

  ##############################. scATAC-6features-ENCODE-E2G  ##################################################

  tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562_Xu_scATAC_6features_genomewide_predictions.tsv.gz")))
  tabb_pre = tabb_pre[which(tabb_pre$ENCODE.rE2G.Score > 0.154), ]
  tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,6], tabb_pre$ENCODE.rE2G.Score)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(tabb2, paste0("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/K562_Xu_scATAC_6features_genomewide_predictions.txt"),
         sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

  ##############################. scATAC-6features-intx-ENCODE-E2G  ##################################################

  tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562_Xu_scATAC_6features_intx_genomewide_predictions.tsv.gz")))
  tabb_pre = tabb_pre[which(tabb_pre$ENCODE.rE2G.Score > 0.22), ]
  tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,6], tabb_pre$ENCODE.rE2G.Score)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(tabb2, paste0("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/K562_Xu_scATAC_6features_intx_genomewide_predictions.txt"),
         sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


  ##############################. scATAC-ABC.  ##################################################

  tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562_Xu_scATAC_6features_genomewide_predictions.tsv.gz")))
  tabb_pre = tabb_pre[which(tabb_pre$ABC.Score > 0.015), ]
  tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,6], tabb_pre$ABC.Score)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(tabb2, paste0("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/K562_Xu_scATAC_ABC_6features_genomewide_predictions.txt"),
         sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


  ##################### Kendall correlations ##################################

  tabb_pre = data.frame(fread(paste0(pred_cell, "/", "Pairs.Kendall.tsv.gz")))
  tabb_pre = tabb_pre[which(tabb_pre$Kendall > 0.011), ]
  tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,4], tabb_pre$Kendall)
  colnames(tabb2) = c("chr", "start", "end", "TargetGene", "Score")
  fwrite(tabb2, paste0("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions/K562_Pairs_Kendall.txt"),
         sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}


out2 = scE2G_2024_ext_EG_preds()

ll = list.files("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions")
oo = cbind.data.frame(ll, rep(c("BLD"), length(ll)),  paste0(ll, ".prec_recall"))
write.table(oo, file = "/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/finemap_precision_recall.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)


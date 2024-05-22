
ENCODE_E2G_2024_DNase_calc <- function(scores,
                              pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_Only/thresholded_predictions",
                              tissuenames = c("ALL", "BLD", "BRN", "K562", "GM12878", "Mono"),
                              output_cell,
                              output_bed = "temp.bed"){
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  for(numtt in 1:length(tissuenames)){
    if(tissuenames[numtt] == "ALL"){
      keep_eids = logreg_dnase_biosamples
    }else if (tissuenames[numtt] == "BLD"){
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
    }else if (tissuenames[numtt] == "BRN"){
      keep_eids = unique(logreg_dnase_biosamples[c(grep("_astro", logreg_dnase_biosamples),
                                            grep("neur", logreg_dnase_biosamples),
                                            grep("cerebell", logreg_dnase_biosamples),
                                            grep("cortex", logreg_dnase_biosamples),
                                            grep("brain", logreg_dnase_biosamples),
                                            grep("glia", logreg_dnase_biosamples))])
    }else if (tissuenames[numtt] == "K562"){
      keep_eids = logreg_dnase_biosamples[grep("K562", logreg_dnase_biosamples)]
    }else if (tissuenames[numtt] == "Mono"){
      keep_eids = logreg_dnase_biosamples[grep("mono", logreg_dnase_biosamples)]
    }else if (tissuenames[numtt] == "GM12878"){
      keep_eids = logreg_dnase_biosamples[grep("GM12878", logreg_dnase_biosamples)]
    }

    pooled_tabb = c()
    for(numl in 1:length(keep_eids)){
      preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
      tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene)
      colnames(tabb2) = c("chr", "start", "end", "TargetGene")
      pooled_tabb = rbind(pooled_tabb, tabb2)
    }

    matched_ids = match(pooled_tabb$TargetGene, names(scores))
    temp = as.numeric(scores)[matched_ids]
    temp[is.na(temp)] = 0
    final_bed1 = cbind(pooled_tabb[,c(1:3)], temp)
    write.table(final_bed1, paste0(output_cell, "/", output_bed, "_", tissuenames[numtt], ".bed"),
                sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}



ENCODE_E2G_2024_calc <- function(scores,
                        pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/DNase_HiC/thresholded_predictions",
                        tissuenames= c("ALL", "BLD", "BRN", "K562", "GM12878", "Mono"),
                        output_cell,
                        output_bed = "temp.bed"){
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  for(numtt in 1:length(tissuenames)){
    if(tissuenames[numtt] == "ALL"){
      keep_eids = logreg_dnase_biosamples
    }else if (tissuenames[numtt] == "BLD"){
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
    }else if (tissuenames[numtt] == "BRN"){
      keep_eids = unique(logreg_dnase_biosamples[c(grep("_astro", logreg_dnase_biosamples),
                                                   grep("neur", logreg_dnase_biosamples),
                                                   grep("cerebell", logreg_dnase_biosamples),
                                                   grep("cortex", logreg_dnase_biosamples),
                                                   grep("brain", logreg_dnase_biosamples),
                                                   grep("glia", logreg_dnase_biosamples))])
    }else if (tissuenames[numtt] == "K562"){
      keep_eids = logreg_dnase_biosamples[grep("K562", logreg_dnase_biosamples)]
    }else if (tissuenames[numtt] == "Mono"){
      keep_eids = logreg_dnase_biosamples[grep("mono", logreg_dnase_biosamples)]
    }else if (tissuenames[numtt] == "GM12878"){
      keep_eids = logreg_dnase_biosamples[grep("GM12878", logreg_dnase_biosamples)]
    }

    pooled_tabb = c()
    for(numl in 1:length(keep_eids)){
      preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
      tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene)
      colnames(tabb2) = c("chr", "start", "end", "TargetGene")
      pooled_tabb = rbind(pooled_tabb, tabb2)
    }

    matched_ids = match(pooled_tabb$TargetGene, names(scores))
    temp = as.numeric(scores)[matched_ids]
    temp[is.na(temp)] = 0
    final_bed1 = cbind(pooled_tabb[,c(1:3)], temp)
    write.table(final_bed1, paste0(output_cell, "/", output_bed, "_", tissuenames[numtt], ".bed"),
                sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}



ABC_2024_DNase_calc <- function(scores,
                                       pred_cell="/data/deyk/ENCODE/ABC_2024/thresholded_predictions",
                                       tissuenames = c("ALL", "BLD", "BRN", "K562", "GM12878", "Mono"),
                                       output_cell,
                                       output_bed = "temp.bed"){
  library(data.table)
  logreg_dnase_biosamples = list.files(pred_cell, pattern = "thresholded_predictions")
  for(numtt in 1:length(tissuenames)){
    if(tissuenames[numtt] == "ALL"){
      keep_eids = logreg_dnase_biosamples
    }else if (tissuenames[numtt] == "BLD"){
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
    }else if (tissuenames[numtt] == "BRN"){
      keep_eids = unique(logreg_dnase_biosamples[c(grep("_astro", logreg_dnase_biosamples),
                                                   grep("neur", logreg_dnase_biosamples),
                                                   grep("cerebell", logreg_dnase_biosamples),
                                                   grep("cortex", logreg_dnase_biosamples),
                                                   grep("brain", logreg_dnase_biosamples),
                                                   grep("glia", logreg_dnase_biosamples))])
    }else if (tissuenames[numtt] == "K562"){
      keep_eids = logreg_dnase_biosamples[grep("K562", logreg_dnase_biosamples)]
    }else if (tissuenames[numtt] == "Mono"){
      keep_eids = logreg_dnase_biosamples[grep("mono", logreg_dnase_biosamples)]
    }else if (tissuenames[numtt] == "GM12878"){
      keep_eids = logreg_dnase_biosamples[grep("GM12878", logreg_dnase_biosamples)]
    }

    pooled_tabb = c()
    for(numl in 1:length(keep_eids)){
      preds_tabb = data.frame(fread(paste0(pred_cell, "/", keep_eids[numl])))
      tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene)
      colnames(tabb2) = c("chr", "start", "end", "TargetGene")
      pooled_tabb = rbind(pooled_tabb, tabb2)
    }

    matched_ids = match(pooled_tabb$TargetGene, names(scores))
    temp = as.numeric(scores)[matched_ids]
    temp[is.na(temp)] = 0
    final_bed1 = cbind(pooled_tabb[,c(1:3)], temp)
    write.table(final_bed1, paste0(output_cell, "/", output_bed, "_", tissuenames[numtt], ".bed"),
                sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}

Baseline_2024_EG_calc <- function(scores,
                             pred_type = "dist_to_tss",
                             preds_file="/data/deyk/ENCODE/Baseline_Preds_2024/preds.txt",
                             tissuenames = c("ALL", "BLD", "BRN", "K562", "GM12878", "Mono"),
                             output_cell,
                             output_bed = "temp.bed"){
  library(data.table)
  logreg_dnase_biosamples = read.table(preds_file, header = F)[,1]
  for(numtt in 1:length(tissuenames)){
    if(tissuenames[numtt] == "ALL"){
      keep_eids =  logreg_dnase_biosamples
    }else if (tissuenames[numtt] == "BLD"){
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
    }else if (tissuenames[numtt] == "BRN"){
      keep_eids = unique(logreg_dnase_biosamples[c(grep("_astro", logreg_dnase_biosamples),
                                                   grep("neur", logreg_dnase_biosamples),
                                                   grep("cerebell", logreg_dnase_biosamples),
                                                   grep("cortex", logreg_dnase_biosamples),
                                                   grep("brain", logreg_dnase_biosamples),
                                                   grep("glia", logreg_dnase_biosamples))])
    }else if (tissuenames[numtt] == "K562"){
      keep_eids = logreg_dnase_biosamples[grep("K562", logreg_dnase_biosamples)]
    }else if (tissuenames[numtt] == "Mono"){
      keep_eids = logreg_dnase_biosamples[grep("mono", logreg_dnase_biosamples)]
    }else if (tissuenames[numtt] == "GM12878"){
      keep_eids = logreg_dnase_biosamples[grep("GM12878", logreg_dnase_biosamples)]
    }

    pooled_tabb = c()
    for(numl in 1:length(keep_eids)){
      if(file.exists(paste0("/data/deyk/ENCODE/Baseline_Preds_2024/",
                            keep_eids[numl], "/", pred_type, ".tsv.gz"))){
        preds_tabb = data.frame(fread(paste0("/data/deyk/ENCODE/Baseline_Preds_2024/",
                                             keep_eids[numl], "/", pred_type, ".tsv.gz")))
        if(pred_type == "dist_to_tss"){
          preds_tabb = preds_tabb[which(preds_tabb$Score < 53976), ]
        }
        if(pred_type == "dist_to_gene"){
          preds_tabb = preds_tabb[which(preds_tabb$Score < 44840), ]
        }
        if(pred_type == "DHS_reads_by_dist_to_tss"){
          preds_tabb = preds_tabb[which(preds_tabb$Score > 3.149175673932778e-4), ]
        }
        if(pred_type == "DHS_reads_by_dist_to_tss_norm"){
          preds_tabb = preds_tabb[which(preds_tabb$Score > 7.391383676221246e-4), ]
        }
        if(pred_type == "H3K27ac_reads_by_dist_to_tss"){
          preds_tabb = preds_tabb[which(preds_tabb$Score > 2.2513369374164912e-4), ]
        }
        if(pred_type == "H3K27ac_reads_by_dist_to_tss_norm"){
          preds_tabb = preds_tabb[which(preds_tabb$Score > 8.201190258983458e-4), ]
        }
        if(pred_type == "nearest_expressed_gene"){
          preds_tabb = preds_tabb[which(preds_tabb$Score > 0), ]
        }
        if(pred_type == "nearest_expressed_tss"){
          preds_tabb = preds_tabb[which(preds_tabb$Score > 0), ]
        }
        tabb2 = cbind.data.frame(preds_tabb$X.chr, preds_tabb$start, preds_tabb$end, preds_tabb$TargetGene)
        colnames(tabb2) = c("chr", "start", "end", "TargetGene")
        pooled_tabb = rbind(pooled_tabb, tabb2)
      }
    }

    matched_ids = match(pooled_tabb$TargetGene, names(scores))
    temp = as.numeric(scores)[matched_ids]
    temp[is.na(temp)] = 0
    final_bed1 = cbind(pooled_tabb[,c(1:3)], temp)
    write.table(final_bed1, paste0(output_cell, "/", output_bed, "_", tissuenames[numtt], ".bed"),
                sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}

ENCODE_E2G_2024_ext_calc <- function(scores,
                        pred_cell="/data/deyk/ENCODE/ENCODE_E2G_2024/predictions_ext",
                        tissuenames= c("K562", "GM12878"),
                        output_cell,
                        output_bed = "temp.bed"){
  library(data.table)
  for(numtt in 1:length(tissuenames)){
    if(tissuenames[numtt] == "K562"){
      tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562", "/", "encode_e2g_predictions_threshold0.336.tsv.gz")))
    }else if(tissuenames[numtt] == "GM12878"){
      tabb_pre = data.frame(fread(paste0(pred_cell, "/", "GM12878", "/", "encode_e2g_predictions_threshold0.336.tsv.gz")))
    }
    tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,6])
    colnames(tabb2) = c("chr", "start", "end", "TargetGene")

    matched_ids = match(tabb2$TargetGene, names(scores))
    temp = as.numeric(scores)[matched_ids]
    temp[is.na(temp)] = 0
    final_bed1 = cbind(tabb2[,c(1:3)], temp)
    write.table(final_bed1, paste0(output_cell, "/", output_bed, "_", tissuenames[numtt], ".bed"),
                sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}


ABC_2024_ext_calc <- function(scores,
                             pred_cell="/data/deyk/ENCODE/ABC_2024/predictions_ext",
                             tissuenames= c("K562", "GM12878"),
                             output_cell,
                             output_bed = "temp.bed"){
  library(data.table)
  for(numtt in 1:length(tissuenames)){
    if(tissuenames[numtt] == "K562"){
      tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562/Predictions", "/", "ABC_pred_filtered_0.027_selfPromoters.tsv.gz")))
    }else if(tissuenames[numtt] == "GM12878"){
      tabb_pre = data.frame(fread(paste0(pred_cell, "/", "GM12878/Predictions/", "/", "ABC_pred_filtered_0.027_selfPromoters.tsv.gz")))
    }
    tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,11])
    colnames(tabb2) = c("chr", "start", "end", "TargetGene")

    matched_ids = match(tabb2$TargetGene, names(scores))
    temp = as.numeric(scores)[matched_ids]
    temp[is.na(temp)] = 0
    final_bed1 = cbind(tabb2[,c(1:3)], temp)
    write.table(final_bed1, paste0(output_cell, "/", output_bed, "_", tissuenames[numtt], ".bed"),
                sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
}


scE2G_calc <- function(scores,
                       pred_cell="/data/deyk/IGVF/IGVF_scE2G",
                       output_cell,
                       output_bed = "temp.bed"){

  ##############################. multiome-7features-ENCODE-E2G  ##################################################

  tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562_Xu_multiome_7features_genomewide_predictions.tsv.gz")))
  tabb_pre = tabb_pre[which(tabb_pre$ENCODE.rE2G.Score > 0.171), ]
  tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,6])
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")
  pooled_tabb = tabb2
  matched_ids = match(pooled_tabb$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(pooled_tabb[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed, "_", "multiome_7features_K562", ".bed"),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

  ##############################. multiome-7features-intx-ENCODE-E2G  ##################################################

  tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562_Xu_multiome_7features_intx_genomewide_predictions.tsv.gz")))
  tabb_pre = tabb_pre[which(tabb_pre$ENCODE.rE2G.Score > 0.176), ]
  tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,6])
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")
  pooled_tabb = tabb2
  matched_ids = match(pooled_tabb$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(pooled_tabb[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed, "_", "multiome_7features_intx_K562", ".bed"),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

  ##############################. scATAC-6features-ENCODE-E2G  ##################################################

  tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562_Xu_scATAC_6features_genomewide_predictions.tsv.gz")))
  tabb_pre = tabb_pre[which(tabb_pre$ENCODE.rE2G.Score > 0.154), ]
  tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,6])
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")
  pooled_tabb = tabb2
  matched_ids = match(pooled_tabb$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(pooled_tabb[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed, "_", "scATAC_6features_K562", ".bed"),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

  ##############################. scATAC-6features-intx-ENCODE-E2G  ##################################################

  tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562_Xu_scATAC_6features_intx_genomewide_predictions.tsv.gz")))
  tabb_pre = tabb_pre[which(tabb_pre$ENCODE.rE2G.Score > 0.22), ]
  tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,6])
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")
  pooled_tabb = tabb2
  matched_ids = match(pooled_tabb$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(pooled_tabb[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed, "_", "scATAC_6features_intx_K562", ".bed"),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

  ##############################. scATAC-ABC.  ##################################################

  tabb_pre = data.frame(fread(paste0(pred_cell, "/", "K562_Xu_scATAC_6features_genomewide_predictions.tsv.gz")))
  tabb_pre = tabb_pre[which(tabb_pre$ABC.Score > 0.015), ]
  tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,6])
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")
  pooled_tabb = tabb2
  matched_ids = match(pooled_tabb$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(pooled_tabb[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed, "_", "scATAC_ABC_K562", ".bed"),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

  ##################### Kendall correlations ##################################

  tabb_pre = data.frame(fread(paste0(pred_cell, "/", "Pairs.Kendall.tsv.gz")))
  tabb_pre = tabb_pre[which(tabb_pre$Kendall > 0.011), ]
  tabb2 = cbind.data.frame(tabb_pre[,1], tabb_pre[,2], tabb_pre[,3], tabb_pre[,4])
  colnames(tabb2) = c("chr", "start", "end", "TargetGene")
  pooled_tabb = tabb2
  matched_ids = match(pooled_tabb$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(pooled_tabb[,c(1:3)], temp)
  write.table(final_bed1, paste0(output_cell, "/", output_bed, "_", "Kendall_K562", ".bed"),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

}

########################     Perform GWAS gene prioritization using PoPS and enhancer-gene linking strategy      #################################

library(GenomicRanges)
library(dplyr)
library(data.table)

benchmark_EG_finemapUKBB_gene = function(EGpred_file,
                                      tissue = c("BLD", "ALL"),
                                      use_pops = FALSE,
                                      num_pops = 3,
                                      num_EGgenes = 1){

  eg_preds = data.frame(fread(EGpred_file))
  eg_preds = eg_preds[, c("chr", "start", "end", "TargetGene", "Score")]
 
  ###############################  generate all finemapped variant set   #######################################################

  credvarset = list()
  diseases = list.files("/data/deyk/GWAS/UKBiobank/finemapping/Jamboree_GWAS_Traits/191010_UKBB_SUSIE/CredibleSets/")
  for(numd in 1:length(diseases)){
    vv = read.table(paste0("/data/deyk/GWAS/UKBiobank/finemapping/Jamboree_GWAS_Traits/191010_UKBB_SUSIE/CredibleSets/",
                           diseases[numd], "/", "variant.list.txt"), header = T)
    credvarset[[numd]] = vv
    cat("We are processing disease:", numd, "\n")
  }
  names(credvarset) = diseases

  ###############################  generate Variant universe from BIM files   #######################################################

  bimlist = list()
  for(numchr in 1:22){
    xx = data.frame(fread(paste0("/data/deyk/kushal/extras/BIMS_hg38/", "1000G.EUR.QC.", numchr, ".bim")))
    colnames(xx) = c("CHR", "SNP", "CM", "BP", "REF", "ALT")
    bimlist[[numchr]] = xx
    cat("Processing bimfile for chr:", numchr, "\n")
  }
  names(bimlist) = paste0("chr", 1:22)

  ##########################  PoPS causal gene file   #########################################################

  pops_silver = data.frame(fread("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/UKBiobank.ABCGene.anyabc.tsv"))
  if(tissue == "ALL"){
    pops_silver = pops_silver
  }else if (tissue == "BLD"){
    traits = c("Eosino", "Lym", "Mono", "Neutro", "Plt", "RBC", "MCV")
    pops_silver = pops_silver[which(pops_silver$Disease %in% traits == T), ]
  }

  ######################  Overlap EG predictions on PoPS causal gene linked regions  ############################################

  gr1 = GRanges(seqnames = eg_preds$chr,
                ranges = IRanges(start=eg_preds$start, end = eg_preds$end))

  chrs = as.character(sapply(pops_silver$CredibleSet, function(x) return(strsplit(x, ":")[[1]][1])))
  cred_start =  as.numeric(sapply(pops_silver$CredibleSet, function(x) return(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1])))
  cred_end =  as.numeric(sapply(pops_silver$CredibleSet, function(x) return(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2])))

  pops_silver$chrs = chrs
  pops_silver$cred_start = cred_start
  pops_silver$cred_end = cred_end

  ucreds = unique(pops_silver$CredibleSet)

  true_positives = c()
  positives = c()
  relevant = c()
  creds_vec = c()

  for(numu in 1:length(ucreds)){
    temp = pops_silver[which(pops_silver$CredibleSet == ucreds[numu]), ]
    causal_gene = temp$TargetGene[which(temp$truth == 1)]
    disease_name = unique(temp$Disease)

    if(use_pops){
      pops_genes = temp$TargetGene[order(temp$POPS.Rank, decreasing = F)[1:num_pops]]
    }

    chr_cred = as.character(strsplit(ucreds[numu], ":")[[1]][1])
    rsids = c()
    pips = c()
    for(num_trait in 1:length(credvarset)){
      idx =  which(credvarset[[disease_name[num_trait]]]$CredibleSet == ucreds[numu])
      rsids = c(rsids, credvarset[[disease_name[num_trait]]]$rsid[idx])
      pips = c(pips, credvarset[[disease_name[num_trait]]]$pip[idx])
    }
    rsids = rsids[which(pips > 0.10)]
    coding_promoter_snps = read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Promoter_Coding_variants.txt", header=F)[,1]
    rsids = setdiff(rsids, coding_promoter_snps)
    
    if(length(intersect(rsids,bimlist[[chr_cred]]$SNP)) == 0){
      next
    }
    starts = bimlist[[chr_cred]]$BP[match(intersect(rsids,bimlist[[chr_cred]]$SNP), bimlist[[chr_cred]]$SNP)]-1
    ends = bimlist[[chr_cred]]$BP[match(intersect(rsids,bimlist[[chr_cred]]$SNP), bimlist[[chr_cred]]$SNP)]+1


    gr_cred = GRanges(seqnames = chr_cred,
                      ranges = IRanges(start = starts,
                                       end = ends))
    cc=findOverlaps(gr_cred, gr1, type = "any", select = "all")

    if(length(cc) > 0){
      any = c(any, 1)
      eg_preds3 = eg_preds[unique(subjectHits(cc)), c("TargetGene", "Score")]
      gene_scores = tapply(eg_preds3$Score, eg_preds3$TargetGene, max)
      gene_scores2 = rep(0, length(temp$TargetGene))
      names(gene_scores2) = temp$TargetGene
      gene_scores2[match(intersect(names(gene_scores), temp$TargetGene), temp$TargetGene)] =
        gene_scores[match(intersect(names(gene_scores), temp$TargetGene), names(gene_scores))]

      if(max(gene_scores2) == 0){
        true_positives = c(true_positives, 0)
        positives = c(positives, 0)
      }else{
        genes_selected = names(gene_scores2)[order(gene_scores2, decreasing = T)[1:pmin(num_EGgenes, length(gene_scores2))]]
        if(use_pops){
          genes_selected = intersect(genes_selected, pops_genes)
        }
        true_positives = c(true_positives, length(intersect(genes_selected, causal_gene)))
        positives = c(positives, length(genes_selected))
      }
    }else{
      any = c(any, 0)
      true_positives = c(true_positives, 0)
      positives = c(positives, 0)
    }
    relevant = c(relevant, length(causal_gene))
    cat("Processing credible set:", numu, "out of ", length(ucreds), "\n")
  }

  ###################3  Report Precision and Recall   ################################################################

  precision = sum(true_positives[which(positives !=  0)])/length(positives[which(positives !=  0)])
  recall = sum(true_positives[which(relevant !=  0)])/length(relevant[which(relevant !=  0)])

  precision2 = sum(true_positives[which(positives !=  0)])/sum(positives[which(positives !=  0)])
  recall2 = sum(true_positives[which(relevant !=  0)])/sum(relevant[which(relevant !=  0)])

  precision_boot = c()
  precision2_boot = c()
  recall_boot = c()
  recall2_boot = c()

  for(nboot in 1:100){
   idx = sample(1:length(positives), length(positives), replace = T)
   true_positives_boot = true_positives[idx]
   positives_boot = positives[idx]
   relevant_boot =- relevant[idx]

   precision_boot = c(precision_boot, sum(true_positives_boot[which(positives_boot !=  0)])/length(positives_boot[which(positives_boot !=  0)]))
   recall_boot = c(recall_boot, sum(true_positives_boot[which(relevant_boot !=  0)])/length(relevant_boot[which(relevant_boot !=  0)]))

   precision2_boot = c(precision2_boot, sum(true_positives_boot[which(positives_boot !=  0)])/sum(positives_boot[which(positives_boot !=  0)]))
   recall2_boot = c(recall2_boot, sum(true_positives_boot[which(relevant_boot !=  0)])/sum(relevant_boot[which(relevant_boot !=  0)]))
 }

 sd_precision = sd(precision_boot)
 sd_recall = sd(recall_boot)
 sd_precision2 = sd(precision2_boot)
 sd_recall2 = sd(recall2_boot)

  out_ll = list("precision" = precision,
                "sd.precision" = sd_precision,
                "recall" = recall,
                "sd.recall" = sd_recall,
                "precision2" = precision2,
                "sd.precision2" = sd_precision2,
                "recall2" = recall2,
                "sd.recall2" = sd_recall2,
                "tissue" = tissue,
                "use_pops" = use_pops,
                "num_pops" = num_pops)
  return(out_ll)
}







library(data.table)
library(R.utils)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
annot_cell = toString(args[1])
annot_name = toString(args[2])
output_cell = toString(args[3])

# annot_cell = "/n/groups/price/kushal/ENCODE/data/ANNOTATIONS_hg38/Primary_Freeze2022/ALL/ALL"
# annot_name = "Ramil_BLD"

finemap1 = read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/GWAS_traits/Other_finemapped/FINEMAP_Ulirsch_ALL_WBC_10")[,1]
finemap2 = read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/GWAS_traits/Other_finemapped/FINEMAP_Ulirsch_ALL_WBC_50")[,1]
finemap3 = read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/GWAS_traits/Other_finemapped/FINEMAP_Ulirsch_ALL_Mono_10")[,1]
finemap4 = read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/GWAS_traits/Other_finemapped/FINEMAP_Ulirsch_ALL_Mono_50")[,1]
finemap5 = read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/GWAS_traits/Other_finemapped/FINEMAP_Ulirsch_ALL_Lym_10")[,1]
finemap6 = read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/GWAS_traits/Other_finemapped/FINEMAP_Ulirsch_ALL_Lym_50")[,1]
finemap7= read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/GWAS_traits/Other_finemapped/FINEMAP_Ulirsch_ALL_Eosino_10")[,1]
finemap8 = read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/GWAS_traits/Other_finemapped/FINEMAP_Ulirsch_ALL_Eosino_50")[,1]
finemap9= read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/GWAS_traits/Other_finemapped/FINEMAP_Ulirsch_ALL_AID_Combined_10")[,1]
finemap10 = read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/GWAS_traits/Other_finemapped/FINEMAP_Ulirsch_ALL_AID_Combined_50")[,1]
finemap11= read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/GWAS_traits/Other_finemapped/FINEMAP_Ulirsch_ALL_IBD_10")[,1]
finemap12 = read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/GWAS_traits/Other_finemapped/FINEMAP_Ulirsch_ALL_IBD_50")[,1]

coding_promoter_snps = read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Promoter_Coding_variants.txt", header=F)[,1]

#finemap1 = setdiff(finemap1, coding_promoter_snps)
#finemap2 = setdiff(finemap2, coding_promoter_snps)
#finemap3 = setdiff(finemap3, coding_promoter_snps)
#finemap4 = setdiff(finemap4, coding_promoter_snps)
#finemap5 = setdiff(finemap5, coding_promoter_snps)
#finemap6 = setdiff(finemap6, coding_promoter_snps)
#finemap7 = setdiff(finemap7, coding_promoter_snps)
#finemap8 = setdiff(finemap8, coding_promoter_snps)
#finemap9 = setdiff(finemap9, coding_promoter_snps)
#finemap10 = setdiff(finemap10, coding_promoter_snps)
#finemap11 = setdiff(finemap11, coding_promoter_snps)
#finemap12 = setdiff(finemap12, coding_promoter_snps)

sum_farh_chr = rep(0, 12)
sum_finemap = rep(0, 12)
sum_binannot = 0
sum_all = 0
sum_aggr_annot = matrix(0, 100, 12)
per_chr_1 = matrix(0, 22, 12)
per_chr_2 = array(0, 22)
per_chr_3 = matrix(0, 22, 12)


for(numchr in 1:22){
  df1 = data.frame(fread(paste0(annot_cell, "/",
                                annot_name, "/",
                                annot_name, ".", numchr, ".annot.gz")))

  finemap1_annot = rep(0, nrow(df1))
  finemap2_annot = rep(0, nrow(df1))
  finemap3_annot = rep(0, nrow(df1))
  finemap4_annot = rep(0, nrow(df1))
  finemap5_annot = rep(0, nrow(df1))
  finemap6_annot = rep(0, nrow(df1))
  finemap7_annot = rep(0, nrow(df1))
  finemap8_annot = rep(0, nrow(df1))
  finemap9_annot = rep(0, nrow(df1))
  finemap10_annot = rep(0, nrow(df1))
  finemap11_annot = rep(0, nrow(df1))
  finemap12_annot = rep(0, nrow(df1))

  finemap1_annot[match(intersect(df1$SNP, finemap1), df1$SNP)] = 1
  finemap2_annot[match(intersect(df1$SNP, finemap2), df1$SNP)] = 1
  finemap3_annot[match(intersect(df1$SNP, finemap3), df1$SNP)] = 1
  finemap4_annot[match(intersect(df1$SNP, finemap4), df1$SNP)] = 1
  finemap5_annot[match(intersect(df1$SNP, finemap5), df1$SNP)] = 1
  finemap6_annot[match(intersect(df1$SNP, finemap6), df1$SNP)] = 1
  finemap7_annot[match(intersect(df1$SNP, finemap7), df1$SNP)] = 1
  finemap8_annot[match(intersect(df1$SNP, finemap8), df1$SNP)] = 1
  finemap9_annot[match(intersect(df1$SNP, finemap9), df1$SNP)] = 1
  finemap10_annot[match(intersect(df1$SNP, finemap10), df1$SNP)] = 1
  finemap11_annot[match(intersect(df1$SNP, finemap11), df1$SNP)] = 1
  finemap12_annot[match(intersect(df1$SNP, finemap12), df1$SNP)] = 1

  binannot = df1[,5]
#  coding_promoter_snps = read.table("/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Promoter_Coding_variants.txt", header=F)[,1]
# binannot[match(intersect(coding_promoter_snps, df1$SNP), df1$SNP)] = 0

  per_chr_1[numchr, ] = c(sum(binannot[finemap1_annot == 1]),
                          sum(binannot[finemap2_annot == 1]),
                          sum(binannot[finemap3_annot == 1]),
                          sum(binannot[finemap4_annot == 1]),
                          sum(binannot[finemap5_annot == 1]),
                          sum(binannot[finemap6_annot == 1]),
                          sum(binannot[finemap7_annot == 1]),
                          sum(binannot[finemap8_annot == 1]),
                          sum(binannot[finemap9_annot == 1]),
                          sum(binannot[finemap10_annot == 1]),
                          sum(binannot[finemap11_annot == 1]),
                          sum(binannot[finemap12_annot == 1]))

  sum_farh_chr = sum_farh_chr + per_chr_1[numchr, ]
  per_chr_2[numchr] = sum(binannot)
  sum_binannot = sum_binannot + per_chr_2[numchr]
  sum_all = sum_all + nrow(df1)

  per_chr_3[numchr, ] = c(sum(finemap1_annot),
                          sum(finemap2_annot),
                          sum(finemap3_annot),
                          sum(finemap4_annot),
                          sum(finemap5_annot),
                          sum(finemap6_annot),
                          sum(finemap7_annot),
                          sum(finemap8_annot),
                          sum(finemap9_annot),
                          sum(finemap10_annot),
                          sum(finemap11_annot),
                          sum(finemap12_annot))
  sum_finemap = sum_finemap + per_chr_3[numchr, ]

  aggr_annot = c()

  for(nboot in 1:100){
    binannot_samp = sample(binannot, length(binannot), replace = F)
    aggr_annot = rbind(aggr_annot, c(sum(binannot_samp[finemap1_annot == 1]),
                                     sum(binannot_samp[finemap2_annot == 1]),
                                     sum(binannot_samp[finemap3_annot == 1]),
                                     sum(binannot_samp[finemap4_annot == 1]),
                                     sum(binannot_samp[finemap5_annot == 1]),
                                     sum(binannot_samp[finemap6_annot == 1]),
                                     sum(binannot_samp[finemap7_annot == 1]),
                                     sum(binannot_samp[finemap8_annot == 1]),
                                     sum(binannot_samp[finemap9_annot == 1]),
                                     sum(binannot_samp[finemap10_annot == 1]),
                                     sum(binannot_samp[finemap11_annot == 1]),
                                     sum(binannot_samp[finemap12_annot == 1])))
  }
  sum_aggr_annot = sum_aggr_annot + aggr_annot
  cat("We are at chr:", numchr, "\n")
}

enr_bootmat=matrix(0, 100, 12)
prec_bootmat=matrix(0, 100, 12)
recall_bootmat=matrix(0, 100, 12)

for(nboot in 1:100){
  idx = sample(1:22, 22, replace=TRUE)
  enr_bootmat[nboot, ] = (colSums(per_chr_1[idx, ])/colSums(per_chr_3[idx, ]))/(sum(per_chr_2[idx])/sum_all)
  prec_bootmat[nboot, ] = (colSums(per_chr_1[idx, ])/sum(per_chr_2[idx]))
  recall_bootmat[nboot, ] = (colSums(per_chr_1[idx, ])/colSums(per_chr_3[idx, ]))
}

enr = (sum_farh_chr/sum_binannot)/(sum_finemap/sum_all)
precision = (sum_farh_chr/sum_binannot)
recall = (sum_farh_chr/sum_finemap)

senr = apply(enr_bootmat, 2, sd)
sprec = apply(prec_bootmat, 2, sd)
srecall = apply(recall_bootmat, 2, sd)

boot_enr = matrix(0, 100, 12)
boot_prec = matrix(0, 100, 12)
boot_recall = matrix(0, 100, 12)

for(mm  in 1:100){
  boot_enr[mm, ] = (sum_aggr_annot[mm, ]/sum_binannot)/(sum_finemap/sum_all)
  boot_prec[mm, ] = (sum_aggr_annot[mm, ]/sum_binannot)
  boot_recall[mm, ] = (sum_aggr_annot[mm, ]/sum_finemap)
}

penr = c()
pprec = c()
precall = c()
for(mm in 1:12){
  penr = c(penr, pnorm(enr[mm], mean(boot_enr[,mm]), sd(boot_enr[,mm]), lower.tail = F))
  pprec = c(pprec, pnorm(precision[mm], mean(boot_prec[,mm]), sd(boot_prec[,mm]), lower.tail = F))
  precall = c(precall, pnorm(recall[mm], mean(boot_recall[,mm]), sd(boot_recall[,mm]), lower.tail = F))
}

outdf = cbind(enr, senr, penr, precision, sprec, pprec, recall, srecall, precall)
colnames(outdf) = c("ENR", "sENR", "pENR", "Precision", "s.Precision", "p.Precision", "Recall", "s.Recall", "p.Recall")
rownames(outdf) = c("Ulirsch_ALL_WBC_10", "Ulirsch_ALL_WBC_50",
                    "Ulirsch_ALL_Mono_10", "Ulirsch_ALL_Mono_50",
                    "Ulirsch_ALL_Lym_10", "Ulirsch_ALL_Lym_50",
                    "Ulirsch_ALL_Eosino_10", "Ulirsch_ALL_Eosino_50",
                    "Ulirsch_ALL_AID_Combined_10", "Ulirsch_ALL_AID_Combined_50",
                    "Ulirsch_ALL_IBD_10", "Ulirsch_ALL_IBD_50")

write.table(outdf, file = paste0(output_cell, "/", annot_name, ".FINEMAP.RES.txt"),
            col.names = T, row.names = T, sep = "\t", quote=F)


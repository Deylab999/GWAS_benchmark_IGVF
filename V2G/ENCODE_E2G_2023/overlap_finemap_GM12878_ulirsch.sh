annot_cell=/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/ANNOTATIONS_hg38/ALL/ALL
output_cell=/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/FINEMAP_ENR/Primary_Freeze2022_Ulirsch_GM12878

module load gcc/10.2.0
module load R

IFS="
"

TASKFILE=/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Primary_EG_GM12878.txt
                                                                                               
for line in `cat $TASKFILE | awk '{print $1}' | sort | uniq`;
do
    annot_name=`echo $line | awk '{print $1}'`
    echo $annot_name
    if [ ! -f $output_cell/$annot_name.txt ]
    then
    cmd="Rscript overlap_finemap_GM12878_ulirsch.R  $annot_cell $annot_name $output_cell"
    bsub -W 100 -R "rusage[mem=35]" -e chunks.err -o chunks.out -n 1 "$cmd"
    fi
done


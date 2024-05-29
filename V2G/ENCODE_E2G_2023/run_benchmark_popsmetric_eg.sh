annot_cell=/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Blood_EG_predictions
output_cell=/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/FINEMAP_POPS_GENE3

module load gcc/10.2.0
module load R

IFS="
"

TASKFILE=/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/finemap_precision_recall.txt
                                                                                               
for line in `cat $TASKFILE`;
do
    annot_name=`echo $line | awk '{print $1}'`
    tissue=`echo $line | awk '{print $2}'`
    output_file=`echo $line | awk '{print $3}'`
    echo $annot_name $tissue $output_file
#    if [ ! -f $output_cell/$output_file ]
#    then
    cmd="Rscript run_benchmark_popsmetric_eg.R  $annot_cell/$annot_name $tissue $output_cell/$output_file"
    bsub -W 300 -R "rusage[mem=100]" -e chunks.err -o chunks.out -n 1 "$cmd"
#    fi
done


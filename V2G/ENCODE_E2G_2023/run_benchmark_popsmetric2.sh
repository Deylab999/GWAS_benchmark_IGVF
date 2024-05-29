annot_cell=~/temp/
output_cell=~/output/
module load gcc/6.2.0
module load R

IFS="
"

TASKFILE=finemap_precision_recall3.txt
                                                                                               
for line in `cat $TASKFILE`;
do
    annot_name=`echo $line | awk '{print $1}'`
    tissue=`echo $line | awk '{print $2}'`
    output_file=`echo $line | awk '{print $3}'`
    echo $annot_name $tissue $output_file
#    if [ ! -f $output_cell/$output_file ]
#    then
    cmd="Rscript run_benchmark_popsmetric2.R  $annot_cell/$annot_name $tissue $output_cell/$output_file"
    sbatch --time=200:00 --mem=120000 --output=finemap.out --error=finemap.err -p short -c 1 --wrap="$cmd"
#    fi
done


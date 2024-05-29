annot_cell=~/temp/
output_cell=~/output/

module load gcc/6.2.0
module load R

IFS="
"

TASKFILE=task.txt
                                                                                               
for line in `cat $TASKFILE | awk '{print $1}' | sort | uniq`;
do
    annot_name=`echo $line | awk '{print $1}'`
    echo $annot_name
    if [ ! -f $output_cell/$annot_name.txt ]
    then
    cmd="Rscript overlap_finemap_bloodbio_ulirsch.R  $annot_cell $annot_name $output_cell"
    sbatch --time=90:00 --mem=20000 --output=finemap.out --error=finemap.err -p short -c 1 --wrap="$cmd"
    fi
done


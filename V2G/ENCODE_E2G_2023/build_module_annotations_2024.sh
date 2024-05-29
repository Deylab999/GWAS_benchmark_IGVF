genescore_cell=/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Gene_Scores
bed_cell=/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/BEDFILES

module load gcc/10.2.0
module load R

IFS="
"

TASKFILE=/data/deyk/kushal/ENCODE_E2G_GWAS_benchmark/data/Tnames.txt  ## ALL

for line in `cat $TASKFILE | awk '{print $1}' | sort | uniq`;
do
   temp=`echo $line | awk '{print $1}'`
   bed_dir=$bed_cell/$temp
   if [ ! -d $bed_dir ]
   then
       mkdir $bed_dir
   fi
   genescore_dir=$genescore_cell/$temp
   if [ ! -d $genescore_dir ]
   then
       mkdir $genescore_dir
   fi
   for ll in `ls -1 $genescore_dir | sed 's/\.txt//g' | awk '{print $1}' | sort | uniq`;
   do
      annot_name=`echo $ll | awk '{print $1}'`
      echo $temp $annot_name
#     if [ ! -f $bed_dir/$annot_name/100kb.bed ]
#      then
	  cmd="Rscript build_module_annotations_2024.R  $genescore_dir $bed_dir $annot_name"
          bsub -W 600 -R "rusage[mem=35]" -e chunks.err -o chunks.out -n 1 "$cmd"
#      fi
   done
done


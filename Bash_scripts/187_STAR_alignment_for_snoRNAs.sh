#!/bin/bash

eval "$(conda shell.bash hook)"
 

MASTER_ROUTE=$1
analysis=$2
reference_fasta=$3
reference_features=$4
feature=$5

trimmed_dir=$(echo "/scratch/manuel.tardaguila/my_own_RNA_pipeline/new_trimming/trimmed_sequences/")


module load STAR/2.7.10a

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

mkdir -p $output_dir

Log_files=$(echo "$output_dir""/""Log_files/")


mkdir -p $Log_files


counts_dir=$(echo "/group/soranzo/manuel.tardaguila/small_rna/""$feature""/")

mkdir -p $counts_dir


sample_array=$(echo 'MCO_01334,MCO_01335,MCO_01338,MCO_01339,MCO_01342,MCO_01343')
equivalence_array=$(echo 'S1,S2,S3,S4,S5,S6')
genotype_array=$(echo 'wt_1,rs62237617_1,wt_2,rs62237617_2,wt_3,rs62237617_3')


a=($(echo "$sample_array" | tr "," '\n'))
b=($(echo "$equivalence_array" | tr "," '\n'))
c=($(echo "$genotype_array" | tr "," '\n'))

declare -a arr

array_2_length=${#a[@]}

for (( i=0; i<${array_2_length}; i=i+1 ));
do

    sample_array_sel=${a[$i]}
    echo "$sample_array_sel"

    prefix=$(echo "$output_dir""$sample_array_sel""_")


    equivalence_sel=${b[$i]}

    echo "$equivalence_sel"

    genotype_sel=${c[$i]}

    echo "$genotype_sel"

    raw_sample=$(echo "$sample_array_sel""_smLIB_""$equivalence_sel""_""L004")

    echo "1 $raw_sample"

    r1_TRIMMED=$(echo "$trimmed_dir""$raw_sample""_R1_trimmed.fastq.gz")


    ################## STAR alingment  #######################################################################################

     type=$(echo "$sample_array_sel""_STAR_alingment")
     outfile_STAR_alingment=$(echo "$Log_files""outfile_1_""$type"".log")
     touch $outfile_STAR_alingment
     echo -n "" > $outfile_STAR_alingment
     name_STAR_alingment=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")

     processors=$(echo '8')
     mem=$(echo '4096')
     total_memory=$(echo "scale=0; ($mem / 1) * $processors" | bc)
     total_memory=$(echo "$total_memory""M")

     echo "$total_memory"

    myjobid_STAR_alingment=$(sbatch --job-name=$name_STAR_alingment --output=$outfile_STAR_alingment --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="STAR  \
    							    				 		  		  	    				  		      --runMode alignReads \
                                                                                                                                                                                      --genomeDir $reference_fasta \
                                                                                                                                                                                      --readFilesIn $r1_TRIMMED \
    																						      --runThreadN $processors \
    																						      --outFileNamePrefix $prefix \
    																						       --alignIntronMax 1 \
    																						       --alignEndsType EndToEnd \
    																						       --outFilterMultimapNmax 20 \
    																						       --outFilterMismatchNmax 1 \
    																							--outFilterScoreMinOverLread 0.9 \
    																						      --outFilterMatchNminOverLread 0.9 \
    																						      --outSAMtype BAM SortedByCoordinate \
    																						      --outBAMcompression 10 \
    																						      --readFilesCommand zcat")

    myjobid_seff_STAR_alingment=$(sbatch --dependency=afterany:$myjobid_STAR_alingment --open-mode=append --output=$outfile_STAR_alingment --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_STAR_alingment >> $outfile_STAR_alingment")


     ################## featureCounts_step  #######################################################################################

     type=$(echo "$sample_array_sel""_featureCounts_step")
     outfile_featureCounts_step=$(echo "$Log_files""outfile_2_""$type"".log")
     touch $outfile_featureCounts_step
     echo -n "" > $outfile_featureCounts_step
     name_featureCounts_step=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")

     processors=$(echo '2')
     mem=$(echo '1024')
     total_memory=$(echo "scale=0; ($mem / 1) * $processors" | bc)
     total_memory=$(echo "$total_memory""M")

     echo "$total_memory"

     output_counts=$(echo "$counts_dir""$sample_array_sel""_raw_counts.txt")
     input_bam=$(echo "$prefix""Aligned.sortedByCoord.out.bam")


     # --dependency=afterany:$myjobid_STAR_alingment

     conda activate /home/manuel.tardaguila/conda_envs/featureCounts

    myjobid_featureCounts_step=$(sbatch --dependency=afterany:$myjobid_STAR_alingment --job-name=$name_featureCounts_step --output=$outfile_featureCounts_step --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="featureCounts \
    									    					 		  		  	    			       -F SAF \
																						       -a  $reference_features \
																						       -o $output_counts \
																						      $input_bam")

    myjobid_seff_featureCounts_step=$(sbatch --dependency=afterany:$myjobid_featureCounts_step --open-mode=append --output=$outfile_featureCounts_step --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_featureCounts_step >> $outfile_featureCounts_step")

    conda deactivate



done



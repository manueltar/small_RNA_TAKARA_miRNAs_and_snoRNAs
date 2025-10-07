#!/bin/bash

eval "$(conda shell.bash hook)"
 

MASTER_ROUTE=$1
analysis=$2

path_fastq=$3

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")


mkdir -p $Log_files


trimmed_dir=$(echo "$output_dir""trimmed_sequences/")


mkdir -p $trimmed_dir

conda activate /home/manuel.tardaguila/conda_envs/cutadapt

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


    equivalence_sel=${b[$i]}

    echo "$equivalence_sel"

    genotype_sel=${c[$i]}

    echo "$genotype_sel"

    raw_sample=$(echo "$sample_array_sel""_smLIB_""$equivalence_sel""_""L004")

    echo "1 $raw_sample"

    r1=$(echo "$path_fastq""$raw_sample""_R1_001.fastq.gz")


    r1_TRIMMED=$(echo "$trimmed_dir""$raw_sample""_R1_trimmed.fastq.gz")


     ################################################ TRIM R1 to eliminate 1) sequences without polyA, 2) everything downstream of polyA and 3) the first three bases (from TSO, https://www.takarabio.com/documents/User%20Manual/SMARTer%20smRNA/SMARTer%20smRNA-Seq%20Kit%20for%20Illumina%20User%20Manual.pdf)

     type=$(echo "$sample_array_sel""_trimmed_R1")
     outfile_trimmed_R1=$(echo "$Log_files""outfile_1_""$type"".log")
     touch $outfile_trimmed_R1
     echo -n "" > $outfile_trimmed_R1
     name_trimmed_R1=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")

     processors=$(echo '4')
     mem=$(echo '4096')
     total_memory=$(echo "scale=0; ($mem / 1) * $processors" | bc)
     total_memory=$(echo "$total_memory""M")

     module unload fastqc/0.11.9

    myjobid_trimmed_R1=$(sbatch --job-name=$name_trimmed_R1 --output=$outfile_trimmed_R1 --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="cutadapt \
    							    				 		  		  	    				  		      -j $processors \
                                                                                                                                                                                      -m 15 \
                                                                                                                                                                                      -u 3 \
                                                                                                                                                                                      -a AAAAAAAAAA \
                                                                                                                                                                                      --discard-untrimmed \
                                                                                                                                                                                      -o $r1_TRIMMED \
    																						      $r1")
    myjobid_seff_trimmed_R1=$(sbatch --dependency=afterany:$myjobid_trimmed_R1 --open-mode=append --output=$outfile_trimmed_R1 --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_trimmed_R1 >> $outfile_trimmed_R1")

     ################################################ FASTQC to check trimming -------------------------------------------------------

     type=$(echo "$sample_array_sel""_FASTQC_trimmed_R1")
     outfile_FASTQC_trimmed_R1=$(echo "$Log_files""outfile_2_""$type"".log")
     touch $outfile_FASTQC_trimmed_R1
     echo -n "" > $outfile_FASTQC_trimmed_R1
     name_FASTQC_trimmed_R1=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")

     processors=$(echo '4')
     mem=$(echo '4096')
     total_memory=$(echo "scale=0; ($mem / 1) * $processors" | bc)
     total_memory=$(echo "$total_memory""M")

     
     module unload openjdk
     module load fastqc

     # --dependency=afterany:$myjobid_trimmed_R1

    myjobid_FASTQC_trimmed_R1=$(sbatch --dependency=afterany:$myjobid_trimmed_R1 --job-name=$name_FASTQC_trimmed_R1 --output=$outfile_FASTQC_trimmed_R1 --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="fastqc --outdir $trimmed_dir $r1_TRIMMED")
    myjobid_seff_FASTQC_trimmed_R1=$(sbatch --dependency=afterany:$myjobid_FASTQC_trimmed_R1 --open-mode=append --output=$outfile_FASTQC_trimmed_R1 --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_FASTQC_trimmed_R1 >> $outfile_FASTQC_trimmed_R1")




done



#!/bin/bash
#SBATCH --job-name=snoRNAs_STAR_indexing
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.tardaguila@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --output=snoRNAs_STAR_indexing_%j.log
#SBATCH --mem=32G
#SBATCH --time=36:00:00

module load STAR

echo "========================"
echo "Initiated: $(date)"

output_dir=$1
reference_genome=$2

STAR --runMode genomeGenerate \
     --genomeDir $output_dir \
     --genomeFastaFiles $reference_genome \
     --runThreadN 8 \
     --genomeSAindexNbases 11



echo "========================"
echo "Completed: $(date)"




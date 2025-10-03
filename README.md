# Consult TAKARA guidelines for processing the squencinf results  https://www.takarabio.com/documents/User%20Manual/SMARTer%20smRNA/SMARTer%20smRNA-Seq%20Kit%20for%20Illumina%20User%20Manual.pdf

# Download the necessary files

## download mature miRNA sequence and keep the Human miRNAs

$ wget https://www.mirbase.org/download/mature.fa
$ grep -A 1 "^>hsa-" mature.fa > human_mature_miRNAs.fa

## git clone the nf-core/smrnaseq

$ git clone https://github.com/nf-core/smrnaseq.git

# Run the script to trimm the raw reads (single end)

$ bash ~/Scripts/Wraper_scripts/186_trimmomatic_for_smallRNAs_from_TAKARA.sh /scratch/manuel.tardaguila/my_own_RNA_pipeline/ new_trimming /group/soranzo/paola.benaglio/small_rna/250818_A01481_0333_AHGVYHDSXF/fastq_raw/

# Run the nf-core/smrnaseq skipping the trimming (fastp) step

$ cd smrnaseq/

$ module load nextflow/24.10.4

$ module load singularity/3.8.5

$ nextflow run main.nf   -profile humantechnopole,singularity  -w /scratch/manuel.tardaguila/nf-work   --input /group/soranzo/manuel.tardaguila/small_rna/new_trimming/input.csv --mirtrace_species hsa  --outdir /group/soranzo/manuel.tardaguila/small_rna/new_trimming/  --skip_fastp --mature /scratch/manuel.tardaguila/my_own_RNA_pipeline/reference_sequences/human_mature_miRNAs.fa

# Jupyter notebook to find DE results in the edgeR mature miRNA counts file

------------------> DE_miRNA.ipynb



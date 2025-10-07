# 1. First download the Gencode reference

     wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz
      wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.tRNAs.gtf.gz

## 1.1 Remove the 'chr' from gencode 

$ awk 'BEGIN {OFS="\t"} {
    # Check if the first field starts with "chr"
    if ($1 ~ /^chr/) {
        # If it does, substitute "chr" with an empty string, starting at the beginning
        $1 = substr($1, 4)
    }
    # Print the modified line
    print
}' gencode.v49.primary_assembly.annotation.gtf > gencode.v49.primary_assembly.annotation_no_chr.gtf


## 1.2

$ awk 'BEGIN {OFS="\t"} {
    # Check if the first field starts with "chr"
    if ($1 ~ /^chr/) {
        # If it does, substitute "chr" with an empty string, starting at the beginning
        $1 = substr($1, 4)
    }
    # Print the modified line
    print
}' gencode.v49.tRNAs.gtf > gencode.v49.tRNAs_no_chr.gtf



# 2. Derive all the rRNA decoy files

$ awk 'BEGIN {OFS="\t"}
     {
         if ($3 ~ /^(gene|transcript|exon)$/ && $12 ~ /rRNA/) {
             print
         }
     }' /scratch/manuel.tardaguila/gencode.v49.primary_assembly.annotation_no_chr.gtf > rRNA_only.gtf


$ awk -F'\t' '
BEGIN {
    OFS="\t";
    print "GeneID", "Chr", "Start", "End", "Strand";
}
{
    # Check for valid feature line
    if ($1 !~ /^#/ && NF >= 14) {

        gene_name = $14;

        # Remove the leading/trailing quotes (") and the trailing semicolon (;)
        gsub(/"/, "", gene_name);
        gsub(/;/, "", gene_name);

        # Print the required fields: Gene Name (from $14), Chr ($1), Start ($4), End ($5), Strand ($7)
        print gene_name, $1, $4, $5, $7;
    }
}' rRNA_only.gtf > rRNA_converted.tsv

$ awk '
BEGIN {
    OFS="\t"
    print "#chrom", "chromStart", "chromEnd", "name", "score", "strand"
}
NR > 1 {
    # Input Fields: $1=GeneID, $2=Chr, $3=Start(1-based), $4=End(1-based), $5=Strand

    # BED requires 0-based start, so subtract 1 from the input Start ($3)
    BED_start = $3 - 1

    # Print: Chr($2), Start(0-based), End($4), Name($1), Score(0), Strand($5)
    print $2, BED_start, $4, $1, 0, $5
}' rRNA_converted.tsv  > rRNA_converted.bed


$ module load bedtools
$ bedtools getfasta -fi /scratch/manuel.tardaguila/GRCH38_ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed rRNA_converted.bed -s -fo rRNA_converted.fasta

# 2. Derive all the tRNA decoy files


$ echo -e "GeneID\tChr\tStart\tEnd\tStrand" > tRNA_converted.tsv


$ awk 'BEGIN {OFS="\t"} /^#/ {next} {gene_id = $14;chrom = $1;start = $4;end = $5;strand = $7; gsub(/\";$/, "", gene_id);gsub(/^\"/, "", gene_id) ; print gene_id, chrom, start, end, strand}' gencode.v49.tRNAs_no_chr.gtf >> tRNA_converted.tsv


$ awk '
BEGIN {
    OFS="\t"
    print "#chrom", "chromStart", "chromEnd", "name", "score", "strand"
}
NR > 1 {
    # Input Fields: $1=GeneID, $2=Chr, $3=Start(1-based), $4=End(1-based), $5=Strand

    # BED requires 0-based start, so subtract 1 from the input Start ($3)
    BED_start = $3 - 1

    # Print: Chr($2), Start(0-based), End($4), Name($1), Score(0), Strand($5)
    print $2, BED_start, $4, $1, 0, $5
}' tRNA_converted.tsv  > tRNA_converted.bed


$ module load bedtools
$ bedtools getfasta -fi /scratch/manuel.tardaguila/GRCH38_ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed tRNA_converted.bed -s -fo tRNA_converted.fasta
     


# 4. Obtain the snoRNA reference, convert to the correct format, get a bed file and finally a fasta file

$ awk 'BEGIN {OFS="\t"}
    {
        if ($3 ~ /^(gene|transcript|exon)$/ && $12 ~ /"snoRNA"/) {
            print
        }
    }' gencode.v49.primary_assembly.annotation_no_chr.gtf > SNORD_only.gtf

$ awk -F'\t' '
BEGIN {
    OFS="\t";
    print "GeneID", "Chr", "Start", "End", "Strand";
}
{
    # Check for valid feature line
    if ($1 !~ /^#/ && NF >= 14) {

        gene_name = $14;

        # Remove the leading/trailing quotes (") and the trailing semicolon (;)
        gsub(/"/, "", gene_name);
        gsub(/;/, "", gene_name);

        # Print the required fields: Gene Name (from $14), Chr ($1), Start ($4), End ($5), Strand ($7)
        print gene_name, $1, $4, $5, $7;
    }
}' SNORD_only.gtf > snord_converted.tsv


$ awk '
BEGIN {
    OFS="\t"
    print "#chrom", "chromStart", "chromEnd", "name", "score", "strand"
}
NR > 1 {
    # Input Fields: $1=GeneID, $2=Chr, $3=Start(1-based), $4=End(1-based), $5=Strand
    
    # BED requires 0-based start, so subtract 1 from the input Start ($3)
    BED_start = $3 - 1
    
    # Print: Chr($2), Start(0-based), End($4), Name($1), Score(0), Strand($5)
    print $2, BED_start, $4, $1, 0, $5
}' snord_converted.tsv > snord_converted.bed

$ module load bedtools

$ bedtools getfasta -fi /scratch/manuel.tardaguila/GRCH38_ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed snord_converted.bed -s -fo snord_converted.fasta



# 5. Obtain the miRNA reference, convert to the correct format, get a bed file and finally a fasta file

$ awk 'BEGIN {OFS="\t"}
     {
         if ($3 ~ /^(gene|transcript|exon)$/ && $12 ~ /"miRNA"/) {
             print
         }
     }' /scratch/manuel.tardaguila/gencode.v49.primary_assembly.annotation_no_chr.gtf > miRNA_only.gtf


$ awk -F'\t' '
BEGIN {
    OFS="\t";
    print "GeneID", "Chr", "Start", "End", "Strand";
}
{
    # Check for valid feature line
    if ($1 !~ /^#/ && NF >= 14) {

        gene_name = $14;

        # Remove the leading/trailing quotes (") and the trailing semicolon (;)
        gsub(/"/, "", gene_name);
        gsub(/;/, "", gene_name);

        # Print the required fields: Gene Name (from $14), Chr ($1), Start ($4), End ($5), Strand ($7)
        print gene_name, $1, $4, $5, $7;
    }
}' miRNA_only.gtf > miRNA_converted.tsv


$ awk '
BEGIN {
    OFS="\t"
    print "#chrom", "chromStart", "chromEnd", "name", "score", "strand"
}
NR > 1 {
    # Input Fields: $1=GeneID, $2=Chr, $3=Start(1-based), $4=End(1-based), $5=Strand
    
    # BED requires 0-based start, so subtract 1 from the input Start ($3)
    BED_start = $3 - 1
    
    # Print: Chr($2), Start(0-based), End($4), Name($1), Score(0), Strand($5)
    print $2, BED_start, $4, $1, 0, $5
}' miRNA_converted.tsv  > miRNA_converted.bed

$ module load bedtools

$ bedtools getfasta -fi /scratch/manuel.tardaguila/GRCH38_ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed miRNA_converted.bed -s -fo miRNA_converted.fasta

# 6. Create the decoy references:

  # 6.1 for miRNAs 
  
  $ cat Homo_sapiens.GRCh38.dna.primary_assembly.fa /group/soranzo/manuel.tardaguila/small_rna/Dependencies/rRNA_converted.fasta /group/soranzo/manuel.tardaguila/small_rna/Dependencies/tRNA_converted.fasta /group/soranzo/manuel.tardaguila/small_rna/Dependencies/snord_converted.fasta > /scratch/manuel.tardaguila/GRCH38_ensembl/combined_decoy_reference_miRNAs.fa


  # 6.2 for snoRNAs

  $ cat Homo_sapiens.GRCh38.dna.primary_assembly.fa /group/soranzo/manuel.tardaguila/small_rna/Dependencies/rRNA_converted.fasta /group/soranzo/manuel.tardaguila/small_rna/Dependencies/tRNA_converted.fasta /group/soranzo/manuel.tardaguila/small_rna/Dependencies/miRNA_converted.fasta > /scratch/manuel.tardaguila/GRCH38_ensembl/combined_decoy_reference_snoRNAs.fa
  


# 7. STAR indexing for small reads

  7.1 for miRNAs

  $ sbatch ~/Scripts/sbatch/12_snoRNAs_STAR_indexing.sh /scratch/manuel.tardaguila/STAR_indexed_genome_for_miRNAs/ /scratch/manuel.tardaguila/GRCH38_ensembl/combined_decoy_reference_miRNAs.fa

  7.2 for snoRNAs

  $ sbatch ~/Scripts/sbatch/12_snoRNAs_STAR_indexing.sh /scratch/manuel.tardaguila/STAR_indexed_genome_for_snoRNAs/ /scratch/manuel.tardaguila/GRCH38_ensembl/combined_decoy_reference_snoRNAs.fa

# 8. Trim in the right way the TAKARA reads from the small RNA library prep

$ bash ~/Scripts/Wraper_scripts/186_trimmomatic_for_smallRNAs_from_TAKARA.sh /scratch/manuel.tardaguila/my_own_RNA_pipeline/ new_trimming /group/soranzo/paola.benaglio/small_rna/250818_A01481_0333_AHGVYHDSXF/fastq_raw/


# 9. Create a hierarchical SAF file for feature counts

  # 9.1 for miRNAs

$  head -n 1 miRNA_converted.tsv > prioritized_miRNA_Decoys_master.saf

$  awk 'NR>1' miRNA_converted.tsv >> prioritized_miRNA_Decoys_master.saf

$  awk 'NR>1' rRNA_converted.tsv >> prioritized_miRNA_Decoys_master.saf

$  awk 'NR>1' tRNA_converted.tsv >> prioritized_miRNA_Decoys_master.saf

$  awk 'NR>1' snord_converted.tsv >> prioritized_miRNA_Decoys_master.saf

  # 9.2 for snoRNAs


$  head -n 1 snord_converted.tsv > prioritized_snoRNA_Decoys_master.saf

$  awk 'NR>1' snord_converted.tsv >> prioritized_snoRNA_Decoys_master.saf

$  awk 'NR>1' rRNA_converted.tsv >> prioritized_snoRNA_Decoys_master.saf

$  awk 'NR>1' tRNA_converted.tsv >> prioritized_snoRNA_Decoys_master.saf

$  awk 'NR>1' miRNA_converted.tsv >> prioritized_snoRNA_Decoys_master.saf


# 10. Alingment with STAR and counts with feature counts

  # 10.1 for miRNAs
  
  $ bash ~/Scripts/Wraper_scripts/187_STAR_alignment_for_snoRNAs.sh /scratch/manuel.tardaguila/ miRNAs_STAR /scratch/manuel.tardaguila/STAR_indexed_genome_for_miRNAs/ /group/soranzo/manuel.tardaguila/small_rna/Dependencies/prioritized_miRNA_Decoys_master.saf miRNAs

  # 10.2 for snoRNAs

  $ bash ~/Scripts/Wraper_scripts/187_STAR_alignment_for_snoRNAs.sh /scratch/manuel.tardaguila/ snoRNAs_STAR /scratch/manuel.tardaguila/STAR_indexed_genome_for_snoRNAs/ /group/soranzo/manuel.tardaguila/small_rna/Dependencies/prioritized_snoRNA_Decoys_master.saf snoRNAs

# 11. DE Jupyte notebooks

  # 11.1 for miRNAs

  # 11.2 for snoRNAs

----> /group/soranzo/manuel.tardaguila/small_rna/snoRNAs/DE_snoRNA_STAR_and_feactureCounts.ipynb








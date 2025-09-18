#!/bin/bash

samples_file=$1
bams_dir=$2
coordinates=$3
declare -a sample_ids
declare -a sample_sexes

# Read samples and sexes into arrays
while read -r id sex; do
    [[ "$id" =~ ^#.* ]] && continue
    [ -z "$id" ] && continue
    sample_ids+=("$id")
    sample_sexes+=("$sex")
done < "$samples_file"

# Process each sample
for i in "${!sample_ids[@]}"; do
    sample="${sample_ids[$i]}"
    sex="${sample_sexes[$i]}"
    echo "Processing sample: $sample (Sex: $sex)"
    
    # Find the actual output directory using glob expansion
    output_dirs=(./whole-genome-OPSIN-contigs/${sample}*)
    if [ ${#output_dirs[@]} -eq 0 ] || [ ! -d "${output_dirs[0]}" ]; then
        echo "ERROR: No output directory found for sample $sample"
        continue
    fi
    output_dir="${output_dirs[0]}"
    echo "Using output directory: $output_dir"
    
    # Find the actual BAM files using glob expansion
    hap1_bams=(${bams_dir}/${sample}.ONT.R10/${sample}*.hap1.p_ctg_only_primary.bam)
    hap2_bams=(${bams_dir}/${sample}.ONT.R10/${sample}*.hap2.p_ctg_only_primary.bam)
    
    if [ ${#hap1_bams[@]} -eq 0 ] || [ ! -f "${hap1_bams[0]}" ]; then
        echo "ERROR: No hap1 BAM file found for sample $sample"
        echo "Looked for: ${bams_dir}/${sample}.ONT.R10/${sample}*.hap1.p_ctg_only_primary.bam"
        continue
    fi
    
    if [ ${#hap2_bams[@]} -eq 0 ] || [ ! -f "${hap2_bams[0]}" ]; then
        echo "ERROR: No hap2 BAM file found for sample $sample"
        echo "Looked for: ${bams_dir}/${sample}.ONT.R10/${sample}*.hap2.p_ctg_only_primary.bam"
        continue
    fi
    
    primary_reads_bam_file_hap1="${hap1_bams[0]}"
    primary_reads_bam_file_hap2="${hap2_bams[0]}"
    
    # Extract sample tag from the actual filename
    sample_tag=$(basename "$primary_reads_bam_file_hap1" .hap1.p_ctg_only_primary.bam)
    
    echo "Found hap1 BAM: $primary_reads_bam_file_hap1"
    echo "Found hap2 BAM: $primary_reads_bam_file_hap2"
    echo "Using sample tag: $sample_tag"
    
    echo "Extracting contigs from the OPSIN region for hap1"
    samtools view "$primary_reads_bam_file_hap1" "$coordinates" | awk '{print ">"$1; print $10}' > "$output_dir/${sample_tag}_${sex}.OPSIN-whole-genome-contigs.asm.bp.hap1.p_ctg.fa"

    echo "Extracting contigs from the OPSIN region for hap2"
    samtools view "$primary_reads_bam_file_hap2" "$coordinates" | awk '{print ">"$1; print $10}' > "$output_dir/${sample_tag}_${sex}.OPSIN-whole-genome-contigs.asm.bp.hap2.p_ctg.fa"

    echo "Completed processing for sample $sample"
    echo "---"
done

    
    





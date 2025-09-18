#! /bin/bash

# Parse command line options
while getopts "s:a:n:o:d:r:" opt; do
    case $opt in
        s) sample_list=$OPTARG;;
        a) assembly_dir=$OPTARG;;
        n) annotation_dir=$OPTARG;;
        o) output_base_dir=$OPTARG;;
        d) descriptor=$OPTARG;;
        r) reference_genome=$OPTARG;;
    esac
done

gff_to_bed_path="scripts/gff_to_bed.py"
# Function to log commands before executing them
execute_and_log() {
    local log_file="$1"
    shift
    echo "$(date '+%Y-%m-%d %H:%M:%S') - EXECUTING: $*" >> "$log_file"
    if "$@" 2>> "$log_file"; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') - SUCCESS: Command completed" >> "$log_file"
        return 0
    else
        local exit_code=$?
        echo "$(date '+%Y-%m-%d %H:%M:%S') - ERROR: Command failed with exit code $exit_code" >> "$log_file"
        return $exit_code
    fi
}

# Function for general logging
log_info() {
    local log_file="$1"
    local message="$2"
    echo "$(date '+%Y-%m-%d %H:%M:%S') - INFO: $message" >> "$log_file"
}

# Function for error logging
log_error() {
    local log_file="$1"
    local message="$2"
    echo "$(date '+%Y-%m-%d %H:%M:%S') - ERROR: $message" >> "$log_file"
}

# Validate required parameters
if [[ -z "$sample_list" || -z "$assembly_dir" || -z "$annotation_dir" || -z "$output_base_dir" || -z "$descriptor" || -z "$reference_genome" ]]; then
    echo "Error: Missing required parameters"
    echo "Usage: $0 -s sample_list -a assembly_dir -n annotation_dir -o output_base_dir -d descriptor -r reference_genome"
    exit 1
fi

# Process each sample
while read -r sample sex; do
    # Skip empty lines or comments
    [[ -z "$sample" || "$sample" =~ ^#.*$ ]] && continue
    
    # Create descriptor parent directory and sample-specific subdirectory
    descriptor_dir="$output_base_dir/$descriptor"
    sample_output_dir="$descriptor_dir/${sample}-${sex}-${descriptor}"
    mkdir -p "$sample_output_dir"
    SAMPLE_LOG="$sample_output_dir/${sample}-${sex}-${descriptor}-full.log"
    
    # Initialize comprehensive log
    {
        echo "========================================"
        echo "SAMPLE PROCESSING LOG"
        echo "========================================"
        echo "Sample: $sample"
        echo "Sex: $sex"
        echo "Start time: $(date)"
        echo "Script: $0"
        echo "Process ID: $$"
        echo "User: $(whoami)"
        echo "Working directory: $(pwd)"
        echo "Parameters:"
        echo "  - Sample list: $sample_list"
        echo "  - Assembly dir: $assembly_dir"
        echo "  - Annotation dir: $annotation_dir"
        echo "  - Output base dir: $output_base_dir"
        echo "  - Descriptor: $descriptor"
        echo "  - Reference genome: $reference_genome"
        echo "========================================"
    } >> "$SAMPLE_LOG"

    log_info "$SAMPLE_LOG" "Starting processing for sample: $sample (sex: $sex)"

    if [ "$sex" == "XX" ]; then
        log_info "$SAMPLE_LOG" "Processing diploid sample (XX)"
        
        # Find haplotype 1 files
        log_info "$SAMPLE_LOG" "Searching for haplotype 1 assembly in: $assembly_dir"
        hap1_assembly=$(find "$assembly_dir" -name "$sample*asm.hap1.p_ctg.fa" | head -1)
        if [[ -n "$hap1_assembly" ]]; then
            log_info "$SAMPLE_LOG" "Found hap1 assembly: $hap1_assembly"
        else
            log_error "$SAMPLE_LOG" "No hap1 assembly found for pattern: $sample*asm.hap1.p_ctg.fa"
            continue
        fi
        
        log_info "$SAMPLE_LOG" "Searching for haplotype 1 GFF in: $annotation_dir"
        hap1_gff=$(find "$annotation_dir" -name "$sample*_hap1.gff" | head -1)
        if [[ -n "$hap1_gff" ]]; then
            log_info "$SAMPLE_LOG" "Found hap1 GFF: $hap1_gff"
        else
            log_error "$SAMPLE_LOG" "No hap1 GFF found for pattern: $sample*_hap1.gff"
            continue
        fi

        # Find haplotype 2 files
        log_info "$SAMPLE_LOG" "Searching for haplotype 2 assembly in: $assembly_dir"
        hap2_assembly=$(find "$assembly_dir" -name "$sample*asm.hap2.p_ctg.fa" | head -1)
        if [[ -n "$hap2_assembly" ]]; then
            log_info "$SAMPLE_LOG" "Found hap2 assembly: $hap2_assembly"
        else
            log_error "$SAMPLE_LOG" "No hap2 assembly found for pattern: $sample*asm.hap2.p_ctg.fa"
            continue
        fi
        
        log_info "$SAMPLE_LOG" "Searching for haplotype 2 GFF in: $annotation_dir"
        hap2_gff=$(find "$annotation_dir" -name "$sample*_hap2.gff" | head -1)
        if [[ -n "$hap2_gff" ]]; then
            log_info "$SAMPLE_LOG" "Found hap2 GFF: $hap2_gff"
        else
            log_error "$SAMPLE_LOG" "No hap2 GFF found for pattern: $sample*_hap2.gff"
            continue
        fi

        # Find FASTQ files
        log_info "$SAMPLE_LOG" "Searching for FASTQ files in: $assembly_dir"
        fastq=$(find "$assembly_dir" -name "$sample*fastq" | head -1)
        if [[ -n "$fastq" ]]; then
            log_info "$SAMPLE_LOG" "Found FASTQ: $fastq"
        else
            log_error "$SAMPLE_LOG" "No FASTQ found for pattern: $sample*fastq"
            continue
        fi

        # Extract filename and setup output paths
        filename=$(basename "$hap1_assembly" | cut -d '.' -f1)
        log_info "$SAMPLE_LOG" "Extracted filename: $filename"
        
        output_dir="$descriptor_dir/${filename}-${sex}-${descriptor}"
        log_info "$SAMPLE_LOG" "Creating output directory: $output_dir"
        execute_and_log "$SAMPLE_LOG" mkdir -p "$output_dir"

        # Combine GFF files
        combined_gff="$output_dir/${filename}-${sex}-${descriptor}-combined-diploid.gff"
        log_info "$SAMPLE_LOG" "Combining GFF files into: $combined_gff"
        execute_and_log "$SAMPLE_LOG" sh -c "cat '$hap1_gff' '$hap2_gff' > '$combined_gff'"
        
        # Combine assembly files
        combined_assembly="$output_dir/${filename}-${sex}-${descriptor}-combined-diploid.fa"
        log_info "$SAMPLE_LOG" "Combining assembly files into: $combined_assembly"
        execute_and_log "$SAMPLE_LOG" sh -c "cat '$hap1_assembly' '$hap2_assembly' > '$combined_assembly'"

        # Convert GFF to BED
        log_info "$SAMPLE_LOG" "Converting GFF to BED format"
        python_output=$(python3 $gff_to_bed_path --gff_file "$combined_gff" --sample_name "$sample" --sex "$sex" --output_dir "$output_dir" 2>&1)
        python_exit_code=$?
        
        # Log Python output regardless of success/failure
        echo "$(date '+%Y-%m-%d %H:%M:%S') - PYTHON OUTPUT: $python_output" >> "$SAMPLE_LOG"
        
        if [ $python_exit_code -ne 0 ]; then
            log_error "$SAMPLE_LOG" "GFF to BED conversion failed with exit code $python_exit_code"
            log_error "$SAMPLE_LOG" "Python error output: $python_output"
            echo "ERROR: GFF to BED conversion failed for $sample"
            echo "Python output: $python_output"
        else
            log_info "$SAMPLE_LOG" "GFF to BED conversion completed successfully"
        fi

        # Alignment process
        output_bam="$output_dir/${filename}-${sex}-${descriptor}-combined-diploid.bam"
        log_info "$SAMPLE_LOG" "Starting alignment process, output BAM: $output_bam"
        
        # Minimap2 alignment
        temp_sam="${sample}-temp.sam"
        log_info "$SAMPLE_LOG" "Running minimap2 alignment"
        if execute_and_log "$SAMPLE_LOG" minimap2 -t 10 -ax map-ont -secondary=no "$combined_assembly" "$fastq" > "$temp_sam"; then
            log_info "$SAMPLE_LOG" "Minimap2 alignment completed successfully"
        else
            log_error "$SAMPLE_LOG" "Minimap2 alignment failed, skipping sample"
            continue
        fi

        # Convert SAM to BAM
        temp_bam="${sample}-temp.bam"
        log_info "$SAMPLE_LOG" "Converting SAM to BAM"
        if execute_and_log "$SAMPLE_LOG" samtools view -bS "$temp_sam" -o "$temp_bam"; then
            log_info "$SAMPLE_LOG" "SAM to BAM conversion completed"
        else
            log_error "$SAMPLE_LOG" "SAM to BAM conversion failed, skipping sample"
            continue
        fi

        # Sort BAM
        log_info "$SAMPLE_LOG" "Sorting BAM file"
        if execute_and_log "$SAMPLE_LOG" samtools sort -o "$output_bam" "$temp_bam"; then
            log_info "$SAMPLE_LOG" "BAM sorting completed"
        else
            log_error "$SAMPLE_LOG" "BAM sorting failed, skipping sample"
            continue
        fi

        # Index BAM
        log_info "$SAMPLE_LOG" "Indexing BAM file"
        if execute_and_log "$SAMPLE_LOG" samtools index "$output_bam"; then
            log_info "$SAMPLE_LOG" "BAM indexing completed"
        else
            log_error "$SAMPLE_LOG" "BAM indexing failed"
        fi

        # Align combined diploid assembly to reference genome
        ref_sam="$output_dir/${filename}-${sex}-${descriptor}-combined-diploid.aln-to-ref.sam"
        ref_bam="$output_dir/${filename}-${sex}-${descriptor}-combined-diploid.aln-to-ref.bam"
        log_info "$SAMPLE_LOG" "Aligning combined diploid assembly to reference genome"
        if execute_and_log "$SAMPLE_LOG" minimap2 -t 10 -x asm10 -a -secondary=no "$reference_genome" "$combined_assembly" > "$ref_sam"; then
            log_info "$SAMPLE_LOG" "Assembly to reference SAM generated"
        else
            log_error "$SAMPLE_LOG" "Assembly to reference alignment failed, skipping sample"
            continue
        fi

        # Convert, sort, and index reference alignment BAM
        log_info "$SAMPLE_LOG" "Converting SAM to BAM for reference alignment"
        if execute_and_log "$SAMPLE_LOG" samtools view -bS "$ref_sam" -o "$ref_bam"; then
            log_info "$SAMPLE_LOG" "Reference SAM to BAM conversion completed"
        else
            log_error "$SAMPLE_LOG" "Reference SAM to BAM conversion failed, skipping sample"
            continue
        fi

        log_info "$SAMPLE_LOG" "Sorting reference alignment BAM"
        if execute_and_log "$SAMPLE_LOG" samtools sort -o "$ref_bam" "$ref_bam"; then
            log_info "$SAMPLE_LOG" "Reference BAM sorting completed"
        else
            log_error "$SAMPLE_LOG" "Reference BAM sorting failed, skipping sample"
            continue
        fi

        log_info "$SAMPLE_LOG" "Indexing reference alignment BAM"
        if execute_and_log "$SAMPLE_LOG" samtools index "$ref_bam"; then
            log_info "$SAMPLE_LOG" "Reference BAM indexing completed"
        else
            log_error "$SAMPLE_LOG" "Reference BAM indexing failed"
        fi

        # Clean up reference SAM
        log_info "$SAMPLE_LOG" "Cleaning up reference SAM file"
        execute_and_log "$SAMPLE_LOG" rm -f "$ref_sam"

        # Clean up temporary files
        log_info "$SAMPLE_LOG" "Cleaning up temporary files"
        execute_and_log "$SAMPLE_LOG" rm -f "$temp_sam" "$temp_bam"

    elif [ "$sex" == "XY" ]; then
        log_info "$SAMPLE_LOG" "Processing haploid sample (XY)"
        
        # Find primary assembly
        log_info "$SAMPLE_LOG" "Searching for primary assembly in: $assembly_dir"
        primary_assembly=$(find "$assembly_dir" -name "$sample*asm.bp.p_ctg.fa" | head -1)
        if [[ -n "$primary_assembly" ]]; then
            log_info "$SAMPLE_LOG" "Found primary assembly: $primary_assembly"
        else
            log_error "$SAMPLE_LOG" "No primary assembly found for pattern: $sample*asm.bp.p_ctg.fa"
            continue
        fi
        
        # Find primary GFF
        log_info "$SAMPLE_LOG" "Searching for primary GFF in: $annotation_dir"
        query_gff=$(find "$annotation_dir" -name "$sample*_primary.gff" | head -1)
        if [[ -n "$query_gff" ]]; then
            log_info "$SAMPLE_LOG" "Found primary GFF: $query_gff"
        else
            log_error "$SAMPLE_LOG" "No primary GFF found for pattern: $sample*_primary.gff"
            continue
        fi

        # Find FASTQ files
        log_info "$SAMPLE_LOG" "Searching for FASTQ files in: $assembly_dir"
        fastq=$(find "$assembly_dir" -name "$sample*fastq" | head -1)
        if [[ -n "$fastq" ]]; then
            log_info "$SAMPLE_LOG" "Found FASTQ: $fastq"
        else
            log_error "$SAMPLE_LOG" "No FASTQ found for pattern: $sample*fastq"
            continue
        fi

        # Setup output paths
        filename=$(basename "$primary_assembly" | cut -d '.' -f1)
        log_info "$SAMPLE_LOG" "Extracted filename: $filename"
        
        output_dir="$descriptor_dir/${filename}-${sex}-${descriptor}"
        log_info "$SAMPLE_LOG" "Creating output directory: $output_dir"
        execute_and_log "$SAMPLE_LOG" mkdir -p "$output_dir"

        # Convert GFF to BED
        log_info "$SAMPLE_LOG" "Converting GFF to BED format"
        python_output=$(python3 $gff_to_bed_path --gff_file "$query_gff" --sample_name "$sample" --sex "$sex" --output_dir "$output_dir" 2>&1)
        python_exit_code=$?
        
        # Log Python output regardless of success/failure
        echo "$(date '+%Y-%m-%d %H:%M:%S') - PYTHON OUTPUT: $python_output" >> "$SAMPLE_LOG"
        
        if [ $python_exit_code -ne 0 ]; then
            log_error "$SAMPLE_LOG" "GFF to BED conversion failed with exit code $python_exit_code"
            log_error "$SAMPLE_LOG" "Python error output: $python_output"
            echo "ERROR: GFF to BED conversion failed for $sample"
            echo "Python output: $python_output"
        else
            log_info "$SAMPLE_LOG" "GFF to BED conversion completed successfully"
        fi

        # Alignment process
        output_bam="$output_dir/${filename}-${sex}-${descriptor}.bam"
        log_info "$SAMPLE_LOG" "Starting alignment process, output BAM: $output_bam"
        
        # Minimap2 alignment
        temp_sam="${sample}-temp.sam"
        log_info "$SAMPLE_LOG" "Running minimap2 alignment"
        if execute_and_log "$SAMPLE_LOG" minimap2 -t 10 -ax map-ont -secondary=no "$primary_assembly" "$fastq" > "$temp_sam"; then
            log_info "$SAMPLE_LOG" "Minimap2 alignment completed successfully"
        else
            log_error "$SAMPLE_LOG" "Minimap2 alignment failed, skipping sample"
            continue
        fi

        # Convert SAM to BAM
        temp_bam="${sample}-temp.bam"
        log_info "$SAMPLE_LOG" "Converting SAM to BAM"
        if execute_and_log "$SAMPLE_LOG" samtools view -bS "$temp_sam" -o "$temp_bam"; then
            log_info "$SAMPLE_LOG" "SAM to BAM conversion completed"
        else
            log_error "$SAMPLE_LOG" "SAM to BAM conversion failed, skipping sample"
            continue
        fi

        # Sort BAM
        log_info "$SAMPLE_LOG" "Sorting BAM file"
        if execute_and_log "$SAMPLE_LOG" samtools sort -o "$output_bam" "$temp_bam"; then
            log_info "$SAMPLE_LOG" "BAM sorting completed"
        else
            log_error "$SAMPLE_LOG" "BAM sorting failed, skipping sample"
            continue
        fi

        # Index BAM
        log_info "$SAMPLE_LOG" "Indexing BAM file"
        if execute_and_log "$SAMPLE_LOG" samtools index "$output_bam"; then
            log_info "$SAMPLE_LOG" "BAM indexing completed"
        else
            log_error "$SAMPLE_LOG" "BAM indexing failed"
        fi

        # Clean up temporary files
        log_info "$SAMPLE_LOG" "Cleaning up temporary files"
        execute_and_log "$SAMPLE_LOG" rm -f "$temp_sam" "$temp_bam"

        # Align primary assembly to reference genome
        xy_ref_sam="$output_dir/${filename}-${sex}-${descriptor}-aln-to-ref.sam"
        xy_ref_bam="$output_dir/${filename}-${sex}-${descriptor}-aln-to-ref.bam"
        log_info "$SAMPLE_LOG" "Aligning primary assembly to reference genome"
        if execute_and_log "$SAMPLE_LOG" minimap2 -t 10 -x asm10 -a -secondary=no "$reference_genome" "$primary_assembly" > "$xy_ref_sam"; then
            log_info "$SAMPLE_LOG" "Primary assembly to reference SAM generated"
        else
            log_error "$SAMPLE_LOG" "Primary assembly to reference alignment failed, skipping sample"
            continue
        fi

        log_info "$SAMPLE_LOG" "Converting SAM to BAM for reference alignment"
        if execute_and_log "$SAMPLE_LOG" samtools view -bS "$xy_ref_sam" -o "$xy_ref_bam"; then
            log_info "$SAMPLE_LOG" "Reference SAM to BAM conversion completed"
        else
            log_error "$SAMPLE_LOG" "Reference SAM to BAM conversion failed, skipping sample"
            continue
        fi

        log_info "$SAMPLE_LOG" "Sorting reference alignment BAM"
        if execute_and_log "$SAMPLE_LOG" samtools sort -o "$xy_ref_bam" "$xy_ref_bam"; then
            log_info "$SAMPLE_LOG" "Reference BAM sorting completed"
        else
            log_error "$SAMPLE_LOG" "Reference BAM sorting failed, skipping sample"
            continue
        fi

        log_info "$SAMPLE_LOG" "Indexing reference alignment BAM"
        if execute_and_log "$SAMPLE_LOG" samtools index "$xy_ref_bam"; then
            log_info "$SAMPLE_LOG" "Reference BAM indexing completed"
        else
            log_error "$SAMPLE_LOG" "Reference BAM indexing failed"
        fi

        # Clean up reference SAM
        log_info "$SAMPLE_LOG" "Cleaning up reference SAM file"
        execute_and_log "$SAMPLE_LOG" rm -f "$xy_ref_sam"

    else
        log_error "$SAMPLE_LOG" "Unknown sex designation: $sex (expected XX or XY)"
        continue
    fi

    # Final log entry
    {
        echo "========================================"
        echo "PROCESSING COMPLETED FOR SAMPLE: $sample"
        echo "End time: $(date)"
        echo "========================================"
        echo ""
    } >> "$SAMPLE_LOG"

done < "$sample_list"

echo "All samples processed. Check individual log files for details."
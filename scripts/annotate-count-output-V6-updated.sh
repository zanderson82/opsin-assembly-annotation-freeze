#!/bin/bash

# OPSIN Gene Copy Number and Array Structure Analysis Pipeline V4 (CORRECTED)
# Processes samples based on metadata table
# Uses TSV output from Python analysis script

# =============================================================================
# Input Parameters
# =============================================================================

METADATA_TABLE=$1       # Table with sample IDs and sex info
ASSEMBLIES_DIR=$2       # Directory containing all assemblies
OPSIN_REFERENCE=$3      # OPSIN reference sequences
OUTPUT_DIR=$4           # Output directory
DESCRIPTION=$5          # Description for this analysis run
BAM_DIR=$6              # Directory containing BAM files (optional)

# Scripts directory (adjust path as needed)
SCRIPTS_DIR="./scripts"
ANALYSIS_SCRIPT="${SCRIPTS_DIR}/python_opsin_processing_V6.py"
STATS_SCRIPT="${SCRIPTS_DIR}/bed-file-contig-alignment-stats.sh"
GFF_TO_BED_SCRIPT="${SCRIPTS_DIR}/gff_to_bed2.py"

CHM13_PATH="/n/dat/chm13/T2T-CHM13v2.fasta"
HG38_PATH="/n/dat/hg38/hg38.no_alt.fa"

# =============================================================================
# Validation
# =============================================================================

if [ -z "$METADATA_TABLE" ] || [ -z "$ASSEMBLIES_DIR" ] || [ -z "$OPSIN_REFERENCE" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$DESCRIPTION" ]; then
    echo "Usage: $0 <metadata_table.tsv> <assemblies_directory> <opsin_reference.fasta> <output_directory> <description> [bam_directory]"
    echo ""
    echo "Metadata table format: SampleID<tab>Sex (XX or XY)"
    echo ""
    echo "bam_directory - location of BAM files for alignment stats"
    exit 1
fi

# Check analysis script exists
if [ ! -f "$ANALYSIS_SCRIPT" ]; then
    echo "ERROR: Analysis script not found: $ANALYSIS_SCRIPT"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$BAM_DIR"

echo "=== Starting OPSIN Gene Analysis Pipeline V4 ==="
echo "Description: $DESCRIPTION"
echo "Output directory: $OUTPUT_DIR"
echo ""

# =============================================================================
# Summary File Setup
# =============================================================================

SUMMARY_FILE="${OUTPUT_DIR}/opsin_array_summary_${DESCRIPTION}.tsv"

# Create summary file with header if it doesn't exist or is empty
if [ ! -s "$SUMMARY_FILE" ]; then
    echo "Creating summary file: $SUMMARY_FILE"
    python3 "$ANALYSIS_SCRIPT" --header > "$SUMMARY_FILE"
fi

# =============================================================================
# Functions
# =============================================================================

run_exonerate() {
    local ASSEMBLY=$1
    local ANNOTATION_OUTPUT_PREFIX=$2
    
    echo "  Running exonerate annotation..."
    
    exonerate --model est2genome --bestn 10 --showalignment no --showvulgar yes \
        --percent 98 \
        --ryo ">%qi|%ti|%qab-%qae|%tab-%tae|%s\n%tas\n" \
        --query "$OPSIN_REFERENCE" --target "$ASSEMBLY" > "${ANNOTATION_OUTPUT_PREFIX}.exonerate" 2>/dev/null
    
    # Extract vulgar format
    grep "^vulgar:" "${ANNOTATION_OUTPUT_PREFIX}.exonerate" > "${ANNOTATION_OUTPUT_PREFIX}.vulgar" 2>/dev/null || true
    
    local VULGAR_LINES=$(wc -l < "${ANNOTATION_OUTPUT_PREFIX}.vulgar" 2>/dev/null || echo "0")
    echo "  Found $VULGAR_LINES vulgar lines"
    
    return 0
}

convert_vulgar_to_gff() {
    local VULGAR_FILE=$1
    local GFF_FILE=$2
    
    echo "  Converting vulgar to GFF format..."
    
    awk '{
        query_name = $2
        target_name = $6
        target_start = $7
        target_end = $8
        target_strand = $9
        score = $10
        
        query_lower = tolower(query_name)
        
        if (query_lower ~ /lcr/) {
            type = "LCR"
        } else if (query_lower ~ /lw/ && query_lower ~ /exon/) {
            type = "OPN1LW_exon5"
        } else if (query_lower ~ /mw/ && query_lower ~ /exon/) {
            type = "OPN1MW_exon5"
        } else {
            next
        }
        
        if (target_start > target_end) {
            temp = target_start
            target_start = target_end
            target_end = temp
        }
        
        print target_name "\t" "exonerate:est2genome" "\t" "similarity" "\t" target_start "\t" target_end "\t" score "\t" target_strand "\t" "." "\t" "query=" query_name ";score=" score
    }' "$VULGAR_FILE" > "$GFF_FILE"
    
    local GFF_LINES=$(wc -l < "$GFF_FILE" 2>/dev/null || echo "0")
    echo "  Created GFF with $GFF_LINES annotations"
    
    return 0
}

align_to_assembly() {
    local ASSEMBLY=$1
    local FASTQ=$2
    local BAM_OUTPUT=$3
    
    echo "  Aligning reads to assembly..."
    echo "    Assembly: $ASSEMBLY"
    echo "    FASTQ: $FASTQ"
    echo "    Output BAM: $BAM_OUTPUT"

    # Check inputs exist
    if [ ! -f "$ASSEMBLY" ]; then
        echo "  ERROR: Assembly file not found: $ASSEMBLY"
        return 1
    fi
    
    if [ ! -f "$FASTQ" ]; then
        echo "  ERROR: FASTQ file not found: $FASTQ"
        return 1
    fi

    # Index assembly if needed
    if [ ! -f "${ASSEMBLY}.fai" ]; then
        echo "  Indexing assembly..."
        samtools faidx "$ASSEMBLY"
    fi

    # Minimap2 alignment, sort, and index in one pipeline
    minimap2 -t 10 -ax map-ont -secondary=no "$ASSEMBLY" "$FASTQ" 2>/dev/null | \
        samtools sort -o "$BAM_OUTPUT" -
    
    samtools index "$BAM_OUTPUT"
    
    echo "  Alignment complete: $BAM_OUTPUT"
    return 0
}

add_alignment_stats() {
    local INPUT_BED=$1
    local SAMPLE_DIR=$2
    local BAM_FILE=$3
    local OUTPUT_GFF=$4
    
    echo "  Adding alignment statistics..."
    echo "    Input BED: $INPUT_BED"
    echo "    BAM file: $BAM_FILE"
    echo "    Output GFF: $OUTPUT_GFF"
    
    # Create output GFF with header
    echo -e "contig\tstart\tend\ttype\tstrand\treads\tMQ0\tratio" > "$OUTPUT_GFF"
    
    # Check input BED exists
    if [ ! -f "$INPUT_BED" ]; then
        echo "  WARNING: Input BED file not found: $INPUT_BED"
        return 1
    fi
    
    # Check BAM file exists
    if [ ! -f "$BAM_FILE" ]; then
        echo "  WARNING: BAM file not found: $BAM_FILE"
        # Copy BED to GFF with placeholder stats
        # BED format from gff_to_bed.py: contig, start, end, type, score, strand
        tail -n +2 "$INPUT_BED" | while IFS=$'\t' read -r contig start end type score strand; do
            echo -e "${contig}\t${start}\t${end}\t${type}\t${strand}\t0\t0\t1.0" >> "$OUTPUT_GFF"
        done
        return 0
    fi
    
    # Process each annotation (skip header)
    # BED format from gff_to_bed.py: contig, start, end, type, score, strand
    tail -n +2 "$INPUT_BED" | while IFS=$'\t' read -r contig start end type score strand; do
        # Skip empty lines
        [ -z "$contig" ] && continue
        
        # Sanitize type for filename (remove special chars)
        local safe_type=$(echo "$type" | tr '/' '_')
        
        # Create temp files for this region
        local REGION_PREFIX="${SAMPLE_DIR}/temp_${contig}-${start}-${end}-${safe_type}"
        
        # BED is 0-based, samtools expects 1-based regions
        local samtools_start=$((start + 1))
        
        # Extract reads for this region
        samtools view -b -F 0x900 "$BAM_FILE" "${contig}:${samtools_start}-${end}" > "${REGION_PREFIX}.bam" 2>/dev/null
        
        if [ -s "${REGION_PREFIX}.bam" ]; then
            samtools index "${REGION_PREFIX}.bam" 2>/dev/null
            samtools stats "${REGION_PREFIX}.bam" > "${REGION_PREFIX}.bam.stats" 2>/dev/null
            
            # Extract stats
            reads=$(grep "^SN" "${REGION_PREFIX}.bam.stats" | grep "reads mapped:" | cut -f 3)
            MQ0=$(grep "^SN" "${REGION_PREFIX}.bam.stats" | grep "reads MQ0:" | cut -f 3)
            
            # Calculate ratio
            ratio="NA"
            if [ -n "$reads" ] && [ "$reads" -ne 0 ] && [ -n "$MQ0" ]; then
                ratio=$(echo "scale=10; $MQ0 / $reads" | bc -l)
            fi
            
            # Clean up temp files
            rm -f "${REGION_PREFIX}.bam" "${REGION_PREFIX}.bam.bai" "${REGION_PREFIX}.bam.stats"
        else
            reads="0"
            MQ0="0"
            ratio="NA"
        fi
        
        # Output GFF line (keep BED 0-based coords for consistency)
        echo -e "${contig}\t${start}\t${end}\t${type}\t${strand}\t${reads}\t${MQ0}\t${ratio}" >> "$OUTPUT_GFF"
        
    done
    
    local GFF_LINES=$(tail -n +2 "$OUTPUT_GFF" | wc -l)
    echo "  Created GFF with stats: $GFF_LINES annotations"
    
    return 0
}

convert_gff_to_bed() {
    local GFF_FILE=$1
    local SAMPLE_ID=$2
    local SEX=$3
    local HAP_NAME=$4
    local SAMPLE_DIR=$5
    
    echo "  Converting GFF to BED format..."
    
    python3 "$GFF_TO_BED_SCRIPT" --gff_file "$GFF_FILE" --sample_name "$SAMPLE_ID" --sex "$SEX" --hap_name "$HAP_NAME" --output_dir "$SAMPLE_DIR"

    return 0
}

find_bam_file() {
    local SAMPLE_ID=$1
    local BAM_SAMPLE_DIR=$2
    local BAM_PATTERN=$3
    
    local BAM=""
    
    if [ -n "$BAM_SAMPLE_DIR" ] && [ -d "$BAM_SAMPLE_DIR" ]; then
        BAM=$(find "$BAM_SAMPLE_DIR" -name "*${BAM_PATTERN}*.bam" -type f 2>/dev/null | head -1)
    fi
    
    echo "$BAM"
}

find_fastq_file() {
    local SAMPLE_ID=$1
    local SEARCH_DIR=$2
    
    local FASTQ=""
    
    if [ -n "$SEARCH_DIR" ]; then
        # Try .fastq.gz first
        FASTQ=$(find "$SEARCH_DIR" -name "*${SAMPLE_ID}*.fastq.gz" -type f 2>/dev/null | head -1)
        
        # Try .fastq if .gz not found
        [ -z "$FASTQ" ] && FASTQ=$(find "$SEARCH_DIR" -name "*${SAMPLE_ID}*.fastq" -type f 2>/dev/null | head -1)
        
        # Try .fq.gz
        [ -z "$FASTQ" ] && FASTQ=$(find "$SEARCH_DIR" -name "*${SAMPLE_ID}*.fq.gz" -type f 2>/dev/null | head -1)
    fi
    
    echo "$FASTQ"
}

find_assembly() {
    local ASSEMBLY_SAMPLE_DIR=$1
    local ASSEMBLY_PATTERN=$2
    
    find "$ASSEMBLY_SAMPLE_DIR" -name "$ASSEMBLY_PATTERN" 2>/dev/null | head -1
}

analyze_haplotype() {
    local GFF_FILE=$1
    local SAMPLE_ID=$2
    local SEX=$3
    local HAP_NAME=$4
    
    echo "  Analyzing array structure..."
    
    python3 "$ANALYSIS_SCRIPT" "$GFF_FILE" "$SAMPLE_ID" "$SEX" "$HAP_NAME" >> "$SUMMARY_FILE"
    
    local EXIT_CODE=$?
    if [ $EXIT_CODE -ne 0 ]; then
        echo "  ERROR: Analysis failed with exit code $EXIT_CODE"
        return 1
    fi
    
    return 0
}

process_haplotype() {
    local ASSEMBLY=$1
    local SAMPLE_ID=$2
    local HAP_NAME=$3
    local SEX=$4
    local BAM_FILE=$5  # Pre-computed BAM file path
    
    local ANNOTATION_SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE_ID}_${SEX}_annotation_outputs"
    local ANNOTATION_OUTPUT_PREFIX="${ANNOTATION_SAMPLE_DIR}/${SAMPLE_ID}_${HAP_NAME}"
    
    mkdir -p "$ANNOTATION_SAMPLE_DIR"
    
    echo "Processing $SAMPLE_ID ($SEX, $HAP_NAME)..."
    echo "  Assembly: $ASSEMBLY"
    echo "  BAM file: $BAM_FILE"
    
    # Check if already processed (has GFF with stats)
    if [ -f "${ANNOTATION_OUTPUT_PREFIX}_stats.gff" ]; then
        echo "  GFF with stats already exists, running analysis only..."
        analyze_haplotype "${ANNOTATION_OUTPUT_PREFIX}_stats.gff" "$SAMPLE_ID" "$SEX" "$HAP_NAME"
        return $?
    fi
    
    # Check assembly exists
    if [ ! -f "$ASSEMBLY" ]; then
        echo "  ERROR: Assembly file not found: $ASSEMBLY"
        python3 "$ANALYSIS_SCRIPT" --na-line "$SAMPLE_ID" "$SEX" "$HAP_NAME" >> "$SUMMARY_FILE"
        return 1
    fi
    
    # Step 1: Run exonerate
    run_exonerate "$ASSEMBLY" "$ANNOTATION_OUTPUT_PREFIX"
    
    # Step 2: Convert vulgar to GFF
    convert_vulgar_to_gff "${ANNOTATION_OUTPUT_PREFIX}.vulgar" "${ANNOTATION_OUTPUT_PREFIX}.gff"
    
    # Check if we got any annotations
    if [ ! -s "${ANNOTATION_OUTPUT_PREFIX}.gff" ]; then
        echo "  WARNING: No OPSIN annotations found"
        python3 "$ANALYSIS_SCRIPT" --na-line "$SAMPLE_ID" "$SEX" "$HAP_NAME" >> "$SUMMARY_FILE"
        return 1
    fi

    # Step 3: Convert GFF to BED
    convert_gff_to_bed "${ANNOTATION_OUTPUT_PREFIX}.gff" "$SAMPLE_ID" "$SEX" "$HAP_NAME" "$ANNOTATION_SAMPLE_DIR"

    # Step 4: Add alignment stats
    local BED_FILE="${ANNOTATION_SAMPLE_DIR}/${SAMPLE_ID}_${HAP_NAME}.bed"
    
    if [ -n "$BAM_FILE" ] && [ -f "$BAM_FILE" ]; then
        echo "  Using BAM file for stats: $BAM_FILE"
        add_alignment_stats "$BED_FILE" "$ANNOTATION_SAMPLE_DIR" "$BAM_FILE" "${ANNOTATION_OUTPUT_PREFIX}_stats.gff"
    else
        echo "  WARNING: No BAM file available, using placeholder stats"
        # Create GFF with placeholder stats
        echo -e "contig\tstart\tend\ttype\tstrand\treads\tMQ0\tratio" > "${ANNOTATION_OUTPUT_PREFIX}_stats.gff"
        if [ -f "$BED_FILE" ]; then
            tail -n +2 "$BED_FILE" | while IFS=$'\t' read -r contig start end type strand rest; do
                echo -e "${contig}\t${start}\t${end}\t${type}\t${strand}\t0\t0\t1.0" >> "${ANNOTATION_OUTPUT_PREFIX}_stats.gff"
            done
        fi
    fi
    
    # Step 5: Analyze
    analyze_haplotype "${ANNOTATION_OUTPUT_PREFIX}_stats.gff" "$SAMPLE_ID" "$SEX" "$HAP_NAME"
    
    return $?
}

# =============================================================================
# Main Processing Loop
# =============================================================================

echo ""
echo "Processing samples from: $METADATA_TABLE"
echo ""

while IFS=$'\t' read -r SAMPLE_ID SEX EXTRA_COLS; do
    # Skip header, comments, empty lines
    [[ "$SAMPLE_ID" == "SampleID" || "$SAMPLE_ID" == \#* || -z "$SAMPLE_ID" ]] && continue
    
    echo "========================================"
    echo "Sample: $SAMPLE_ID (Sex: $SEX)"
    echo "========================================"
    
    # Find sample assembly directory
    ASSEMBLY_SAMPLE_DIR=$(find "$ASSEMBLIES_DIR" -maxdepth 2 -type d -name "*${SAMPLE_ID}*" 2>/dev/null | head -1)
    
    if [ -z "$ASSEMBLY_SAMPLE_DIR" ]; then
        echo "ERROR: Could not find directory for sample $SAMPLE_ID"
        if [[ "$SEX" == "XY" ]]; then
            python3 "$ANALYSIS_SCRIPT" --na-line "$SAMPLE_ID" "$SEX" "primary" >> "$SUMMARY_FILE"
        elif [[ "$SEX" == "XX" ]]; then
            python3 "$ANALYSIS_SCRIPT" --na-line "$SAMPLE_ID" "$SEX" "hap1" >> "$SUMMARY_FILE"
            python3 "$ANALYSIS_SCRIPT" --na-line "$SAMPLE_ID" "$SEX" "hap2" >> "$SUMMARY_FILE"
        else
            python3 "$ANALYSIS_SCRIPT" --na-line "$SAMPLE_ID" "$SEX" "NA" >> "$SUMMARY_FILE"
        fi
        continue
    fi
    
    echo "Found assembly directory: $ASSEMBLY_SAMPLE_DIR"
    
    # Define BAM output directory for this sample
    BAM_OUTPUT_DIR="${BAM_DIR}/${SAMPLE_ID}_${SEX}-${DESCRIPTION}"
    mkdir -p "$BAM_OUTPUT_DIR"
    
    # Find FASTQ file for this sample
    FASTQ_FILE=$(find_fastq_file "$SAMPLE_ID" "$ASSEMBLY_SAMPLE_DIR")
    [ -z "$FASTQ_FILE" ] && FASTQ_FILE=$(find_fastq_file "$SAMPLE_ID" "$ASSEMBLIES_DIR")
    
    echo "FASTQ file: ${FASTQ_FILE:-NOT FOUND}"
    
    # Process based on sex
    if [[ "$SEX" == "XY" ]]; then
        # =====================================================================
        # XY: Primary assembly only
        # =====================================================================
        echo "Processing XY sample (primary assembly only)..."
        
        # Find primary assembly
        PRIMARY_ASSEMBLY=$(find_assembly "$ASSEMBLY_SAMPLE_DIR" "*.asm.bp.p_ctg.fa")
        [ -z "$PRIMARY_ASSEMBLY" ] && PRIMARY_ASSEMBLY=$(find_assembly "$ASSEMBLY_SAMPLE_DIR" "*.p_ctg.fa")
        [ -z "$PRIMARY_ASSEMBLY" ] && PRIMARY_ASSEMBLY=$(find_assembly "$ASSEMBLY_SAMPLE_DIR" "*primary*.fa")
        
        if [ -z "$PRIMARY_ASSEMBLY" ] || [ ! -f "$PRIMARY_ASSEMBLY" ]; then
            echo "ERROR: Primary assembly not found for $SAMPLE_ID"
            python3 "$ANALYSIS_SCRIPT" --na-line "$SAMPLE_ID" "$SEX" "primary" >> "$SUMMARY_FILE"
            continue
        fi
        
        echo "Primary assembly: $PRIMARY_ASSEMBLY"
        
        # Define BAM output path
        PRIMARY_BAM="${BAM_OUTPUT_DIR}/${SAMPLE_ID}_reads-to-assembly.bam"
        
        # Step 1: Align reads to primary assembly (if not already done)
        if [ ! -f "$PRIMARY_BAM" ]; then
            if [ -n "$FASTQ_FILE" ] && [ -f "$FASTQ_FILE" ]; then
                align_to_assembly "$PRIMARY_ASSEMBLY" "$FASTQ_FILE" "$PRIMARY_BAM"
            else
                echo "WARNING: No FASTQ file found for $SAMPLE_ID, skipping alignment"
            fi
        else
            echo "BAM file already exists: $PRIMARY_BAM"
        fi
        
        # Step 2: Process haplotype with BAM
        process_haplotype "$PRIMARY_ASSEMBLY" "$SAMPLE_ID" "primary" "$SEX" "$PRIMARY_BAM"
        
    elif [[ "$SEX" == "XX" ]]; then
        # =====================================================================
        # XX: Both haplotypes with combined diploid alignment
        # =====================================================================
        echo "Processing XX sample (both haplotypes)..."
        
        # Find haplotype assemblies
        HAP1_ASSEMBLY=$(find_assembly "$ASSEMBLY_SAMPLE_DIR" "*hap1*.fa")
        HAP2_ASSEMBLY=$(find_assembly "$ASSEMBLY_SAMPLE_DIR" "*hap2*.fa")
        
        echo "Hap1 assembly: ${HAP1_ASSEMBLY:-NOT FOUND}"
        echo "Hap2 assembly: ${HAP2_ASSEMBLY:-NOT FOUND}"
        
        # Check both haplotypes exist
        if [ -z "$HAP1_ASSEMBLY" ] || [ ! -f "$HAP1_ASSEMBLY" ]; then
            echo "ERROR: Haplotype 1 assembly not found for $SAMPLE_ID"
            python3 "$ANALYSIS_SCRIPT" --na-line "$SAMPLE_ID" "$SEX" "hap1" >> "$SUMMARY_FILE"
            HAP1_ASSEMBLY=""
        fi
        
        if [ -z "$HAP2_ASSEMBLY" ] || [ ! -f "$HAP2_ASSEMBLY" ]; then
            echo "ERROR: Haplotype 2 assembly not found for $SAMPLE_ID"
            python3 "$ANALYSIS_SCRIPT" --na-line "$SAMPLE_ID" "$SEX" "hap2" >> "$SUMMARY_FILE"
            HAP2_ASSEMBLY=""
        fi
        
        # Skip if neither haplotype found
        if [ -z "$HAP1_ASSEMBLY" ] && [ -z "$HAP2_ASSEMBLY" ]; then
            echo "ERROR: No haplotype assemblies found, skipping sample"
            continue
        fi
        
        # Step 1: Create combined diploid assembly (if both haplotypes exist)
        COMBINED_DIPLOID_ASSEMBLY="${BAM_OUTPUT_DIR}/${SAMPLE_ID}_${SEX}_combined-diploid.fa"
        
        if [ -n "$HAP1_ASSEMBLY" ] && [ -n "$HAP2_ASSEMBLY" ]; then
            if [ ! -f "$COMBINED_DIPLOID_ASSEMBLY" ]; then
                echo "Creating combined diploid assembly..."
                cat "$HAP1_ASSEMBLY" "$HAP2_ASSEMBLY" > "$COMBINED_DIPLOID_ASSEMBLY"
                samtools faidx "$COMBINED_DIPLOID_ASSEMBLY"
            else
                echo "Combined diploid assembly already exists: $COMBINED_DIPLOID_ASSEMBLY"
            fi
        fi
        
        # Step 2: Align reads to combined diploid assembly BEFORE processing haplotypes
        COMBINED_BAM="${BAM_OUTPUT_DIR}/${SAMPLE_ID}_reads-to-combined-diploid.bam"
        
        if [ ! -f "$COMBINED_BAM" ]; then
            if [ -n "$FASTQ_FILE" ] && [ -f "$FASTQ_FILE" ] && [ -f "$COMBINED_DIPLOID_ASSEMBLY" ]; then
                echo "Aligning reads to combined diploid assembly..."
                align_to_assembly "$COMBINED_DIPLOID_ASSEMBLY" "$FASTQ_FILE" "$COMBINED_BAM"
            else
                echo "WARNING: Cannot create combined BAM - missing FASTQ or combined assembly"
            fi
        else
            echo "Combined diploid BAM already exists: $COMBINED_BAM"
        fi
        
        # Step 3: Process haplotype 1 (using combined BAM for stats)
        if [ -n "$HAP1_ASSEMBLY" ] && [ -f "$HAP1_ASSEMBLY" ]; then
            process_haplotype "$HAP1_ASSEMBLY" "$SAMPLE_ID" "hap1" "$SEX" "$COMBINED_BAM"
        fi
        
        # Step 4: Process haplotype 2 (using combined BAM for stats)
        if [ -n "$HAP2_ASSEMBLY" ] && [ -f "$HAP2_ASSEMBLY" ]; then
            process_haplotype "$HAP2_ASSEMBLY" "$SAMPLE_ID" "hap2" "$SEX" "$COMBINED_BAM"
        fi
        
    else
        echo "WARNING: Unknown sex '$SEX' for sample $SAMPLE_ID (must be XX or XY)"
        python3 "$ANALYSIS_SCRIPT" --na-line "$SAMPLE_ID" "$SEX" "NA" >> "$SUMMARY_FILE"
    fi
    
    echo ""
    
done < "$METADATA_TABLE"

# =============================================================================
# Summary
# =============================================================================

echo "========================================"
echo "Pipeline completed!"
echo "========================================"
echo ""
echo "Summary file: $SUMMARY_FILE"
echo ""
echo "Summary preview (first 10 lines):"
head -10 "$SUMMARY_FILE" | column -t -s $'\t'
echo ""

TOTAL_ENTRIES=$(tail -n +2 "$SUMMARY_FILE" | wc -l)
SUCCESSFUL=$(tail -n +2 "$SUMMARY_FILE" | grep -v $'\tNA\t' | wc -l)
echo "Total entries: $TOTAL_ENTRIES"
echo "Successful analyses: $SUCCESSFUL"
echo "Failed/NA entries: $((TOTAL_ENTRIES - SUCCESSFUL))"
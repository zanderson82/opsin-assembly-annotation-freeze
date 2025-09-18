#!/bin/bash

# OPSIN Gene Copy Number and Array Structure Analysis Pipeline V3
# Processes samples based on metadata table, focuses only on exon5 annotations
# Updated to handle multiple LCRs and group by contig properly

# Input parameters
METADATA_TABLE=$1       # Table with sample IDs and sex info
ASSEMBLIES_DIR=$2       # Directory containing all assemblies
OPSIN_REFERENCE=$3      # OPSIN reference sequences
OUTPUT_DIR=$4           # Output directory

# Check if required parameters are provided
if [ -z "$METADATA_TABLE" ] || [ -z "$ASSEMBLIES_DIR" ] || [ -z "$OPSIN_REFERENCE" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <metadata_table.tsv> <assemblies_directory> <opsin_reference.fasta> <output_directory>"
    echo "Metadata table format: SampleID<tab>Sex (XX or XY)"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

echo "=== Starting OPSIN gene analysis pipeline V3 (Updated for contig-aware analysis) ==="

# Function to process each haplotype
process_haplotype() {
    local ASSEMBLY=$1
    local SAMPLE_ID=$2
    local HAP_NAME=$3
    local SEX=$4
    local SAMPLE_DIR="$OUTPUT_DIR/${SAMPLE_ID}_${SEX}_annotation_outputs"
    local OUTPUT_PREFIX="$SAMPLE_DIR/${SAMPLE_ID}_${HAP_NAME}"
    
    echo "Processing $SAMPLE_ID ($HAP_NAME)..."
    
    # Check if annotation files already exist
    if [ -f "${OUTPUT_PREFIX}_analysis.txt" ] && [ -f "${OUTPUT_PREFIX}.gff" ]; then
        echo "  Annotation files already exist for $SAMPLE_ID ($HAP_NAME), skipping..."
        return 0
    fi
    
    # Check if assembly file exists
    if [ ! -f "$ASSEMBLY" ]; then
        echo "ERROR: Assembly file not found: $ASSEMBLY"
        return 1
    fi
    
    # Run Exonerate to annotate OPSIN exon5 and LCR
    echo "  Annotating OPSIN exon5 and LCR..."
    exonerate --model est2genome --bestn 10 --showalignment no --showvulgar yes \
        --percent 98 \
        --ryo ">%qi|%ti|%qab-%qae|%tab-%tae|%s\n%tas\n" \
        --query $OPSIN_REFERENCE --target $ASSEMBLY > ${OUTPUT_PREFIX}.exonerate
    
    # Extract vulgar format for processing and convert to GFF-like format for our analysis
    echo "  Extracting vulgar format and creating contig-aware GFF file..."
    grep "^vulgar:" ${OUTPUT_PREFIX}.exonerate > ${OUTPUT_PREFIX}.vulgar || echo "Warning: No vulgar lines found in exonerate output"
    
    # Debug: show the number of lines in the vulgar file
    VULGAR_LINES=$(wc -l < ${OUTPUT_PREFIX}.vulgar || echo "0")
    echo "  Found $VULGAR_LINES vulgar lines for conversion to GFF"
    
    # Show the first few lines for debugging if any exist
    if [ "$VULGAR_LINES" -gt 0 ]; then
        echo "  First 3 vulgar lines:"
        head -n 3 ${OUTPUT_PREFIX}.vulgar
    fi
    
    # Convert vulgar to GFF-like format with proper contig information
    # Vulgar format: vulgar: query_name query_start query_end query_strand target_name target_start target_end target_strand score
    # We want: contig_name, source, feature_type, start, end, score, strand, frame, attributes
    awk '{
        # Extract information from vulgar line
        query_name = $2
        target_name = $6  # This is the contig name
        target_start = $7
        target_end = $8
        target_strand = $9
        score = $10
        
        # Print in GFF format: contig_name, source, feature_type, start, end, score, strand, frame, attributes
        print target_name "\t" "exonerate:est2genome" "\t" "similarity" "\t" target_start "\t" target_end "\t" score "\t" target_strand "\t" "." "\t" "query=" query_name ";score=" score
    }' ${OUTPUT_PREFIX}.vulgar > ${OUTPUT_PREFIX}.gff
    
    # Debug: show the number of lines in the GFF file
    GFF_LINES=$(wc -l < ${OUTPUT_PREFIX}.gff || echo "0")
    echo "  Created $GFF_LINES GFF lines for analysis"
    
    # Show the first few GFF lines for debugging if any exist
    if [ "$GFF_LINES" -gt 0 ]; then
        echo "  First 3 GFF lines:"
        head -n 3 ${OUTPUT_PREFIX}.gff
    fi
    
    # Analyze gene copies and array structure focusing only on exon5
    echo "  Analyzing gene copies and array structure by contig..."
    python3 $OUTPUT_DIR/analyze_opsin_exon5_updated_V3.py ${OUTPUT_PREFIX}.gff $SAMPLE_ID $HAP_NAME > ${OUTPUT_PREFIX}_analysis.txt
    
    # Return success
    return 0
}

# Create the updated Python analysis script that handles contigs properly
cat > $OUTPUT_DIR/analyze_opsin_exon5_updated_V3.py << 'EOF'
#!/usr/bin/env python3
import re
import sys
from collections import OrderedDict, defaultdict

def analyze_opsin_exon5(gff_file, sample_id, hap_name):
    # Initialize counters for exon5
    lw_exon5_count = 0
    mw_exon5_count = 0
    
    # Store annotations grouped by contig
    contigs = defaultdict(list)
    
    # Debug all contents of GFF file
    print(f"Analyzing GFF file: {gff_file}")
    
    # Parse the GFF file to find the annotations
    try:
        with open(gff_file, 'r') as f:
            content = f.read()
            print(f"File loaded, total characters: {len(content)}")
            
            # If file is empty
            if not content.strip():
                print("WARNING: GFF file is empty!")
                return 0, 0, "Unknown (Empty GFF)", 0
            
            lines = content.splitlines()
            print(f"Total lines in file: {len(lines)}")
            
            # reading the gff file line by line
            for line_count, line in enumerate(lines, 1):
                if line.startswith('#') or line.strip() == '':
                    continue
                
                # Print the raw line for debugging
                print(f"Line {line_count}: {line.strip()}")
                
                # split the line into fields
                fields = line.strip().split('\t')
                # check if the line has fewer than 9 fields
                if len(fields) < 9:
                    print(f"  Warning: Line {line_count} has fewer than 9 fields ({len(fields)} fields)")
                    continue
                
                # Extract information from GFF format
                contig_name = fields[0]  # Contig/target name
                start_pos = int(fields[3]) if fields[3].isdigit() else 0
                end_pos = int(fields[4]) if fields[4].isdigit() else 0
                strand = fields[6]
                attributes = fields[8]
                
                # Extract query name from attributes
                query_name = ""
                if "query=" in attributes:
                    query_match = re.search(r'query=([^;]+)', attributes)
                    if query_match:
                        query_name = query_match.group(1)
                
                # Debug info
                print(f"  Contig: {contig_name}, Query: {query_name}, Position: {start_pos}-{end_pos}, Strand: {strand}")
                
                # Only process OPSIN-related annotations
                if query_name in ['LCR', 'OPN1LW_exon5', 'OPN1MW_exon5']:
                    # Create annotation object
                    annotation = {
                        'type': query_name,
                        'position': start_pos,
                        'end_position': end_pos,
                        'strand': strand,
                        'contig': contig_name
                    }
                    
                    # Add to contig-specific list
                    contigs[contig_name].append(annotation)
                    
                    # Count exon5 annotations for totals
                    if query_name == 'OPN1LW_exon5':
                        lw_exon5_count += 1
                        print(f"    Found LW exon5 on contig {contig_name} at position {start_pos}")
                    elif query_name == 'OPN1MW_exon5':
                        mw_exon5_count += 1
                        print(f"    Found MW exon5 on contig {contig_name} at position {start_pos}")
                    elif query_name == 'LCR':
                        print(f"    Found LCR on contig {contig_name} at position {start_pos}")
                
    except Exception as e:
        print(f"Error reading GFF file: {e}")
        return 0, 0, "Error", 0
    
    # Calculate total LCR count BEFORE printing results
    total_lcr_count = 0
    for contig_name, annotations in contigs.items():    # Iterate through each contig and its annotations
        lcr_annotations = [a for a in annotations if a['type'] == 'LCR']    # Extract LCR annotations from the contig
        total_lcr_count += len(lcr_annotations)    # Add the number of LCR annotations to the total count
    
    # Print results
    print(f"\n=== Sample {sample_id} ({hap_name}) OPSIN Gene Analysis ===")
    
    print(f"  OPN1LW: {lw_exon5_count} copies (based on exon 5 count)")
    print(f"  OPN1MW: {mw_exon5_count} copies (based on exon 5 count)")
    print(f"  Total OPSIN genes: {lw_exon5_count + mw_exon5_count} copies")
    print(f"  Total LCR count: {total_lcr_count} LCRs")
    
    # Array structure analysis by contig
    print("\n=== Array Structure Analysis by Contig ===")
    
    # Check if we have any annotations
    if not contigs:
        print("WARNING: No OPSIN annotations found in the GFF file!")
        return 0, 0, "Unknown (No annotations)", 0
    
    print(f"Found OPSIN annotations on {len(contigs)} contig(s):")
    for contig_name, annotations in contigs.items():
        lcr_count = len([a for a in annotations if a['type'] == 'LCR'])
        exon5_count = len([a for a in annotations if a['type'] in ['OPN1LW_exon5', 'OPN1MW_exon5']])
        print(f"  Contig {contig_name}: {len(annotations)} annotations ({lcr_count} LCRs, {exon5_count} exon5s)")
    
    # Separate contigs with and without LCRs
    lcr_contigs = []
    non_lcr_contigs = []
    
    for contig_name, annotations in contigs.items():
        lcr_annotations = [a for a in annotations if a['type'] == 'LCR']
        exon5_annotations = [a for a in annotations if a['type'] in ['OPN1LW_exon5', 'OPN1MW_exon5']]
        # Split contigs into LCR and non-LCR contigs and store in separate lists with their exon5 counts and annotations
        if lcr_annotations:
            lcr_contigs.append({
                'name': contig_name,
                'annotations': annotations,
                'lcr_count': len(lcr_annotations),
                'exon5_count': len(exon5_annotations),
                'lcr_annotations': lcr_annotations,
                'exon5_annotations': exon5_annotations
            })
        elif exon5_annotations:  # Only add if it has exon5 annotations
            non_lcr_contigs.append({
                'name': contig_name,
                'annotations': annotations,
                'exon5_count': len(exon5_annotations),
                'exon5_annotations': exon5_annotations
            })
    
    print(f"\nContigs with LCRs: {len(lcr_contigs)}")
    print(f"Contigs with exon5s but no LCRs: {len(non_lcr_contigs)}")
    
    if not lcr_contigs:
        print("WARNING: No contigs with LCR annotations found!")
        return lw_exon5_count, mw_exon5_count, "Unknown (No LCR found)", total_lcr_count
    
    # Sort LCR contigs by number of exon5s (descending) to prioritize the one with most genes
    lcr_contigs.sort(key=lambda x: x['exon5_count'], reverse=True)
    primary_contig = lcr_contigs[0]
    
    print(f"\nPrimary contig (most genes): {primary_contig['name']} with {primary_contig['exon5_count']} exon5s and {primary_contig['lcr_count']} LCRs")
    
    # Analyze primary contig structure
    print(f"\n--- Analyzing Primary Contig {primary_contig['name']} ---")
    primary_structure = analyze_contig_structure(primary_contig['name'], primary_contig['lcr_annotations'], primary_contig['exon5_annotations'], primary_contig['annotations'])
    
    # Collect additional genes from non-LCR contigs
    additional_genes = []
    if non_lcr_contigs:
        print(f"\n--- Processing Additional Contigs without LCRs ---")
        for contig in non_lcr_contigs:
            print(f"Contig {contig['name']}: {contig['exon5_count']} exon5 genes")
            
            # Sort genes by position for consistent ordering
            sorted_genes = sorted(contig['exon5_annotations'], key=lambda x: x['position'])
            additional_genes.extend(sorted_genes)
            
            for i, gene in enumerate(sorted_genes):
                print(f"  {i+1}. {gene['type']} at position {gene['position']}")
    
    # Build final combined structure
    if primary_structure and additional_genes:
        # Extract the gene portion from primary structure (remove orientation info)
        base_structure = primary_structure
        orientation_info = ""
        if "(reverse compliment)" in primary_structure:
            base_structure = primary_structure.replace(" (reverse compliment)", "")
            orientation_info = " (reverse compliment)"
        
        # Add additional genes
        additional_structure = ""
        for gene in additional_genes:
            if gene['type'] == 'OPN1LW_exon5':
                additional_structure += "-L"
            elif gene['type'] == 'OPN1MW_exon5':
                additional_structure += "-M"
        
        final_structure = base_structure + additional_structure + orientation_info
        
        print(f"\nCombined Array Structure:")
        print(f"  Primary contig: {primary_structure}")
        if additional_genes:
            print(f"  Additional genes: {additional_structure}")
        print(f"  Final structure: {final_structure}")
        
        return lw_exon5_count, mw_exon5_count, final_structure, total_lcr_count
    elif primary_structure:
        print(f"\nFinal structure: {primary_structure}")
        return lw_exon5_count, mw_exon5_count, primary_structure, total_lcr_count
    else:
        print("No valid array structure could be determined")
        return lw_exon5_count, mw_exon5_count, "Unknown (No valid arrays found)", total_lcr_count

def analyze_contig_structure(contig_name, lcr_annotations, contig_exon5, all_annotations):
    """Analyze the array structure for a single contig"""
    
    if not lcr_annotations or not contig_exon5:
        return None
    
    # Determine the dominant strand orientation for this contig
    minus_strand_count = sum(1 for a in all_annotations if a['strand'] == '-')
    plus_strand_count = sum(1 for a in all_annotations if a['strand'] == '+')
    is_reverse_orientation = minus_strand_count > plus_strand_count
    
    print(f"Contig {contig_name} strand orientation: {'+' if not is_reverse_orientation else '-'} strand dominant " +
          f"({plus_strand_count} + strands, {minus_strand_count} - strands)")
    
    # Sort annotations on this contig by position
    if is_reverse_orientation:
        # For reverse orientation, sort in descending order
        all_annotations.sort(key=lambda x: x['position'], reverse=True)
        contig_exon5.sort(key=lambda x: x['position'], reverse=True)
        lcr_annotations.sort(key=lambda x: x['position'], reverse=True)
        print(f"Contig {contig_name} is in REVERSE orientation - sorting positions in descending order")
    else:
        # For forward orientation, sort in ascending order
        all_annotations.sort(key=lambda x: x['position'])
        contig_exon5.sort(key=lambda x: x['position'])
        lcr_annotations.sort(key=lambda x: x['position'])
        print(f"Contig {contig_name} is in FORWARD orientation - sorting positions in ascending order")
    
    # Debug: Print sorted annotations for this contig
    relevant_annotations = [a for a in all_annotations if a['type'] in ['LCR', 'OPN1LW_exon5', 'OPN1MW_exon5']]
    print(f"Sorted relevant annotations on contig {contig_name}:")
    for i, a in enumerate(relevant_annotations):
        print(f"  {i+1}. {a['type']} at position {a['position']} on strand {a['strand']}")
    
    # Handle multiple LCRs on this contig by grouping genes with their closest LCR
    if len(lcr_annotations) > 1:
        print(f"Multiple LCR annotations detected on contig {contig_name} ({len(lcr_annotations)}). Analyzing separate arrays:")
        
        # Group each exon5 with its closest LCR on this contig
        arrays = []
        for lcr in lcr_annotations:
            lcr_position = lcr['position']
            # Find exon5s closest to this LCR
            closest_exon5s = []
            for exon5 in contig_exon5:
                # Find which LCR this exon5 is closest to
                distances = [(abs(exon5['position'] - l['position']), l) for l in lcr_annotations]
                closest_lcr = min(distances, key=lambda x: x[0])[1]
                if closest_lcr['position'] == lcr_position:
                    closest_exon5s.append(exon5)
            
            if closest_exon5s:
                arrays.append({
                    'lcr': lcr,
                    'exon5s': closest_exon5s
                })
        
        # Print array information
        for i, array in enumerate(arrays):
            print(f"  Array {i+1} on contig {contig_name}: LCR at {array['lcr']['position']} with {len(array['exon5s'])} genes")
            for j, exon5 in enumerate(array['exon5s']):
                print(f"    Gene {j+1}: {exon5['type']} at position {exon5['position']}")
        
        # Use the array with the most genes
        if arrays:
            primary_array = max(arrays, key=lambda x: len(x['exon5s']))
            primary_exon5s = primary_array['exon5s']
            
            print(f"Using primary array on contig {contig_name} (most genes) for structure determination:")
            print(f"  LCR at position {primary_array['lcr']['position']}")
            print(f"  {len(primary_exon5s)} associated genes")
            
            return build_array_structure(primary_exon5s, is_reverse_orientation)
    else:
        # Single LCR on this contig - use original logic
        first_lcr = lcr_annotations[0]
        lcr_position = first_lcr['position']
        
        print(f"Single LCR found on contig {contig_name} at position {lcr_position}")
        
        # Find the closest exon5 to the LCR based on absolute distance
        closest_exon5 = min(contig_exon5, key=lambda x: abs(x['position'] - lcr_position))
        
        print(f"Closest exon5 to LCR on contig {contig_name}: {closest_exon5['type']} at position {closest_exon5['position']}")
        
        # Reorder contig_exon5 to put the closest one first
        reordered_exon5 = [closest_exon5] + [e for e in contig_exon5 if e != closest_exon5]
        
        # Print the reordered annotations
        print(f"Reordered exon5 annotations on contig {contig_name}:")
        for i, ann in enumerate(reordered_exon5):
            print(f"  {i+1}. {ann['type']} at position {ann['position']} on strand {ann['strand']}")
        
        return build_array_structure(reordered_exon5, is_reverse_orientation)
    
    return None

def build_array_structure(exon5_list, is_reverse_orientation):
    """Build the array structure string from a list of exon5 annotations"""
    
    if not exon5_list:
        return "Unknown (no exon5 annotations found)"
    
    # Build array structure
    full_structure = "[LCR]-"
    for a in exon5_list:
        if a['type'] == 'OPN1LW_exon5':
            full_structure += "L-"
        elif a['type'] == 'OPN1MW_exon5':
            full_structure += "M-"
    
    # Remove trailing dash
    full_structure = full_structure.rstrip('-')
    
    # Add orientation indicator
    orientation_indicator = " (reverse compliment)" if is_reverse_orientation else ""
    
    # Build gene structure with consistent dashes
    gene_parts = []
    for a in exon5_list:
        if a['type'] == 'OPN1LW_exon5':
            gene_parts.append("L")
        elif a['type'] == 'OPN1MW_exon5':
            gene_parts.append("M")
    
    # Join with dashes
    gene_structure = "-".join(gene_parts)
    
    print(f"  Full array structure: {full_structure}")
    
    return f"{gene_structure}{orientation_indicator}"

if __name__ == "__main__":
    gff_file = sys.argv[1]
    sample_id = sys.argv[2]
    hap_name = sys.argv[3]
    analyze_opsin_exon5(gff_file, sample_id, hap_name)
EOF

chmod +x $OUTPUT_DIR/analyze_opsin_exon5_updated_V3.py

# Create a summary file with array structure information
if [[ "$2" == *"with-hg-size"* ]]; then 
    SUMMARY_FILE="${OUTPUT_DIR}/opsin_array_summary-with-hg-size-updated-V3.tsv"
else
    SUMMARY_FILE="${OUTPUT_DIR}/opsin_array_summary-without-hg-size-updated-V3.tsv"
fi
if [ ! -s "$SUMMARY_FILE" ]; then
    echo -e "SampleID\tSex\tHaplotype\tOPN1LW_copies\tOPN1MW_copies\tTotal_OPSIN_copies\tArray_Structure\tLCR_count" > "$SUMMARY_FILE"
fi

# Process each sample in the metadata table
while IFS=$'\t' read -r SAMPLE_ID SEX EXTRA_COLS; do
    # Skip header line, comments, or empty lines
    [[ "$SAMPLE_ID" == "SampleID" || "$SAMPLE_ID" == \#* || -z "$SAMPLE_ID" ]] && continue
    
    echo -e "\nProcessing sample: $SAMPLE_ID (Sex: $SEX)"
    
    # Create sample-specific directory
    SAMPLE_OUTPUT_DIR="${OUTPUT_DIR}/${SAMPLE_ID}_${SEX}_annotation_outputs"
    mkdir -p "$SAMPLE_OUTPUT_DIR"
    
    # Get the folder name for this sample
    SAMPLE_DIR=$(find "$ASSEMBLIES_DIR" -type d -name "*${SAMPLE_ID}*" | head -1)
    
    if [ -z "$SAMPLE_DIR" ]; then
        echo "ERROR: Could not find directory for sample $SAMPLE_ID in $ASSEMBLIES_DIR"
        echo "Searching for directories containing '$SAMPLE_ID':"
        find "$ASSEMBLIES_DIR" -type d | grep -i "$SAMPLE_ID" || echo "  No directories found"
        echo -e "$SAMPLE_ID\t$SEX\tNA\tNA\tNA\tNA\tNA\tNA" >> $SUMMARY_FILE
        continue
    fi
    
    echo "Found sample directory: $SAMPLE_DIR"
    
    # Process based on sex
    if [[ "$SEX" == "XY" ]]; then
        # For XY samples, use primary contig
        ASSEMBLY=$(find "$SAMPLE_DIR" -name "*.asm.bp.p_ctg.fa" | head -1)
        
        if [ -z "$ASSEMBLY" ]; then
            # Try alternate naming patterns
            ASSEMBLY=$(find "$SAMPLE_DIR" -name "*.p_ctg.fa" -not -name "*hap*.p_ctg.fa" | head -1)
        fi
        
        # Debug output
        echo "Searching for primary assembly in $SAMPLE_DIR:"
        find "$SAMPLE_DIR" -name "*.fa" | grep -v -E "hap1|hap2" || echo "  No primary assembly files found"
        
        if [[ -n "$ASSEMBLY" ]]; then
            echo "Found primary assembly: $ASSEMBLY"
        fi
        
        if [ -f "$ASSEMBLY" ]; then
            # Process as a single haplotype
            process_haplotype "$ASSEMBLY" "$SAMPLE_ID" "primary" "$SEX"
            
            # Extract counts and array structure for summary
            LW_COUNT=$(grep "OPN1LW:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_primary_analysis.txt" | awk '{print $2}')
            MW_COUNT=$(grep "OPN1MW:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_primary_analysis.txt" | awk '{print $2}')
            TOTAL_COUNT=$(grep "Total OPSIN genes:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_primary_analysis.txt" | awk '{print $4}')
            # Extract LCR count from the new output format
            LCR_COUNT=$(grep "Total LCR count:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_primary_analysis.txt" | awk '{print $4}')
            
            # If LCR_COUNT is empty, set to 0
            if [ -z "$LCR_COUNT" ]; then
                LCR_COUNT="0"
            fi
            
            # Get array structure - look for "Final Array Structure:" or "Using primary structure" or "Final structure:"
            ARRAY_STRUCTURE=$(grep -E "(Final Array Structure:|Using primary structure|Final structure:)" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_primary_analysis.txt" | tail -n 1 | sed 's/.*: //')
            
            # If array structure is empty, set to "Unknown"
            if [ -z "$ARRAY_STRUCTURE" ]; then
                ARRAY_STRUCTURE="Unknown"
            fi
            
            # Add to summary
            echo -e "$SAMPLE_ID\t$SEX\tprimary\t$LW_COUNT\t$MW_COUNT\t$TOTAL_COUNT\t$ARRAY_STRUCTURE\t$LCR_COUNT" >> $SUMMARY_FILE
        else
            echo "ERROR: Assembly file not found for $SAMPLE_ID"
            echo -e "$SAMPLE_ID\t$SEX\tprimary\tNA\tNA\tNA\tNA\tNA" >> $SUMMARY_FILE
        fi
    elif [[ "$SEX" == "XX" ]]; then
        # For XX samples, use both haplotypes
        HAP1_ASSEMBLY=$(find "$SAMPLE_DIR" -name "*.asm.hap1.p_ctg.fa" -o -name "*.asm.bp.hap1.p_ctg.fa" -o -name "*.hap1.p_ctg.fa" | head -1)
        HAP2_ASSEMBLY=$(find "$SAMPLE_DIR" -name "*.asm.hap2.p_ctg.fa" -o -name "*.asm.bp.hap2.p_ctg.fa" -o -name "*.hap2.p_ctg.fa" | head -1)
        
        # Debug output
        echo "Searching for haplotype files in $SAMPLE_DIR:"
        find "$SAMPLE_DIR" -name "*.fa" | grep -E "hap1|hap2" || echo "  No haplotype files found"
        
        if [[ -n "$HAP1_ASSEMBLY" ]]; then
            echo "Found haplotype 1 assembly: $HAP1_ASSEMBLY"
        fi
        
        if [[ -n "$HAP2_ASSEMBLY" ]]; then
            echo "Found haplotype 2 assembly: $HAP2_ASSEMBLY"
        fi
        
        # Process haplotype 1
        if [ -f "$HAP1_ASSEMBLY" ]; then
            process_haplotype "$HAP1_ASSEMBLY" "$SAMPLE_ID" "hap1" "$SEX"
            
            # Extract counts and array structure for summary
            LW_COUNT=$(grep "OPN1LW:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap1_analysis.txt" | awk '{print $2}')
            MW_COUNT=$(grep "OPN1MW:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap1_analysis.txt" | awk '{print $2}')
            TOTAL_COUNT=$(grep "Total OPSIN genes:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap1_analysis.txt" | awk '{print $4}')
            # Extract LCR count from the new output format
            LCR_COUNT=$(grep "Total LCR count:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap1_analysis.txt" | awk '{print $4}')
            
            # If LCR_COUNT is empty, set to 0
            if [ -z "$LCR_COUNT" ]; then
                LCR_COUNT="0"
            fi
            
            # Get array structure
            ARRAY_STRUCTURE=$(grep -E "(Final Array Structure:|Using primary structure|Final structure:)" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap1_analysis.txt" | tail -n 1 | sed 's/.*: //')
            
            # If array structure is empty, set to "Unknown"
            if [ -z "$ARRAY_STRUCTURE" ]; then
                ARRAY_STRUCTURE="Unknown"
            fi
            
            # Add to summary
            echo -e "$SAMPLE_ID\t$SEX\thap1\t$LW_COUNT\t$MW_COUNT\t$TOTAL_COUNT\t$ARRAY_STRUCTURE\t$LCR_COUNT" >> $SUMMARY_FILE
        else
            echo "ERROR: Haplotype 1 assembly file not found for $SAMPLE_ID"
            echo -e "$SAMPLE_ID\t$SEX\thap1\tNA\tNA\tNA\tNA\tNA" >> $SUMMARY_FILE
        fi
        
        # Process haplotype 2
        if [ -f "$HAP2_ASSEMBLY" ]; then
            process_haplotype "$HAP2_ASSEMBLY" "$SAMPLE_ID" "hap2" "$SEX"
            
            # Extract counts and array structure for summary
            LW_COUNT=$(grep "OPN1LW:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap2_analysis.txt" | awk '{print $2}')
            MW_COUNT=$(grep "OPN1MW:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap2_analysis.txt" | awk '{print $2}')
            TOTAL_COUNT=$(grep "Total OPSIN genes:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap2_analysis.txt" | awk '{print $4}')
            # Extract LCR count from the new output format
            LCR_COUNT=$(grep "Total LCR count:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap2_analysis.txt" | awk '{print $4}')
            
            # If LCR_COUNT is empty, set to 0
            if [ -z "$LCR_COUNT" ]; then
                LCR_COUNT="0"
            fi
            
            # Get array structure
            ARRAY_STRUCTURE=$(grep -E "(Final Array Structure:|Using primary structure|Final structure:)" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap2_analysis.txt" | tail -n 1 | sed 's/.*: //')
            
            # If array structure is empty, set to "Unknown"
            if [ -z "$ARRAY_STRUCTURE" ]; then
                ARRAY_STRUCTURE="Unknown"
            fi
            
            # Add to summary
            echo -e "$SAMPLE_ID\t$SEX\thap2\t$LW_COUNT\t$MW_COUNT\t$TOTAL_COUNT\t$ARRAY_STRUCTURE\t$LCR_COUNT" >> $SUMMARY_FILE
        else
            echo "ERROR: Haplotype 2 assembly file not found for $SAMPLE_ID"
            echo -e "$SAMPLE_ID\t$SEX\thap2\tNA\tNA\tNA\tNA\tNA" >> $SUMMARY_FILE
        fi
    else
        echo "WARNING: Unknown sex '$SEX' for sample $SAMPLE_ID. Must be XX or XY."
        echo -e "$SAMPLE_ID\t$SEX\tNA\tNA\tNA\tNA\tNA\tNA" >> $SUMMARY_FILE
    fi
    
done < "$METADATA_TABLE"

echo -e "\n=== Updated Pipeline V3 completed successfully! ==="
echo "Results are available in $OUTPUT_DIR"
echo "Summary file: $SUMMARY_FILE" 
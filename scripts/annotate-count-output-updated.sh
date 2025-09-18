#!/bin/bash

# OPSIN Gene Copy Number and Array Structure Analysis Pipeline
# Processes samples based on metadata table, focuses only on exon5 annotations
# Updated to handle multiple LCR annotations properly

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

echo "=== Starting OPSIN gene analysis pipeline (Updated for multiple LCRs) ==="

# Function to process each haplotype
process_haplotype() {
    local ASSEMBLY=$1
    local SAMPLE_ID=$2
    local HAP_NAME=$3
    local SEX=$4
    local SAMPLE_DIR="$OUTPUT_DIR/${SAMPLE_ID}_${SEX}_annotation_outputs"
    local OUTPUT_PREFIX="$SAMPLE_DIR/${SAMPLE_ID}_${HAP_NAME}"
    
    echo "Processing $SAMPLE_ID ($HAP_NAME)..."
    
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
    echo "  Extracting vulgar format and creating GFF file..."
    grep "^vulgar:" ${OUTPUT_PREFIX}.exonerate > ${OUTPUT_PREFIX}.vulgar || echo "Warning: No vulgar lines found in exonerate output"
    
    # Debug: show the number of lines in the vulgar file
    VULGAR_LINES=$(wc -l < ${OUTPUT_PREFIX}.vulgar || echo "0")
    echo "  Found $VULGAR_LINES vulgar lines for conversion to GFF"
    
    # Show the first few lines for debugging if any exist
    if [ "$VULGAR_LINES" -gt 0 ]; then
        echo "  First 3 vulgar lines:"
        head -n 3 ${OUTPUT_PREFIX}.vulgar
    fi
    
    # Convert vulgar to GFF-like format our script can process
    # Format: query_name, source, feature_type, start, end, score, strand, frame, attributes
    awk '{print $2"\t""exonerate:est2genome""\t""similarity""\t"$8"\t"$9"\t"$11"\t"$5"\t"".""\t""sequence "$2";score="$11";target="$3";target_start="$4";target_end="$5";query_start="$8";query_end="$9";strand="$5}' ${OUTPUT_PREFIX}.vulgar > ${OUTPUT_PREFIX}.gff
    
    # Debug: show the number of lines in the GFF file
    GFF_LINES=$(wc -l < ${OUTPUT_PREFIX}.gff || echo "0")
    echo "  Created $GFF_LINES GFF lines for analysis"
    
    # Show the first few GFF lines for debugging if any exist
    if [ "$GFF_LINES" -gt 0 ]; then
        echo "  First 3 GFF lines:"
        head -n 3 ${OUTPUT_PREFIX}.gff
    fi
    
    # Analyze gene copies and array structure focusing only on exon5
    echo "  Analyzing gene copies and array structure..."
    python3 $OUTPUT_DIR/analyze_opsin_exon5_updated.py ${OUTPUT_PREFIX}.gff $SAMPLE_ID $HAP_NAME > ${OUTPUT_PREFIX}_analysis.txt
    
    # Return success
    return 0
}

# Create the updated Python analysis script that handles multiple LCRs
cat > $OUTPUT_DIR/analyze_opsin_exon5_updated.py << 'EOF'
#!/usr/bin/env python3
import re
import sys
from collections import OrderedDict

def analyze_opsin_exon5(gff_file, sample_id, hap_name):
    # Initialize counters for exon5
    lw_exon5_count = 0
    mw_exon5_count = 0
    
    # Store annotation positions for array structure analysis
    annotations = []
    
    # Track the strand orientation
    strands = {}
    
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
                
                # The name is in the 1st column
                name = fields[0]
                
                # Extract position from 4th column (should be numeric)
                try:
                    position = int(fields[3]) if fields[3].isdigit() else 0
                except (ValueError, IndexError):
                    print(f"  Warning: Could not parse position from field 4")
                    position = 0
                
                # Extract strand from column 5
                strand = fields[4] if len(fields) > 4 else "+"
                strands[name] = strand
                
                # Debug info
                print(f"  Name: {name}, Position: {position}, Strand: {strand}")
                
                # Check for LCR (exact match)
                if name == 'LCR':
                    print(f"  Found LCR annotation at position {position} on strand {strand}")
                    annotations.append({
                        'type': 'LCR',
                        'position': position,
                        'strand': strand
                    })
                
                # Check for OPN1LW exon5 (exact match)
                elif name == 'OPN1LW_exon5':
                    print(f"  Found LW exon5 annotation at position {position} on strand {strand}")
                    lw_exon5_count += 1
                    annotations.append({
                        'type': 'OPN1LW_exon5',
                        'position': position,
                        'strand': strand
                    })
                
                # Check for OPN1MW exon5 (exact match)
                elif name == 'OPN1MW_exon5':
                    print(f"  Found MW exon5 annotation at position {position} on strand {strand}")
                    mw_exon5_count += 1
                    annotations.append({
                        'type': 'OPN1MW_exon5',
                        'position': position,
                        'strand': strand
                    })
                
    except Exception as e:
        print(f"Error reading GFF file: {e}")
        return 0, 0, "Error", 0
    
    # Print results
    print(f"\n=== Sample {sample_id} ({hap_name}) OPSIN Gene Analysis ===")
    
    print(f"  OPN1LW: {lw_exon5_count} copies (based on exon 5 count)")
    print(f"  OPN1MW: {mw_exon5_count} copies (based on exon 5 count)")
    print(f"  Total OPSIN genes: {lw_exon5_count + mw_exon5_count} copies")
    
    # Array structure analysis
    print("\n=== Array Structure Analysis ===")
    
    # Check if we have any annotations
    if not annotations:
        print("WARNING: No annotations found in the GFF file!")
        return 0, 0, "Unknown (No annotations)", 0
    
    # Debug: Print all annotations in the order they appear in the file
    print(f"Found {len(annotations)} total annotations:")
    for i, a in enumerate(annotations):
        print(f"  {i+1}. {a['type']} at position {a['position']} on strand {a['strand']}")
    
    # Count LCR annotations
    lcr_annotations = [a for a in annotations if a['type'] == 'LCR']
    lcr_count = len(lcr_annotations)
    print(f"\nFound {lcr_count} LCR annotations")
    
    if lcr_count == 0:
        print("WARNING: No LCR annotations found. Cannot determine array structure.")
        return lw_exon5_count, mw_exon5_count, "Unknown (No LCR found)", 0
    
    # Get all exon5 annotations
    all_exon5 = [a for a in annotations 
                if (a['type'] == 'OPN1LW_exon5' or a['type'] == 'OPN1MW_exon5')]
    
    print(f"Found {len(all_exon5)} total exon5 annotations in the haplotype")
    
    if len(all_exon5) == 0:
        print("No exon5 annotations found in the haplotype")
        return lw_exon5_count, mw_exon5_count, "Incomplete (no exon5 annotations found)", lcr_count
    
    # Determine the dominant strand orientation
    minus_strand_count = sum(1 for a in annotations if a['strand'] == '-')
    plus_strand_count = sum(1 for a in annotations if a['strand'] == '+')
    is_reverse_orientation = minus_strand_count > plus_strand_count
    
    print(f"Strand orientation: {'+' if not is_reverse_orientation else '-'} strand dominant " +
          f"({plus_strand_count} + strands, {minus_strand_count} - strands)")
    
    # Sort all annotations by position (differently depending on orientation)
    if is_reverse_orientation:
        # For reverse orientation, sort in descending order
        annotations.sort(key=lambda x: x['position'], reverse=True)
        all_exon5.sort(key=lambda x: x['position'], reverse=True)
        lcr_annotations.sort(key=lambda x: x['position'], reverse=True)
        print("Assembly is in REVERSE orientation - sorting positions in descending order")
    else:
        # For forward orientation, sort in ascending order
        annotations.sort(key=lambda x: x['position'])
        all_exon5.sort(key=lambda x: x['position'])
        lcr_annotations.sort(key=lambda x: x['position'])
        print("Assembly is in FORWARD orientation - sorting positions in ascending order")
    
    # Debug: Print sorted annotations
    print("Sorted annotations:")
    for i, a in enumerate(annotations):
        print(f"  {i+1}. {a['type']} at position {a['position']} on strand {a['strand']}")
    
    # Handle multiple LCRs by grouping genes with their closest LCR
    if lcr_count > 1:
        print(f"\nMultiple LCR annotations detected ({lcr_count}). Analyzing separate arrays:")
        
        # Group each exon5 with its closest LCR
        arrays = []
        for lcr in lcr_annotations:
            lcr_position = lcr['position']
            # Find exon5s closest to this LCR
            closest_exon5s = []
            #record exon5s that are closest to this LCR
            for exon5 in all_exon5:
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
            print(f"  Array {i+1}: LCR at {array['lcr']['position']} with {len(array['exon5s'])} genes")
            for j, exon5 in enumerate(array['exon5s']):
                print(f"    Gene {j+1}: {exon5['type']} at position {exon5['position']}")
        
        # Use the first/primary array for reporting (to avoid double-counting)
        if arrays:
            primary_array = arrays[0]
            primary_lcr = primary_array['lcr']
            primary_exon5s = primary_array['exon5s']
            
            print(f"\nUsing primary array (Array 1) for structure determination:")
            print(f"  LCR at position {primary_lcr['position']}")
            print(f"  {len(primary_exon5s)} associated genes")
            
            # Build array structure for primary array
            if len(primary_exon5s) >= 1:
                # Build array structure
                full_structure = "[LCR]-"
                for a in primary_exon5s:
                    if a['type'] == 'OPN1LW_exon5':
                        full_structure += "L-"
                    elif a['type'] == 'OPN1MW_exon5':
                        full_structure += "M-"
                
                # Remove trailing dash
                full_structure = full_structure.rstrip('-')
                
                # Add orientation indicator
                orientation_indicator = " (reverse compliment)" if is_reverse_orientation else ""
                
                print(f"Primary Array Structure:{orientation_indicator}")
                
                # Determine the first two genes (or just the first if only one exists)
                if len(primary_exon5s) >= 2:
                    first_gene = "L" if primary_exon5s[0]['type'] == 'OPN1LW_exon5' else "M"
                    second_gene = "L" if primary_exon5s[1]['type'] == 'OPN1LW_exon5' else "M"
                    print(f"  First two genes after LCR: {first_gene}-{second_gene}")
                    array_structure = f"{first_gene}-{second_gene}"
                else:
                    first_gene = "L" if primary_exon5s[0]['type'] == 'OPN1LW_exon5' else "M"
                    print(f"  First gene after LCR: {first_gene}")
                    array_structure = first_gene
                
                print(f"  Full primary array structure: {full_structure}")
                
                return lw_exon5_count, mw_exon5_count, f"{array_structure}{orientation_indicator}", lcr_count
        
    else:
        # Single LCR - use original logic
        first_lcr = lcr_annotations[0]
        lcr_position = first_lcr['position']
        
        if is_reverse_orientation:
            print(f"LCR found at position {lcr_position} (reverse orientation)")
        else:
            print(f"LCR found at position {lcr_position} (forward orientation)")
        
        # Find the closest exon5 to the LCR based on absolute distance
        closest_exon5 = min(all_exon5, key=lambda x: abs(x['position'] - lcr_position))
        closest_index = all_exon5.index(closest_exon5)
        
        print(f"Closest exon5 to LCR: {closest_exon5['type']} at position {closest_exon5['position']}")
        
        # Reorder all_exon5 to put the closest one first
        reordered_exon5 = [closest_exon5] + [e for e in all_exon5 if e != closest_exon5]
        
        # Print the reordered annotations
        print("Reordered exon5 annotations:")
        for i, ann in enumerate(reordered_exon5):
            print(f"  {i+1}. {ann['type']} at position {ann['position']} on strand {ann['strand']}")
        
        # Determine array structure based on reordered exon5 annotations
        if len(reordered_exon5) >= 1:
            # Build array structure
            full_structure = "[LCR]-"
            for a in reordered_exon5:
                if a['type'] == 'OPN1LW_exon5':
                    full_structure += "L-"
                elif a['type'] == 'OPN1MW_exon5':
                    full_structure += "M-"
            
            # Remove trailing dash
            full_structure = full_structure.rstrip('-')
            
            # Add orientation indicator
            orientation_indicator = " (reverse compliment)" if is_reverse_orientation else ""
            
            print(f"Array Structure:{orientation_indicator}")
            
            # Determine the first two genes (or just the first if only one exists)
            if len(reordered_exon5) >= 2:
                first_gene = "L" if reordered_exon5[0]['type'] == 'OPN1LW_exon5' else "M"
                second_gene = "L" if reordered_exon5[1]['type'] == 'OPN1LW_exon5' else "M"
                print(f"  First two genes after LCR: {first_gene}-{second_gene}")
                array_structure = f"{first_gene}-{second_gene}"
            else:
                first_gene = "L" if reordered_exon5[0]['type'] == 'OPN1LW_exon5' else "M"
                print(f"  First gene after LCR: {first_gene}")
                array_structure = first_gene
            
            print(f"  Full array structure: {full_structure}")
            
            return lw_exon5_count, mw_exon5_count, f"{array_structure}{orientation_indicator}", lcr_count
    
    print("No exon5 annotations found to determine array structure")
    return lw_exon5_count, mw_exon5_count, "Incomplete (no exon5 annotations found)", lcr_count

if __name__ == "__main__":
    gff_file = sys.argv[1]
    sample_id = sys.argv[2]
    hap_name = sys.argv[3]
    analyze_opsin_exon5(gff_file, sample_id, hap_name)
EOF

chmod +x $OUTPUT_DIR/analyze_opsin_exon5_updated.py

# Create a summary file with array structure information
if [[ "$2" == *"with-hg-size"* ]]; then 
    SUMMARY_FILE="${OUTPUT_DIR}/opsin_array_summary-with-hg-size-updated.tsv"
else
    SUMMARY_FILE="${OUTPUT_DIR}/opsin_array_summary-without-hg-size-updated.tsv"
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
            LCR_COUNT=$(grep "Found [0-9]* LCR annotations" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_primary_analysis.txt" | awk '{print $2}')
            
            # Get array structure - extract the "First two genes after LCR:" line
            ARRAY_STRUCTURE=$(grep -A 1 "First two genes after LCR:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_primary_analysis.txt" | tail -n 1 | awk '{$1=""; $2=""; print $0}' | sed 's/^[ \t]*//')
            
            # If array structure is empty, set to "Unknown"
            if [ -z "$ARRAY_STRUCTURE" ]; then
                ARRAY_STRUCTURE="Unknown"
            fi
            
            # Check if this is a reverse orientation
            if grep -q "reverse compliment" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_primary_analysis.txt"; then
                ORIENTATION="(reverse compliment)"
                # Add orientation info only if it's not already there
                if [[ "$ARRAY_STRUCTURE" != *"reverse compliment"* ]]; then
                    ARRAY_STRUCTURE="$ARRAY_STRUCTURE $ORIENTATION"
                fi
            fi
            
            # Add to summary
            echo -e "$SAMPLE_ID\t$SEX\tprimary\t$LW_COUNT\t$MW_COUNT\t$TOTAL_COUNT\tstructure: [LCR]-$ARRAY_STRUCTURE\t$LCR_COUNT" >> $SUMMARY_FILE
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
            LCR_COUNT=$(grep "Found [0-9]* LCR annotations" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap1_analysis.txt" | awk '{print $2}')
            
            # Get array structure - extract the "First two genes after LCR:" line
            ARRAY_STRUCTURE=$(grep -A 1 "First two genes after LCR:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap1_analysis.txt" | tail -n 1 | awk '{$1=""; $2=""; print $0}' | sed 's/^[ \t]*//')
            
            # If array structure is empty, set to "Unknown"
            if [ -z "$ARRAY_STRUCTURE" ]; then
                ARRAY_STRUCTURE="Unknown"
            fi
            
            # Check if this is a reverse orientation
            if grep -q "reverse compliment" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap1_analysis.txt"; then
                ORIENTATION="(reverse compliment)"
                # Add orientation info only if it's not already there
                if [[ "$ARRAY_STRUCTURE" != *"reverse compliment"* ]]; then
                    ARRAY_STRUCTURE="$ARRAY_STRUCTURE $ORIENTATION"
                fi
            fi
            
            # Add to summary
            echo -e "$SAMPLE_ID\t$SEX\thap1\t$LW_COUNT\t$MW_COUNT\t$TOTAL_COUNT\tstructure: [LCR]-$ARRAY_STRUCTURE\t$LCR_COUNT" >> $SUMMARY_FILE
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
            LCR_COUNT=$(grep "Found [0-9]* LCR annotations" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap2_analysis.txt" | awk '{print $2}')
            
            # Get array structure - extract the "First two genes after LCR:" line
            ARRAY_STRUCTURE=$(grep -A 1 "First two genes after LCR:" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap2_analysis.txt" | tail -n 1 | awk '{$1=""; $2=""; print $0}' | sed 's/^[ \t]*//')
            
            # If array structure is empty, set to "Unknown"
            if [ -z "$ARRAY_STRUCTURE" ]; then
                ARRAY_STRUCTURE="Unknown"
            fi
            
            # Check if this is a reverse orientation
            if grep -q "reverse compliment" "$SAMPLE_OUTPUT_DIR/${SAMPLE_ID}_hap2_analysis.txt"; then
                ORIENTATION="(reverse compliment)"
                # Add orientation info only if it's not already there
                if [[ "$ARRAY_STRUCTURE" != *"reverse compliment"* ]]; then
                    ARRAY_STRUCTURE="$ARRAY_STRUCTURE $ORIENTATION"
                fi
            fi
            
            # Add to summary
            echo -e "$SAMPLE_ID\t$SEX\thap2\t$LW_COUNT\t$MW_COUNT\t$TOTAL_COUNT\tstructure: [LCR]-$ARRAY_STRUCTURE\t$LCR_COUNT" >> $SUMMARY_FILE
        else
            echo "ERROR: Haplotype 2 assembly file not found for $SAMPLE_ID"
            echo -e "$SAMPLE_ID\t$SEX\thap2\tNA\tNA\tNA\tNA\tNA" >> $SUMMARY_FILE
        fi
    else
        echo "WARNING: Unknown sex '$SEX' for sample $SAMPLE_ID. Must be XX or XY."
        echo -e "$SAMPLE_ID\t$SEX\tNA\tNA\tNA\tNA\tNA\tNA" >> $SUMMARY_FILE
    fi
    
done < "$METADATA_TABLE"

echo -e "\n=== Updated Pipeline completed successfully! ==="
echo "Results are available in $OUTPUT_DIR"
echo "Summary file: $SUMMARY_FILE" 
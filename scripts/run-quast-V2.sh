#!/bin/bash

# Usage: ./run_quast.sh <Assembly_dir> <sample_list> <ref_genome>
# Example: ./run_quast.sh assemblies/ samples.txt hg38

usage() {
    echo "Usage: $0 -a <Assembly_dir> -s <sample_list> -r <ref_genome> -o <output_dir_suffix> -b <output_base_dir>"
    echo "Example: $0 -a Assembly_outputs/Assemblies-from-hifiasm-with-hg-size -s samples.txt -r hg38 -o with-hg-size-all-reads -b Assembly_outputs/Quast-outputs"
    exit 1
}
while getopts ":a:s:r:o:b:" opt; do
    case ${opt} in
        a) Assembly_dir=$OPTARG ;;
        s) sample_list=$OPTARG ;;
        r) ref_genome=$OPTARG ;;
        o) output_dir_suffix=$OPTARG ;;
        b) out_base_dir=$OPTARG ;;
        h) usage ;;
        \?) echo "Invalid option: -$OPTARG" >&2
            usage ;;
    esac
done

# Check if all arguments are provided
if [ -z "$Assembly_dir" ] || [ -z "$sample_list" ] || [ -z "$ref_genome" ] || [ -z "$output_dir_suffix" ] || [ -z "$out_base_dir" ]; then
    echo "Usage: $0 -a <Assembly_dir> -s <sample_list> -r <ref_genome> -o <output_dir_suffix> -b <output_base_dir>"
    echo "Example: $0 -a Assembly_outputs/Assemblies-from-hifiasm-with-hg-size -s samples.txt -r hg38 -o with-hg-size-all-reads -b Assembly_outputs/Quast-outputs"
    exit 1
fi

# Check if files/directories exist
if [ ! -d "$Assembly_dir" ]; then
    echo "Error: Assembly directory '$Assembly_dir' does not exist"
    exit 1
fi

if [ ! -f "$sample_list" ]; then
    echo "Error: Sample list file '$sample_list' does not exist"
    exit 1
fi

hg38_filtered_reference="hg38-2MB-assembly-region.fa"
chm13_filtered_reference="chm13-2MB-assembly-region.fa"

if [ "$ref_genome" = "hg38" ]; then
    reference_genome=$hg38_filtered_reference
elif [ "$ref_genome" = "chm13" ]; then
    reference_genome=$chm13_filtered_reference
else
    echo "Error: Invalid reference genome '$ref_genome'. Use 'hg38' or 'chm13'"
    exit 1
fi

# Check if reference genome file exists
if [ ! -f "$reference_genome" ]; then
    echo "Error: Reference genome file '$reference_genome' does not exist"
    exit 1
fi

echo "Starting QUAST analysis with reference: $reference_genome"
echo "Assembly directory: $Assembly_dir"
echo "Sample list: $sample_list"
echo ""

processed_count=0
error_count=0

# Initialize combined transposed report output
combined_report="${out_base_dir}/combined-transposed_report-${output_dir_suffix}-${ref_genome}.tsv"
header_written=false
rm -f "$combined_report"


while read -r sample sex; do
    # Skip comment lines and empty lines
    [[ "$sample" =~ ^#.* ]] && continue
    [ -z "$sample" ] && continue
    
    echo "Processing sample: $sample (Sex: $sex)"
    
    output_dir="${out_base_dir}/${sample}-${sex}-${output_dir_suffix}-${ref_genome}-quast-results"
    mkdir -p "$output_dir"

    if [ "$sex" = "XY" ]; then
        primary_assembly=$(find "$Assembly_dir" -name "${sample}*.bp.p_ctg.fa")
        
        if [ -z "$primary_assembly" ]; then
            echo "Warning: No primary assembly found for XY sample $sample"
            ((error_count++))
            continue
        fi
        
        echo "  Found primary assembly: $primary_assembly"
        
        quast --reference "$reference_genome" \
              --output-dir "$output_dir" \
              --threads 10 \
              --min-contig 1000 \
              --eukaryote \
              --fragmented \
              --labels "${sample}_primary" \
              "$primary_assembly"
              
    elif [ "$sex" = "XX" ]; then
        hap1_assembly=$(find "$Assembly_dir" -name "${sample}*asm.hap1.p_ctg.fa")
        hap2_assembly=$(find "$Assembly_dir" -name "${sample}*asm.hap2.p_ctg.fa")
        
        if [ -z "$hap1_assembly" ] || [ -z "$hap2_assembly" ]; then
            echo "Warning: Missing haplotype assembly for XX sample $sample"
            echo "  hap1: $hap1_assembly"
            echo "  hap2: $hap2_assembly"
            ((error_count++))
            continue
        fi
        
        echo "  Found hap1 assembly: $hap1_assembly"
        echo "  Found hap2 assembly: $hap2_assembly"
        
        quast --reference "$reference_genome" \
              --output-dir "$output_dir" \
              --threads 10 \
              --min-contig 1000 \
              --eukaryote \
              --fragmented \
              --labels "${sample}_hap1,${sample}_hap2" \
              "$hap1_assembly" "$hap2_assembly"
 
    else
        echo "Warning: Unknown sex '$sex' for sample $sample. Skipping."
        ((error_count++))
        continue
    fi
    
    if [ $? -eq 0 ]; then
        echo "  QUAST completed successfully for $sample"
        ((processed_count++))
        
        # Append this sample's transposed report to the combined file
        transposed_file="$output_dir/transposed_report.tsv"
        if [ -f "$transposed_file" ] && [ -s "$transposed_file" ]; then
            if [ "$header_written" = false ]; then
                # Write header and all data rows from the first available sample
                cat "$transposed_file" > "$combined_report"
                header_written=true
            else
                # Append only data rows (skip the header line)
                tail -n +2 "$transposed_file" >> "$combined_report"
            fi
            echo "  Appended transposed_report.tsv to combined file"
        else
            echo "  Warning: transposed_report.tsv not found or empty for $sample"
        fi
    else
        echo "  Error: QUAST failed for $sample"
        ((error_count++))
    fi
    echo ""
    
done < "$sample_list"

echo "QUAST analysis complete!"
echo "Successfully processed: $processed_count samples"
echo "Errors encountered: $error_count samples"

if [ -f "$combined_report" ]; then
    echo "Combined transposed report written to: $combined_report"
else
    echo "No combined transposed report was generated."
fi
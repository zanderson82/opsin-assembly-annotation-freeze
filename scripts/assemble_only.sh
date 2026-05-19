#!/bin/bash
set +e  # Allow script to continue even if a command exits with non-zero status

# Initialize counters for reporting
total_samples=0
processed_samples=0
failed_samples=0
samtools_flag="-F 0"  # Default samtools flag

# Check if arguments are provided
if [ $# -lt 3 ]; then
    echo "Usage: $0 <samples_file> <genomic_region> <region_name> <bam_dir> [options]"
    echo "Example: $0 samples_with_sex.txt chrX:153121316-155216212 (hg38) chrX:151389254-153479422 (T2T) OPSIN-PLUS-1MB /n/zanderson/1KGP_alignments/align-minimap2-2.24-hg38"
    echo ""
    echo "The samples file should have two columns: Sample_ID Sex"
    echo "Example:"
    echo "HG00096 XX"
    echo "HG00097 XY"
    echo ""
    echo "Options:"
    echo "  --size [true|false]    Whether to use the --hg-size flag with hifiasm"
    echo "                         Default: false if not specified"
    echo "                         Use just '--size' (no value) to enable"
    echo "  --output-dir <path>    Directory to store all sample output directories"
    echo "                         Default: current directory"
    echo "  --samtools-flag <flag> Flag to pass to samtools view"
    echo "                         Default: '-F 0' (primary alignments)"
    echo "                         Example: '--samtools-flag \"-F 0x900\"'"
    echo "  --use-subdirs          Look for BAM files in FIRST_100 and 100_PLUS subdirectories"
    echo "                         Default: false (look directly in bam_dir)"
    echo "  --input-suffix <suffix> Default: .phased.bam"
    echo "                         Example: '--input-suffix .fastq'"
    exit 1
fi

samples_file="$1"
region="$2"
region_name="$3"
bam_dir="$4"
HIFIASM="../hifiasm/hifiasm"

# Default to NOT using --hg-size
use_size=false
# Default output directory is current directory
output_base_dir="$(pwd)"
# Default to NOT using subdirectories
use_subdirs=false
# Default input suffix is .phased.bam
input_suffix=".phased.bam"
# Check for optional flags in arguments
shift 3  # Skip the first three required arguments
while [ $# -gt 0 ]; do
    case "$1" in
        --size)
            # Check if the next argument is a flag (starts with --) or doesn't exist
            if [[ "$2" == --* ]] || [ -z "$2" ]; then
                use_size=true  # If --size is provided without a value, assume true
                echo "Will use --hg-size flag with hifiasm when possible"
                shift 1
            else
                if [ "$2" = "true" ]; then
                    use_size=true
                    echo "Will use --hg-size flag with hifiasm when possible"
                else
                    use_size=false
                    echo "Will NOT use --hg-size flag with hifiasm"
                fi
                shift 2
            fi
            ;;
        --output-dir)
            if [ -n "$2" ]; then
                output_base_dir="$2"
                # Create the output directory if it doesn't exist
                mkdir -p "$output_base_dir"
                echo "Will use $output_base_dir as the base output directory"
                shift 2
            else
                echo "ERROR: No directory specified for --output-dir"
                exit 1
            fi
            ;;
        --samtools-flag)
            if [ -n "$2" ]; then
                samtools_flag="$2"
                echo "Will use samtools flag: $samtools_flag"
                shift 2
            else
                echo "ERROR: No flag specified for --samtools-flag"
                exit 1
            fi
            ;;
        --use-subdirs)
            use_subdirs=true
            echo "Will look for BAM files in FIRST_100 and 100_PLUS subdirectories"
            shift 1
            ;;
        --input-suffix)
            if [ -n "$2" ]; then
                input_suffix="$2"
                echo "Will use $input_suffix as the input suffix"
                shift 2
            fi
            ;;
        *)
            echo "Unknown parameter: $1"
            shift
            ;;
    esac
done

# Check if samples file exists
if [ ! -f "$samples_file" ]; then
    echo "ERROR: Sample file $samples_file not found"
    exit 1
fi

# Calculate region size for --hg-size flag
region_size=0
region_size_with_unit=""
if [[ "$region" =~ ([^:]+):([0-9]+)-([0-9]+) ]]; then
    chr="${BASH_REMATCH[1]}"
    start="${BASH_REMATCH[2]}"
    end="${BASH_REMATCH[3]}"
    region_size=$((end - start))
    
    # Format the region size with the appropriate unit (k/m/g) for hifiasm
    # Using integer values and rounding to the nearest integer
    if [ $region_size -ge 1000000000 ]; then
        # Size in gigabases (integer)
        size_in_g=$(( (region_size + 500000000) / 1000000000 ))
        region_size_with_unit="${size_in_g}g"
    elif [ $region_size -ge 1000000 ]; then
        # Size in megabases (integer)
        size_in_m=$(( (region_size + 500000) / 1000000 ))
        region_size_with_unit="${size_in_m}m"
    else
        # Size in kilobases (integer)
        size_in_k=$(( (region_size + 500) / 1000 ))
        # Ensure at least 1k
        if [ $size_in_k -lt 1 ]; then
            size_in_k=1
        fi
        region_size_with_unit="${size_in_k}k"
    fi
    
    echo "Calculated region size: $region_size bp (${region_size_with_unit}) from $start to $end on $chr"
else
    echo "WARNING: Could not parse region format for size calculation. Expected format: chrX:START-END"
    echo "Will not use --hg-size flag."
fi

# Define log file
# Create logs directory if it doesn't exist
mkdir -p "${output_base_dir}/logs"
log_file="${output_base_dir}/logs/assembly_$(date +%Y%m%d_%H%M%S).log"
echo "Starting assembly for region $region ($region_name)" | tee -a "$log_file"
echo "Using samples file: $samples_file" | tee -a "$log_file"
echo "Output base directory: $output_base_dir" | tee -a "$log_file"
echo "Samtools flag: $samtools_flag" | tee -a "$log_file"
echo "Use subdirectories: $use_subdirs" | tee -a "$log_file"
if [ -n "$region_size_with_unit" ] && [ "$use_size" = true ]; then
    echo "Region size: $region_size bp (${region_size_with_unit} for hifiasm)" | tee -a "$log_file"
else
    echo "Will not use --hg-size flag with hifiasm" | tee -a "$log_file"
fi
echo "Logging to $log_file"

# Read the samples file into arrays
declare -a sample_ids
declare -a sample_sexes

# Count lines in sample file
sample_count=$(wc -l < "$samples_file")
echo "Found $sample_count entries in $samples_file" | tee -a "$log_file"

# Read sample IDs and sexes from the file
while read -r id sex; do
    # Skip lines starting with # (headers)
    [[ "$id" =~ ^#.* ]] && continue
    
    # Skip empty lines
    [ -z "$id" ] && continue
    
    # Add to arrays
    sample_ids+=("$id")
    sample_sexes+=("$sex")
done < "$samples_file"

echo "Loaded ${#sample_ids[@]} samples from $samples_file" | tee -a "$log_file"

# Process each sample
for i in "${!sample_ids[@]}"; do
    sample="${sample_ids[$i]}"
    sex="${sample_sexes[$i]}"
    
    ((total_samples++))
    echo "[$total_samples/${#sample_ids[@]}] Processing sample: $sample (Sex: $sex)" | tee -a "$log_file"
    sample_success=false
    
    # Set HiFiasm options based on sex
    hifiasm_opts="--ont -t10"
    if [ "$sex" = "XY" ]; then
        echo "  Using -l0 flag for XY sample" | tee -a "$log_file"
        hifiasm_opts="$hifiasm_opts -l0"
    fi
    
    # Add region size as --hg-size if available and if use_size is true
    if [ -n "$region_size_with_unit" ] && [ "$use_size" = true ]; then
        echo "  Using --hg-size $region_size_with_unit" | tee -a "$log_file"
        hifiasm_opts="$hifiasm_opts --hg-size $region_size_with_unit"
    fi
    
    # Determine search directories based on use_subdirs flag
    if [ "$use_subdirs" = true ]; then
        search_dirs=("FIRST_100" "100_PLUS")
        echo "  Searching in subdirectories: FIRST_100 and 100_PLUS" | tee -a "$log_file"
    else
        search_dirs=(".")
        echo "  Searching directly in bam_dir: $bam_dir" | tee -a "$log_file"
    fi
    
    for dir in "${search_dirs[@]}"; do
        # Set the search path based on whether we're using subdirs or not
        if [ "$use_subdirs" = true ]; then
            search_path="$bam_dir/$dir"
        else
            search_path="$bam_dir"
        fi
        
        # Find the BAM files in the specific directories
        #candidate files based on provided suffix
        input_files=$(find "$search_path" -name "*${sample}*${input_suffix}" 2>/dev/null)
        if [ -z "$input_files" ]; then
            if [ "$use_subdirs" = true ]; then
                echo "  No input files (suffix: $input_suffix) found for $sample in $dir" | tee -a "$log_file"
            else
                echo "  No input files (suffix: $input_suffix) found for $sample in $bam_dir" | tee -a "$log_file"
            fi
            continue
        fi
        
        if [ "$use_subdirs" = true ]; then
            echo "  Found $(echo "$input_files" | wc -l) input files (suffix: $input_suffix) for $sample in $dir" | tee -a "$log_file"
        else
            echo "  Found $(echo "$input_files" | wc -l) input files (suffix: $input_suffix) for $sample" | tee -a "$log_file"
        fi
        
        for input_file in $input_files; do
            # Extract sample name from BAM filename
            #filename=$(basename "$input_file" .phased.bam)
            filename=$(basename "$input_file" "$input_suffix")
            # Normalize: strip any trailing dot (can occur if input_suffix lacks leading '.')
            filename="${filename%.}"
            output_dir="${output_base_dir}/${filename}_${sex}_hifiasm_outputs"
            mkdir -p "$output_dir"
            
            # Determine reads FASTQ depending on input type
            reads_fastq="${output_dir}/${filename}.${sex}.${region_name}.fastq"
            if [ "$input_suffix" = ".fastq" ]; then
                echo "  Using existing FASTQ input: $input_file" | tee -a "$log_file"
                # Create a symlink with the expected name for downstream steps if it doesn't exist
                if [ ! -e "$reads_fastq" ]; then
                    ln -s "$input_file" "$reads_fastq"
                fi
            else
                echo "  Extracting reads from the bam file: $input_file for region $region" | tee -a "$log_file"
                if ! samtools view -@ 15 -b "$samtools_flag" "$input_file" "$region" > "${output_dir}/${filename}.${sex}.${region_name}.temp.bam"; then
                    echo "  ERROR: Failed to extract reads from $input_file" | tee -a "$log_file"
                    continue
                fi
                read_count=$(samtools view -c "${output_dir}/${filename}.${sex}.${region_name}.temp.bam")
                if [ "$read_count" -eq 0 ]; then
                    echo "  WARNING: No reads found in region $region for $filename" | tee -a "$log_file"
                    rm "${output_dir}/${filename}.${sex}.${region_name}.temp.bam"
                    continue
                fi
                echo "  Converting extracted reads to FASTQ" | tee -a "$log_file"
                if ! samtools fastq "${output_dir}/${filename}.${sex}.${region_name}.temp.bam" > "$reads_fastq"; then
                    echo "  ERROR: Failed to convert BAM to FASTQ for $filename" | tee -a "$log_file"
                    rm "${output_dir}/${filename}.${sex}.${region_name}.temp.bam"
                    continue
                fi
                rm "${output_dir}/${filename}.${sex}.${region_name}.temp.bam"
            fi

                # Save current directory before changing to output_dir
                current_dir=$(pwd)
                cd "$output_dir" || { 
                    echo "  ERROR: Failed to change directory to $output_dir" | tee -a "$current_dir/$log_file"
                    continue
                }
                
                # Run hifiasm assembly with appropriate options based on sex
                echo "  Running hifiasm assembly for ${filename} (${sex})" | tee -a "$log_file"
                # Use basename here since we're already inside output_dir
                reads_fastq_basename=$(basename "$reads_fastq")
                echo "  Command: $HIFIASM -o ${filename}.${sex}.${region_name}.asm $hifiasm_opts $reads_fastq_basename" | tee -a "$log_file"
                
                if ! $HIFIASM -o "${filename}.${sex}.${region_name}.asm" $hifiasm_opts "$reads_fastq_basename" 2> "${filename}.${sex}.${region_name}.asm.log"; then
                    echo "  WARNING: hifiasm may have encountered issues, check ${filename}.${sex}.${region_name}.asm.log" | tee -a "$log_file"
                    # Continue anyway as hifiasm might still produce usable output
                fi
                
                # Define assembly files - note that output paths now include sex
                hap1_gfa="${filename}.${sex}.${region_name}.asm.bp.hap1.p_ctg.gfa"
                hap2_gfa="${filename}.${sex}.${region_name}.asm.bp.hap2.p_ctg.gfa"
                primary_gfa="${filename}.${sex}.${region_name}.asm.bp.p_ctg.gfa" # Non-haplotype specific assembly
                
                # Check if haplotype files exist and convert to FASTA
                echo "  Converting assembly GFA files to FASTA" | tee -a "$log_file"
                
                assembly_found=false
                
                # Check for hap1
                if [ -f "$hap1_gfa" ] && [ -s "$hap1_gfa" ]; then
                    echo "    Found haplotype 1 assembly for ${filename} (${sex})" | tee -a "$log_file"
                    awk '/^S/{print ">"$2"\n"$3}' "$hap1_gfa" > "${filename}.${sex}.${region_name}.asm.hap1.p_ctg.fa"
                    samtools faidx "${filename}.${sex}.${region_name}.asm.hap1.p_ctg.fa"
                    assembly_found=true
                else
                    echo "    WARNING: No haplotype 1 assembly found for ${filename} (${sex})" | tee -a "$log_file"
                fi
                
                # Check for hap2
                if [ -f "$hap2_gfa" ] && [ -s "$hap2_gfa" ]; then
                    echo "    Found haplotype 2 assembly for ${filename} (${sex})" | tee -a "$log_file"
                    awk '/^S/{print ">"$2"\n"$3}' "$hap2_gfa" > "${filename}.${sex}.${region_name}.asm.hap2.p_ctg.fa"
                    samtools faidx "${filename}.${sex}.${region_name}.asm.hap2.p_ctg.fa"
                    assembly_found=true
                else
                    echo "    WARNING: No haplotype 2 assembly found for ${filename} (${sex})" | tee -a "$log_file"
                fi
                
                # Always check for primary assembly and process it
                if [ -f "$primary_gfa" ] && [ -s "$primary_gfa" ]; then
                    echo "    Processing primary assembly for ${filename} (${sex})" | tee -a "$log_file"
                    awk '/^S/{print ">"$2"\n"$3}' "$primary_gfa" > "${filename}.${sex}.${region_name}.asm.bp.p_ctg.fa"
                    samtools faidx "${filename}.${sex}.${region_name}.asm.bp.p_ctg.fa"
                    assembly_found=true
                else
                    echo "    WARNING: No primary assembly found for ${filename} (${sex})" | tee -a "$log_file"
                fi
                
                if [ "$assembly_found" = true ]; then
                    echo "  Finished assembly for ${filename} (${sex})" | tee -a "$log_file"
                    ((processed_samples++))
                    sample_success=true
                else
                    echo "  ERROR: No assemblies found for ${filename} (${sex})" | tee -a "$log_file"
                fi
                
                cd "$current_dir" || echo "Failed to return to $current_dir"
            
        done
    done
    
    if [ "$sample_success" = false ]; then
        echo "WARNING: Failed to process sample $sample (${sex})" | tee -a "$log_file"
        ((failed_samples++))
    fi
done

# Summary
echo "==== SUMMARY ====" | tee -a "$log_file"
echo "Total samples: $total_samples" | tee -a "$log_file"
echo "Successfully processed: $processed_samples" | tee -a "$log_file"
echo "Failed samples: $failed_samples" | tee -a "$log_file"
echo "Assembly completed." | tee -a "$log_file"
exit 0 
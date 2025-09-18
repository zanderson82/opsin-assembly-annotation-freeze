#!/bin/bash
set -e  # Exit immediately if a command exits with non-zero status

# Check if arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <genomic_region> <region_name>"
    exit 1
fi

region="$1"
region_name="$2"
HIFIASM="/n/zanderson/OPSIN-carrier-screen/manual-assembly-and-annotate/hifiasm/hifiasm"
OPSIN_FASTA="./OPSIN-gene-sequences.fasta"  # Path to your OPSIN genes
declare -a samples
mapfile -t samples < samples.txt

for sample in "${samples[@]}"; do 
    for dir in "FIRST_100" "100_PLUS"; do
        # Find the BAM files in the specific directories
        bam_files=$(find /waldo/1KGP_alignments/align-minimap2-2.24-hg38/$dir -name "*${sample}*.phased.bam" 2>/dev/null)
        
        for bam_file in $bam_files; do
            # Extract sample name from BAM filename
            filename=$(basename "$bam_file" .phased.bam)
            output_dir="${filename}_hifiasm_outputs"
            mkdir -p "$output_dir"
            
            # First extract the reads to a temporary SAM file
            samtools view -@ 15 -b -F 0 "$bam_file" "$region" > "${output_dir}/${filename}.${region_name}.temp.bam"
            # Then convert to FASTQ
            samtools fastq "${output_dir}/${filename}.${region_name}.temp.bam" > "${output_dir}/${filename}.${region_name}.fastq"
            # Remove the temporary file
            rm "${output_dir}/${filename}.${region_name}.temp.bam"
            
            # Save current directory
            current_dir=$(pwd)
            cd "$output_dir" || { echo "Failed to change directory to $output_dir"; exit 1; }
            
            # Run hifiasm assembly
            $HIFIASM -o "${filename}.${region_name}.asm" --ont -t15 "${filename}.${region_name}.fastq" 2> "${filename}.${region_name}.asm.log"
            
            # Define assembly files
            hap1_gfa="${filename}.${region_name}.asm.bp.hap1.p_ctg.gfa"
            hap2_gfa="${filename}.${region_name}.asm.bp.hap2.p_ctg.gfa"
            primary_gfa="${filename}.${region_name}.asm.bp.p_ctg.gfa" # Non-haplotype specific assembly
            
            # Check if haplotype files exist and convert to FASTA
            hap1_exists=false
            hap2_exists=false
            primary_exists=false
            
            # Check for hap1
            if [ -f "$hap1_gfa" ] && [ -s "$hap1_gfa" ]; then
                echo "Found haplotype 1 assembly for ${filename}"
                awk '/^S/{print ">"$2"\n"$3}' "$hap1_gfa" > "${filename}.${region_name}.asm.hap1.p_ctg.fa"
                hap1_exists=true
            else
                echo "WARNING: No haplotype 1 assembly found for ${filename}"
            fi
            
            # Check for hap2
            if [ -f "$hap2_gfa" ] && [ -s "$hap2_gfa" ]; then
                echo "Found haplotype 2 assembly for ${filename}"
                awk '/^S/{print ">"$2"\n"$3}' "$hap2_gfa" > "${filename}.${region_name}.asm.hap2.p_ctg.fa"
                hap2_exists=true
            else
                echo "WARNING: No haplotype 2 assembly found for ${filename}"
            fi
            
            # Always check for primary assembly and process it
            if [ -f "$primary_gfa" ] && [ -s "$primary_gfa" ]; then
                echo "Processing primary assembly for ${filename}"
                awk '/^S/{print ">"$2"\n"$3}' "$primary_gfa" > "${filename}.${region_name}.asm.bp.p_ctg.fa"
                primary_exists=true
            else
                echo "WARNING: No primary assembly found for ${filename}"
            fi
            
            # Run minimap2 to align OPSIN genes to available assemblies
            if [ "$hap1_exists" = true ]; then
                echo "Aligning OPSIN to haplotype 1 assembly..."
                samtools faidx "${filename}.${region_name}.asm.hap1.p_ctg.fa"
                minimap2 -t8 -ax map-ont --secondary=yes -N10 "${filename}.${region_name}.asm.hap1.p_ctg.fa" "$current_dir/$OPSIN_FASTA" > "${filename}.${region_name}.opsin_alignment_hap1.sam"
                
                # Convert to BAM, sort and index for IGV viewing
                samtools view -bS "${filename}.${region_name}.opsin_alignment_hap1.sam" | \
                samtools sort > "${filename}.${region_name}.opsin_alignment_hap1.bam"
                samtools index "${filename}.${region_name}.opsin_alignment_hap1.bam"
                echo "Haplotype 1 alignment complete."
            fi
            
            if [ "$hap2_exists" = true ]; then
                echo "Aligning OPSIN to haplotype 2 assembly..."
                samtools faidx "${filename}.${region_name}.asm.hap2.p_ctg.fa"
                minimap2 -t8 -ax map-ont --secondary=yes -N10 "${filename}.${region_name}.asm.hap2.p_ctg.fa" "$current_dir/$OPSIN_FASTA" > "${filename}.${region_name}.opsin_alignment_hap2.sam"
                
                # Convert to BAM, sort and index for IGV viewing
                samtools view -bS "${filename}.${region_name}.opsin_alignment_hap2.sam" | \
                samtools sort > "${filename}.${region_name}.opsin_alignment_hap2.bam"
                samtools index "${filename}.${region_name}.opsin_alignment_hap2.bam"
                echo "Haplotype 2 alignment complete."
            fi
            
            # Always align to the primary assembly if it exists
            if [ "$primary_exists" = true ]; then
                echo "Aligning OPSIN to primary assembly..."
                samtools faidx "${filename}.${region_name}.asm.bp.p_ctg.fa"
                minimap2 -t8 -ax map-ont --secondary=yes -N10 "${filename}.${region_name}.asm.bp.p_ctg.fa" "$current_dir/$OPSIN_FASTA" > "${filename}.${region_name}.opsin_alignment_primary.sam"
                
                # Convert to BAM, sort and index for IGV viewing
                samtools view -bS "${filename}.${region_name}.opsin_alignment_primary.sam" | \
                samtools sort > "${filename}.${region_name}.opsin_alignment_primary.bam"
                samtools index "${filename}.${region_name}.opsin_alignment_primary.bam"
                echo "Primary assembly alignment complete."
            fi
            
            if [ "$hap1_exists" = false ] && [ "$hap2_exists" = false ] && [ "$primary_exists" = false ]; then
                echo "ERROR: No assemblies found for ${filename}. Check hifiasm output."
            fi
            
            cd "$current_dir" || { echo "Failed to return to $current_dir"; exit 1; }
        done
    done
done

echo "Processing completed successfully."
exit 0
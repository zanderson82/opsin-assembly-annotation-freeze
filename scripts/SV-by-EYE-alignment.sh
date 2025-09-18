#! /bin/bash

# hg38 coordinates: chrX:154,100,586-154,342,424 ~240kb

#chm13 coordinates: chrX:152,371,986-152,573,691 ~200kb

while getopts "r:a:s:d:o:h" opt; do
    case $opt in
        r) reference_genome="$OPTARG" ;;
        a) assembly_directory="$OPTARG" ;;
        s) sample_file="$OPTARG" ;;
        d) descriptor="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        h) echo "Usage: $0 -r <reference_genome> -a <assembly_base_directory> -s <sample_file> -d <descriptor> -o <output_dir>"
           echo "Example: $0 -r hg38 -a /path/to/assemblies -s /path/to/samples.txt -d 'sample_name' -o /path/to/output"
           exit 0 ;;
    esac
done

if [ -z "$reference_genome" ] || [ -z "$assembly_directory" ] || [ -z "$sample_file" ] || [ -z "$output_dir" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 -r <reference_genome> -a <assembly_directory> -s <sample_file> -d <descriptor> -o <output_dir>"
    exit 1
fi
mkdir -p $output_dir

while read -r sample sex correct_status; do
    if [ "$sex" == "XY" ]; then
        echo "Processing $sample"
        query_fasta=($( find $assembly_directory -name "$sample*asm.bp.p_ctg.fa" ))
        filename=$(basename "$query_fasta" | cut -d '.' -f1)
        output_paf=$output_dir/${filename}-${sex}-${correct_status}-SVbyEye-paf-alignments/${filename}-${sex}-${descriptor}-${correct_status}.paf
        #Generation of PAF alignments
        mkdir -p $output_dir/${filename}-${sex}-${correct_status}-SVbyEye-paf-alignments
        minimap2 -x asm20 -c -eqx -secondary=no $reference_genome $query_fasta > $output_paf
    fi
    if [ "$sex" == "XX" ]; then
        echo "Processing $sample"
        hap1_fasta=($( find $assembly_directory -name "$sample*asm.hap1.p_ctg.fa" ))
        hap2_fasta=($( find $assembly_directory -name "$sample*asm.hap2.p_ctg.fa" ))
        filename=$(basename "$hap1_fasta" | cut -d '.' -f1)
        output_paf_hap1=$output_dir/${filename}-${sex}-${correct_status}-SVbyEye-paf-alignments/${filename}-${sex}-${descriptor}-${correct_status}.hap1.paf
        output_paf_hap2=$output_dir/${filename}-${sex}-${correct_status}-SVbyEye-paf-alignments/${filename}-${sex}-${descriptor}-${correct_status}.hap2.paf
        #Generation of PAF alignments
        mkdir -p $output_dir/${filename}-${sex}-${correct_status}-SVbyEye-paf-alignments
        minimap2 -x asm20 -c -eqx -secondary=no $reference_genome $hap1_fasta > $output_paf_hap1 
        minimap2 -x asm20 -c -eqx -secondary=no $reference_genome $hap2_fasta > $output_paf_hap2
    fi
done < $sample_file
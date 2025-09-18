#!/bin/bash

# Run QUAST on the assemblies

# Check if the assemblies directory exists

SAMPLE_LIST=$1
THREADS=${2:-32}
ASSEMBLY_BASE_DIR=${3:-"/path/to/assemblies"}
OUTPUT_BASE_DIR=${4:-"./quast-results"}

REFERENCE_ASSEMBLY="/n/zanderson/NSC0202_Kennedy_Macaque_assemblies/rheMac10.fa.gz"
# Create the output directory if it doesn't exist
mkdir -p ${OUTPUT_BASE_DIR}

# Run QUAST on each assembly
declare -a samples
mapfile -t samples < "$SAMPLE_LIST"

for sample in "${samples[@]}"; do
    mkdir -p ${OUTPUT_BASE_DIR}/${sample}-quast-results
    primary_assembly="${ASSEMBLY_BASE_DIR}/${sample}/${sample}.p_ctg.fa"
    if [ -f "$primary_assembly" ]; then
        echo "Running QUAST on $sample primary assembly"
        mkdir -p ${OUTPUT_BASE_DIR}/${sample}-quast-results/primary-assembly
        quast --large --x-for-Nx 75 --fragmented --threads ${THREADS} -r $REFERENCE_ASSEMBLY -o ${OUTPUT_BASE_DIR}/${sample}-quast-results/primary-assembly $primary_assembly
    else
        echo "Primary assembly not found for $sample at $primary_assembly"
    fi
    
    hap1_assembly="${ASSEMBLY_BASE_DIR}/${sample}/${sample}.hap1.fa"
    hap2_assembly="${ASSEMBLY_BASE_DIR}/${sample}/${sample}.hap2.fa"
    if [ -f "$hap1_assembly" ] && [ -f "$hap2_assembly" ]; then
        echo "Running QUAST on $sample haplotype assemblies"
        mkdir -p ${OUTPUT_BASE_DIR}/${sample}-quast-results/hap1-assembly
        quast --large --x-for-Nx 75 --fragmented --threads ${THREADS} -r $REFERENCE_ASSEMBLY -o ${OUTPUT_BASE_DIR}/${sample}-quast-results/hap1-assembly $hap1_assembly 
        mkdir -p ${OUTPUT_BASE_DIR}/${sample}-quast-results/hap2-assembly
        quast --large --x-for-Nx 75 --fragmented --threads ${THREADS} -r $REFERENCE_ASSEMBLY -o ${OUTPUT_BASE_DIR}/${sample}-quast-results/hap2-assembly $hap2_assembly 
    else
        echo "Haplotype 1 or 2 assembly not found for $sample"
        echo "  Looking for: $hap1_assembly"
        echo "  Looking for: $hap2_assembly"
    fi
    echo "Combining transposed results for $sample"
    cat ${OUTPUT_BASE_DIR}/${sample}-quast-results/hap1-assembly/transposed_report.tsv > ${OUTPUT_BASE_DIR}/${sample}-quast-results/combined_transposed.tsv
    awk 'NR==2' ${OUTPUT_BASE_DIR}/${sample}-quast-results/hap2-assembly/transposed_report.tsv >> ${OUTPUT_BASE_DIR}/${sample}-quast-results/combined_transposed.tsv
    awk 'NR==2' ${OUTPUT_BASE_DIR}/${sample}-quast-results/primary-assembly/transposed_report.tsv >> ${OUTPUT_BASE_DIR}/${sample}-quast-results/combined_transposed.tsv
done


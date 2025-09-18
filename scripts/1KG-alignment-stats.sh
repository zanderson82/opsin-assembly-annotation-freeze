#!/bin/bash
eval "$(conda shell.bash hook)"
bam_dir=$1
output_file=$2

#chm13 region: chrX:152382510-152571182
#hg38 region: chrX:154111236-154302341

OPSIN_REGION="chrX:152382510-152571182"
if [ ! -f $output_file ]; then
    echo -e "sample\tsex\tmean_depth\tcoverage\tN50\tmean_length" > $output_file
fi 

while read -r sample sex; do
    # Skip comment lines and empty lines
    [[ "$sample" =~ ^#.* ]] && continue
    [ -z "$sample" ] && continue
    
    bam_file=$(find $bam_dir -name "${sample}*_PMDV_FINAL.haplotagged.bam" | head -n 1)
    if [ -f "$bam_file" ]; then
        echo "Calculating alignment stats for $sample"
        filtered_bam="/tmp/${sample}_filtered.bam"
        samtools view -b -h -F 0x900 "$bam_file" "$OPSIN_REGION" > "$filtered_bam"
        samtools index "$filtered_bam"
        conda activate samtools-1.22
        MEANDEPTH=$(samtools coverage -r $OPSIN_REGION "$filtered_bam" | awk 'NR>1 {print $7}')
        COVERAGE=$(samtools coverage -r $OPSIN_REGION "$filtered_bam" | awk 'NR>1 {print $6}')
        conda deactivate

        # Create filtered BAM for cramino
        #filtered_bam="/tmp/${sample}_filtered.bam"
        
        conda activate cramino-0.14.1
        N50=$(cramino -t 10 "$filtered_bam" | grep "N50" | awk '{print $2}')
        MEANLENGTH=$(cramino -t 10 "$filtered_bam" | grep "Mean length" | awk '{print $3}')
        conda deactivate
        
        # Clean up temporary file
        rm -f "$filtered_bam"
        
        echo -e "$sample\t$sex\t$MEANDEPTH\t$COVERAGE\t$N50\t$MEANLENGTH" >> $output_file
    else
        echo "BAM file not found for $sample"
    fi
done < samples.txt

#! /bin/bash

bam_dir=$1

XIST="chrX:73820656-73852714"
GAPDH="chr12:6534517-6538371"

chrX="chrX:1-156040895"
chr8="chr8:1-145138636"
echo -e "sample\tsex\tmean_depth_XIST\tmean_depth_GAPDH\tchr8\tchrx\tGAPDH/XIST\tchr8/chrx" > /n/zanderson/OPSIN-carrier-screen/manual-assembly-and-annotate/alignment-stats/1KG_X_VS_autosome_hg38.txt

mkdir -p /n/zanderson/OPSIN-carrier-screen/manual-assembly-and-annotate/temp_bams
while read -r sample sex; do 
    bam_file=$(find "$bam_dir" -name "${sample}*.phased.bam" -print -quit)
    if [ -f "$bam_file" ]; then
        XIST_primary_read_bam="/n/zanderson/OPSIN-carrier-screen/manual-assembly-and-annotate/temp_bams/${sample}_XIST.bam"
        GAPDH_primary_read_bam="/n/zanderson/OPSIN-carrier-screen/manual-assembly-and-annotate/temp_bams/${sample}_GAPDH.bam"
        chr8_primary_read_bam="/n/zanderson/OPSIN-carrier-screen/manual-assembly-and-annotate/temp_bams/${sample}_chr8.bam"
        chrX_primary_read_bam="/n/zanderson/OPSIN-carrier-screen/manual-assembly-and-annotate/temp_bams/${sample}_chrX.bam"

        samtools view -b -h -F 0x900 "$bam_file" "$XIST" > "$XIST_primary_read_bam"
        samtools view -b -h -F 0x900 "$bam_file" "$GAPDH" > "$GAPDH_primary_read_bam"
        samtools view -b -h -F 0x900 "$bam_file" "$chr8" > "$chr8_primary_read_bam"
        samtools view -b -h -F 0x900 "$bam_file" "$chrX" > "$chrX_primary_read_bam"

        samtools index "$XIST_primary_read_bam"
        samtools index "$GAPDH_primary_read_bam"
        samtools index "$chr8_primary_read_bam"
        samtools index "$chrX_primary_read_bam"

        echo "Calculating coverage information for X chromosome"
        MEANDEPTH_XIST=$(samtools coverage -r "$XIST" "$XIST_primary_read_bam" | awk 'NR>1 {print $7}')
        echo "Calculating coverage information for GAPDH"
        MEANDEPTH_GAPDH=$(samtools coverage -r "$GAPDH" "$GAPDH_primary_read_bam" | awk 'NR>1 {print $7}')
        echo "Calculating coverage information for autosome"
        MEANDEPTH_CHR8=$(samtools coverage -r "$chr8" "$chr8_primary_read_bam" | awk 'NR>1 {print $7}')
        echo "Calculating coverage information for X chromosome"
        MEANDEPTH_CHRX=$(samtools coverage -r "$chrX" "$chrX_primary_read_bam" | awk 'NR>1 {print $7}')
        echo "Done"
        RATIO=$(echo "$MEANDEPTH_GAPDH / $MEANDEPTH_XIST" | bc -l)
        RATIO_CHR8=$(echo "$MEANDEPTH_CHR8 / $MEANDEPTH_CHRX" | bc -l)

        echo -e "$sample\t$sex\t$MEANDEPTH_XIST\t$MEANDEPTH_GAPDH\t$MEANDEPTH_CHR8\t$MEANDEPTH_CHRX\t$RATIO\t$RATIO_CHR8" >> \
            "/n/zanderson/OPSIN-carrier-screen/manual-assembly-and-annotate/alignment-stats/1KG_X_VS_autosome_hg38.txt"
        rm "$XIST_primary_read_bam" "$GAPDH_primary_read_bam" "$chr8_primary_read_bam" "$chrX_primary_read_bam"
        rm "$XIST_primary_read_bam.bai" "$GAPDH_primary_read_bam.bai" "$chr8_primary_read_bam.bai" "$chrX_primary_read_bam.bai"
    else
        echo "BAM file not found for sample: $sample"
    fi
done < /n/zanderson/OPSIN-carrier-screen/manual-assembly-and-annotate/sample_files/samples.txt

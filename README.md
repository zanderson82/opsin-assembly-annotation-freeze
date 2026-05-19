# Opsin assembly and annotation
Targeted assembly and annotation of the opsin locus at Xq28

## Dependencies
Assembly
* Hifiasm https://github.com/chhylp123/hifiasm
* Samtools 1.23.1

Annotation
* Exonerate https://github.com/nathanweeks/exonerate
* Minimap2 https://github.com/lh3/minimap2
* Samtools 1.23.1

## Usage and Scripts
All scripts needed are in ```./scripts```

### Assembly
For the targeted assembly step, use  ```./scripts/assemble_only.sh```
The input file types can be an aligned bam file or a fastq file. For a fastq, it must be specific to the region of interest.
The inputs are:
* ```samples_file``` : This is a tab separated file with the sample identifier and the sex of the sample
* ```region``` : These are the coordinates that you want to be extracted and assembled
* ```region_name``` : This is a descriptor of the region for naming purposes
* ```bam_dir``` : Directory containing aligned bam files

Options:
* --size [true|false]    Whether to use the --hg-size flag with hifiasm"
                          Default: false if not specified"
                          Use just '--size' (no value) to enable"
* --output-dir <path>    Directory to store all sample output directories"
                          Default: current directory"
* --samtools-flag <flag> Flag to pass to samtools view"
                          Default: '-F 0' (primary alignments)"
                          Example: '--samtools-flag \"-F 0x900\"'"
* --input-suffix <suffix> Default: .phased.bam"
                          Example: '--input-suffix .fastq'"

Example useage: 
```
bash scripts/assemble_only.sh \
"$samples_file" \
"$region" \
"$region_name" \
"$bam_dir" \
--size \
--output-dir $output-dir \
--samtools-flag '-F 0x900' \
--input-suffix .bam 
```

### Annotation
For the annotation and reporting step, use ```./scripts/annotate-count-output-V6-updated.sh```
Inputs:
* ```METADATA_TABLE``` : Tab-separated file with sample ID and sex
* ```ASSEMBLIES_DIR``` : Directory containing the assemblies
* ```OPSIN_REFERENCE``` : Opsin gene sequences (these are provided in ```./gene_sequences```)
* ```OUTPUT_DIR``` : Directory for all annotation outputs
* ```DESCRIPTION```: Descriptor for naming purposes
* ```BAM_DIR``` : Output directory for bam files generated for assembly alignment statistics. These are made by aligning the the fastq file to the assemblies

Example usage: 
```
bash scripts/annotate-count-output-V6-updated.sh \
$METADATA_TABLE \
$ASSEMBLIES_DIR \
$OPSIN_REFERENCE \
$OUTPUT_DIR \
$DESCRIPTION \
$BAM_DIR 
```

import pandas as pd
import numpy as np
import argparse
import re
import os

def gff_to_bed(gff_file, sample_name, sex, output_dir):
    try:
        # Check if input file exists
        if not os.path.exists(gff_file):
            raise FileNotFoundError(f"GFF file not found: {gff_file}")
        
        # Read GFF file with proper column names
        df = pd.read_csv(gff_file, sep='\t', header=None, 
                        names=['seqname', 'source', 'similarity', 'start', 'end', 'score', 'strand', 'blank', 'attribute'])
        
        # Remove empty rows
        df = df.dropna(subset=['seqname'])
        
        if df.empty:
            raise ValueError(f"No valid data found in GFF file: {gff_file}")
        # Convert data types
        df['start'] = pd.to_numeric(df['start'], errors='coerce')
        df['end'] = pd.to_numeric(df['end'], errors='coerce')
        df['score'] = pd.to_numeric(df['score'], errors='coerce')
        df['strand'] = df['strand'].astype(str)
        df['blank'] = df['blank'].astype(str)
        
        # Extract query name from attribute column
        df['query_name'] = df['attribute'].str.extract(r'query=([^;]+)')
        
        # Check if we successfully extracted query names
        if df['query_name'].isna().all():
            print(f"Warning: No query names found in attribute column for {gff_file}")
            df['query_name'] = 'unknown'

        # Handle strand-specific coordinate adjustments (if needed)
        # Note: BED format uses 0-based coordinates, GFF uses 1-based
        # Convert GFF 1-based to BED 0-based by subtracting 1 from start
        df['bed_start'] = df['start'] - 1
        df['bed_end'] = df['end']
    
        # Create BED format: chrom, start, end, name, score, strand
        output_bed = df[['seqname', 'bed_start', 'bed_end', 'query_name', 'score', 'strand']].copy()
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Write BED file
        output_file = f"{output_dir}/{sample_name}-{sex}-annotation-file.bed"
        output_bed.to_csv(output_file, sep='\t', header=False, index=False)
        
        print(f"Successfully converted {len(output_bed)} records from {gff_file} to {output_file}")
        
    except Exception as e:
        print(f"ERROR in gff_to_bed: {str(e)}")
        print(f"Failed to process file: {gff_file}")
        raise




def main():
    parser = argparse.ArgumentParser(description="Convert GFF to BED")
    parser.add_argument("--gff_file", type = str, required = True, help = "Path to the GFF file")
    parser.add_argument("--sample_name", type = str, required = True, help = "Name of the sample")
    parser.add_argument("--sex", type = str, required = True, help = "Sex of the sample")
    parser.add_argument("--output_dir", type = str, required = True, help = "Path to the output directory")
    args = parser.parse_args()

    gff_to_bed(args.gff_file, args.sample_name, args.sex, args.output_dir)

if __name__ == "__main__":
    main()







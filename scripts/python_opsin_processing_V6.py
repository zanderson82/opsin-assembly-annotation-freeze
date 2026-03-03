#!/usr/bin/env python3
"""
Opsin Array Structure Analysis V5

Analyzes GFF annotation files to determine opsin gene array structure.
Handles multiple arrays per contig, fragmented assemblies, and orientation detection.

Coordinate handling:
- For reverse strand features, start > end in the GFF (biological 5' > 3')
- We preserve raw coordinates for biological ordering
- We use normalized coordinates (min/max) for distance calculations

Output:
- Default: TSV format for easy aggregation
- --legacy: Human-readable text format (backward compatible)

Usage:
    python script.py <gff_file> <sample_id> <hap_name> <sex> [options]
"""

import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class Annotation:
    """
    Single annotation from GFF file.
    
    Coordinate conventions:
    - raw_start, raw_end: Original coordinates from GFF (start > end for - strand)
    - norm_start, norm_end: Normalized (norm_start <= norm_end always)
    - For biological ordering, use raw_start
    """
    contig: str
    raw_start: int      # Original start from GFF (5' end, higher for - strand)
    raw_end: int        # Original end from GFF (3' end, lower for - strand)
    norm_start: int     # min(raw_start, raw_end) - for distance calculations
    norm_end: int       # max(raw_start, raw_end) - for distance calculations
    type: str           # 'LCR', 'OPN1LW_exon5', 'OPN1MW_exon5'
    strand: str         # '+' or '-'
    reads: int          # Total reads supporting this annotation
    mapq0: int          # Reads with mapping quality 0
    ratio: float        # MQ0/reads ratio (lower is better)
    
    @property
    def midpoint(self) -> float:
        """Midpoint for distance calculations."""
        return (self.norm_start + self.norm_end) / 2
    
    @property
    def biological_position(self) -> int:
        """Position for biological ordering (5' end of feature)."""
        return self.raw_start


@dataclass 
class Array:
    """
    An opsin array consisting of one LCR and its associated exon5 genes.
    
    Each array determines its own orientation from its genes' strand distribution.
    """
    lcr: Annotation
    exon5s: list[Annotation] = field(default_factory=list)
    contig: str = ""
    _is_forward: Optional[bool] = field(default=None, repr=False)
    
    def __post_init__(self):
        if self.exon5s and self._is_forward is None:
            self._calculate_orientation()
    
    def _calculate_orientation(self):
        """Determine orientation from genes' strand distribution."""
        plus = sum(1 for e in self.exon5s if e.strand == '+')
        minus = sum(1 for e in self.exon5s if e.strand == '-')
        
        if plus > minus:
            self._is_forward = True
        elif minus > plus:
            self._is_forward = False
        else:
            # Ambiguous - use LCR strand as tiebreaker
            self._is_forward = (self.lcr.strand == '+')
    
    @property
    def is_forward(self) -> bool:
        if self._is_forward is None:
            self._calculate_orientation()
        return self._is_forward
    
    @property
    def is_reverse(self) -> bool:
        return not self.is_forward
    
    @property
    def plus_count(self) -> int:
        return sum(1 for e in self.exon5s if e.strand == '+')
    
    @property
    def minus_count(self) -> int:
        return sum(1 for e in self.exon5s if e.strand == '-')
    
    @property
    def gene_count(self) -> int:
        return len(self.exon5s)
    
    @property
    def lw_count(self) -> int:
        return sum(1 for e in self.exon5s if 'LW' in e.type)
    
    @property
    def mw_count(self) -> int:
        return sum(1 for e in self.exon5s if 'MW' in e.type)
    
    @property
    def orientation_is_ambiguous(self) -> bool:
        return self.plus_count == self.minus_count and self.gene_count > 0
    
    @property
    def orientation_clarity(self) -> int:
        return abs(self.plus_count - self.minus_count)
    
    def alignment_score(self) -> tuple:
        """
        Score based on alignment quality metrics.
        Returns tuple where HIGHER = BETTER (for sorting with reverse=True).
        """
        return (-self.lcr.ratio, self.lcr.reads, -self.lcr.mapq0)
    
    def full_score(self) -> tuple:
        """Complete scoring: alignment quality > gene count > orientation clarity."""
        align_score = self.alignment_score()
        return (
            align_score[0],          # -ratio (lower ratio = higher score)
            align_score[1],          # reads (higher = better)
            align_score[2],          # -mapq0 (lower mapq0 = higher score)
            self.gene_count,         # more genes = better
            self.orientation_clarity # clearer orientation = better
        )
    
    def get_sorted_genes(self) -> list[Annotation]:
        """
        Get genes sorted in biological order from LCR outward.
        
        Forward array: genes ordered by ascending position (LCR at low coord)
        Reverse array: genes ordered by descending position (LCR at high coord)
        
        Uses biological_position (raw_start = 5' end) for ordering.
        """
        return sorted(
            self.exon5s, 
            key=lambda x: x.biological_position, 
            reverse=self.is_reverse
        )
    
    def build_structure(self) -> str:
        """Build structure string (without orientation suffix)."""
        sorted_genes = self.get_sorted_genes()
        if not sorted_genes:
            return ""
        
        symbols = []
        for gene in sorted_genes:
            if 'LW' in gene.type:
                symbols.append('L')
            elif 'MW' in gene.type:
                symbols.append('M')
        
        return '-'.join(symbols)
    
    def __str__(self) -> str:
        orientation = "+" if self.is_forward else "-"
        ambig = " (ambiguous)" if self.orientation_is_ambiguous else ""
        return (f"Array on {self.contig}: LCR@{self.lcr.norm_start}-{self.lcr.norm_end} ({self.lcr.strand}), "
                f"{self.gene_count} genes ({self.lw_count}L/{self.mw_count}M), "
                f"orientation={orientation}{ambig}, "
                f"ratio={self.lcr.ratio:.4f}, reads={self.lcr.reads}, MQ0={self.lcr.mapq0}")


def parse_gff(gff_file: str, verbose: bool = False) -> list[Annotation]:
    """
    Parse custom GFF format into Annotation objects.
    
    Expected format (tab-separated, with header):
    contig  start  end  type  strand  reads  MQ0  ratio
    
    For reverse strand features, start > end (biological 5' > 3').
    """
    annotations = []
    
    with open(gff_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            # Skip header, comments, and empty lines
            if line.startswith('contig') or line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            
            if len(fields) < 8:
                if verbose:
                    print(f"  Skipping line {line_num}: insufficient fields ({len(fields)})")
                continue
            
            ann_type = fields[3]
            if ann_type not in ['LCR', 'OPN1LW_exon5', 'OPN1MW_exon5']:
                continue
            
            # Parse coordinates (keep raw values)
            try:
                raw_start = int(fields[1])
                raw_end = int(fields[2])
            except ValueError:
                if verbose:
                    print(f"  Skipping line {line_num}: invalid coordinates")
                continue
            
            # Normalized coordinates for distance calculations
            norm_start = min(raw_start, raw_end)
            norm_end = max(raw_start, raw_end)
            
            # Parse strand
            strand = fields[4]
            if strand not in ['+', '-']:
                # Infer from coordinates if not explicit
                strand = '-' if raw_start > raw_end else '+'
            
            # Parse numeric fields with safe defaults
            try:
                reads = int(fields[5])
            except ValueError:
                reads = 0
            
            try:
                mapq0 = int(fields[6])
            except ValueError:
                mapq0 = 0
            
            try:
                ratio = float(fields[7])
            except ValueError:
                ratio = float('inf')
            
            ann = Annotation(
                contig=fields[0],
                raw_start=raw_start,
                raw_end=raw_end,
                norm_start=norm_start,
                norm_end=norm_end,
                type=ann_type,
                strand=strand,
                reads=reads,
                mapq0=mapq0,
                ratio=ratio
            )
            annotations.append(ann)
            
            if verbose:
                print(f"  Parsed: {ann.type} on {ann.contig} @ {ann.norm_start}-{ann.norm_end} ({ann.strand})")
    
    return annotations


def identify_arrays_on_contig(contig_name: str, annotations: list[Annotation], 
                              verbose: bool = False) -> list[Array]:
    """
    Identify all arrays on a single contig.
    Each gene is assigned to its nearest LCR based on midpoint distance.
    """
    lcrs = [a for a in annotations if a.type == 'LCR']
    exon5s = [a for a in annotations if 'exon5' in a.type]
    
    if not lcrs:
        return []
    
    if not exon5s:
        return [Array(lcr=lcr, exon5s=[], contig=contig_name) for lcr in lcrs]
    
    if len(lcrs) == 1:
        return [Array(lcr=lcrs[0], exon5s=exon5s, contig=contig_name)]
    
    # Multiple LCRs - assign each gene to closest LCR by midpoint distance
    if verbose:
        print(f"  Multiple LCRs on {contig_name}: {len(lcrs)}")
    
    arrays = []
    for lcr in lcrs:
        closest_genes = []
        
        for gene in exon5s:
            # Find which LCR this gene is closest to
            min_distance = float('inf')
            closest_lcr = None
            
            for candidate_lcr in lcrs:
                distance = abs(gene.midpoint - candidate_lcr.midpoint)
                if distance < min_distance:
                    min_distance = distance
                    closest_lcr = candidate_lcr
            
            # If this LCR is the closest, claim this gene
            if closest_lcr is not None and closest_lcr.norm_start == lcr.norm_start:
                closest_genes.append(gene)
        
        arr = Array(lcr=lcr, exon5s=closest_genes, contig=contig_name)
        arrays.append(arr)
        
        if verbose:
            print(f"    LCR@{lcr.norm_start} ({lcr.strand}): {len(closest_genes)} genes, ratio={lcr.ratio:.4f}")
    
    return arrays


def build_orphan_structure(orphan_genes: list[Annotation]) -> str:
    """Build '[M-M]' style string for orphan genes (from contigs without LCRs)."""
    if not orphan_genes:
        return ""
    
    # Sort by contig then by normalized start position
    sorted_genes = sorted(orphan_genes, key=lambda x: (x.contig, x.norm_start))
    
    symbols = []
    for gene in sorted_genes:
        if 'LW' in gene.type:
            symbols.append('L')
        elif 'MW' in gene.type:
            symbols.append('M')
    
    if not symbols:
        return ""
    
    return '[' + '-'.join(symbols) + ']'


def determine_array_structure(annotations: list[Annotation], verbose: bool = False) -> dict:
    """
    Determine opsin array structure from annotations.
    
    Returns dict with all analysis results.
    """
    
    # Empty input handling
    if not annotations:
        return {
            'structure': 'Unknown (no annotations)',
            'lw_count': 0,
            'mw_count': 0,
            'lcr_count': 0,
            'total_genes': 0,
            'total_contigs': 0,
            'contigs_with_lcr': 0,
            'contigs_without_lcr': 0,
            'primary_contig': 'NA',
            'primary_lcr_position': 'NA',
            'primary_lcr_ratio': 'NA',
            'primary_lcr_reads': 'NA',
            'primary_lcr_mapq0': 'NA',
            'is_reverse': False,
            'arrays_found': 0,
            'orphan_genes': 0,
            'orientation_ambiguous': False
        }
    
    # Group by contig
    by_contig = defaultdict(list)
    for ann in annotations:
        by_contig[ann.contig].append(ann)
    
    total_contigs = len(by_contig)
    
    if verbose:
        print(f"\nGrouped annotations into {total_contigs} contig(s)")
    
    # Identify arrays on each contig
    all_arrays = []
    orphan_genes = []
    contigs_with_lcr = 0
    contigs_without_lcr = 0
    
    for contig_name, contig_anns in by_contig.items():
        if verbose:
            print(f"\nProcessing contig: {contig_name}")
            lcr_count = sum(1 for a in contig_anns if a.type == 'LCR')
            gene_count = sum(1 for a in contig_anns if 'exon5' in a.type)
            print(f"  LCRs: {lcr_count}, Genes: {gene_count}")
        
        arrays = identify_arrays_on_contig(contig_name, contig_anns, verbose=verbose)
        
        if arrays:
            all_arrays.extend(arrays)
            contigs_with_lcr += 1
            
            if verbose:
                for arr in arrays:
                    print(f"  Array: {arr}")
        else:
            genes = [a for a in contig_anns if 'exon5' in a.type]
            if genes:
                orphan_genes.extend(genes)
                contigs_without_lcr += 1
                if verbose:
                    print(f"  No LCR - {len(genes)} orphan gene(s)")
    
    # Handle no arrays case
    if not all_arrays:
        all_anns = [a for anns in by_contig.values() for a in anns]
        orphan_structure = build_orphan_structure(orphan_genes)
        structure = orphan_structure if orphan_structure else 'Unknown (no LCR found)'
        
        return {
            'structure': structure,
            'lw_count': sum(1 for a in all_anns if a.type == 'OPN1LW_exon5'),
            'mw_count': sum(1 for a in all_anns if a.type == 'OPN1MW_exon5'),
            'lcr_count': 0,
            'total_genes': sum(1 for a in all_anns if 'exon5' in a.type),
            'total_contigs': total_contigs,
            'contigs_with_lcr': 0,
            'contigs_without_lcr': contigs_without_lcr,
            'primary_contig': 'NA',
            'primary_lcr_position': 'NA',
            'primary_lcr_ratio': 'NA',
            'primary_lcr_reads': 'NA',
            'primary_lcr_mapq0': 'NA',
            'is_reverse': False,
            'arrays_found': 0,
            'orphan_genes': len(orphan_genes),
            'orientation_ambiguous': False
        }
    
    # Sort and select best array
    all_arrays.sort(key=lambda a: a.full_score(), reverse=True)
    best_array = all_arrays[0]
    
    if verbose:
        print(f"\n=== Array Ranking (best first) ===")
        for i, arr in enumerate(all_arrays):
            score = arr.full_score()
            marker = " <-- SELECTED" if i == 0 else ""
            print(f"  {i+1}. score={score[:3]} genes={arr.gene_count} {arr.contig} LCR@{arr.lcr.norm_start}{marker}")
    
    # Build structure string
    primary_structure = best_array.build_structure()
    
    if verbose:
        print(f"\nPrimary structure: {primary_structure}")
        print(f"Orientation: {'reverse' if best_array.is_reverse else 'forward'}")
    
    # Add orphans in brackets
    if orphan_genes:
        orphan_structure = build_orphan_structure(orphan_genes)
        if verbose:
            print(f"Orphan structure: {orphan_structure}")
        
        if primary_structure and orphan_structure:
            combined_structure = f"{primary_structure}-{orphan_structure}"
        elif orphan_structure:
            combined_structure = orphan_structure
        else:
            combined_structure = primary_structure
    else:
        combined_structure = primary_structure
    
    # Add orientation suffix
    if best_array.is_reverse:
        final_structure = f"{combined_structure} (reverse complement)"
    else:
        final_structure = combined_structure
    
    if not final_structure:
        final_structure = "Unknown (no genes)"
    
    if verbose:
        print(f"\nFinal structure: {final_structure}")
    
    # Count totals across ALL annotations
    all_anns = [a for anns in by_contig.values() for a in anns]
    
    return {
        'structure': final_structure,
        'lw_count': sum(1 for a in all_anns if a.type == 'OPN1LW_exon5'),
        'mw_count': sum(1 for a in all_anns if a.type == 'OPN1MW_exon5'),
        'lcr_count': sum(1 for a in all_anns if a.type == 'LCR'),
        'total_genes': sum(1 for a in all_anns if 'exon5' in a.type),
        'total_contigs': total_contigs,
        'contigs_with_lcr': contigs_with_lcr,
        'contigs_without_lcr': contigs_without_lcr,
        'primary_contig': best_array.contig,
        'primary_lcr_position': best_array.lcr.norm_start,
        'primary_lcr_ratio': best_array.lcr.ratio,
        'primary_lcr_reads': best_array.lcr.reads,
        'primary_lcr_mapq0': best_array.lcr.mapq0,
        'is_reverse': best_array.is_reverse,
        'arrays_found': len(all_arrays),
        'orphan_genes': len(orphan_genes),
        'orientation_ambiguous': best_array.orientation_is_ambiguous
    }


# =============================================================================
# Output Functions
# =============================================================================

def get_tsv_columns() -> list[str]:
    """Return list of TSV column names in order."""
    return [
        'sample_id', 
        'sex',
        'haplotype', 
        'structure', 
        'lw_count', 
        'mw_count',
        'total_genes', 
        'lcr_count', 
        'total_contigs', 
        'contigs_with_lcr',
        'contigs_without_lcr', 
        'orphan_genes', 
        'arrays_found',
        'is_reverse', 
        'orientation_ambiguous', 
        'primary_contig',
        'primary_lcr_position', 
        'primary_lcr_ratio', 
        'primary_lcr_reads',
        'primary_lcr_mapq0'
    ]


def get_tsv_header() -> str:
    """Return TSV header line."""
    return '\t'.join(get_tsv_columns())


def format_tsv_line(result: dict, sample_id: str, sex: str, hap_name: str) -> str:
    """Format result as TSV line."""
    cols = get_tsv_columns()
    
    # Add sample info to result
    result['sample_id'] = sample_id
    result['sex'] = sex
    result['haplotype'] = hap_name
    
    values = []
    for col in cols:
        val = result.get(col, 'NA')
        if val is None:
            val = 'NA'
        elif isinstance(val, bool):
            val = str(val).lower()  # 'true' or 'false'
        elif isinstance(val, float):
            if val == float('inf'):
                val = 'NA'
            else:
                val = f"{val:.6f}"
        values.append(str(val))
    
    return '\t'.join(values)


def format_na_line(sample_id: str, sex: str, hap_name: str) -> str:
    """Format a line with NA values for failed samples."""
    cols = get_tsv_columns()
    values = []
    for col in cols:
        if col == 'sample_id':
            values.append(sample_id)
        elif col == 'sex':
            values.append(sex)
        elif col == 'haplotype':
            values.append(hap_name)
        else:
            values.append('NA')
    return '\t'.join(values)


def print_legacy_output(result: dict, sample_id: str, sex: str, hap_name: str):
    """
    Print output in legacy text format for backward compatibility.
    
    This format can be parsed by bash scripts using grep/awk.
    """
    print(f"\n=== Sample {sample_id} ({sex}, {hap_name}) OPSIN Gene Analysis ===")
    print(f"  OPN1LW: {result['lw_count']} copies (based on exon 5 count)")
    print(f"  OPN1MW: {result['mw_count']} copies (based on exon 5 count)")
    print(f"  Total OPSIN genes: {result['total_genes']} copies")
    print(f"  Total LCR count: {result['lcr_count']} LCRs")
    
    print(f"\n=== Contig Information ===")
    print(f"  Total contigs: {result['total_contigs']}")
    print(f"  Contigs with LCR: {result['contigs_with_lcr']}")
    print(f"  Contigs without LCR: {result['contigs_without_lcr']}")
    
    print(f"\n=== Array Structure Analysis ===")
    print(f"  Arrays found: {result['arrays_found']}")
    
    if result['primary_contig'] != 'NA':
        print(f"  Primary contig: {result['primary_contig']}")
        print(f"  Primary LCR position: {result['primary_lcr_position']}")
        ratio_str = f"{result['primary_lcr_ratio']:.4f}" if isinstance(result['primary_lcr_ratio'], float) else result['primary_lcr_ratio']
        print(f"  Primary LCR stats: ratio={ratio_str}, reads={result['primary_lcr_reads']}, MQ0={result['primary_lcr_mapq0']}")
        print(f"  Orientation: {'reverse complement' if result['is_reverse'] else 'forward'}")
    
    if result['orphan_genes'] > 0:
        print(f"  Orphan genes: {result['orphan_genes']} (shown in brackets)")
    
    # This line is parsed by bash script
    print(f"\nFinal structure: {result['structure']}")
    
    if result['orientation_ambiguous']:
        print(f"\nWARNING: Orientation was ambiguous (equal + and - strands)")


# =============================================================================
# Main Entry Point
# =============================================================================

def print_usage():
    """Print usage information."""
    print("Opsin Array Structure Analysis V5")
    print("=" * 40)
    print("")
    print("Usage:")
    print("  python script.py <gff_file> <sample_id> <sex> <haplotype> [options]")
    print("")
    print("Arguments:")
    print("  gff_file    Path to GFF annotation file with alignment stats")
    print("  sample_id   Sample identifier (e.g., 'HG002', 'NA12878')")
    print("  sex         Sample sex ('XX' or 'XY')")
    print("  haplotype   Haplotype name ('hap1', 'hap2', or 'primary')")
    print("")
    print("Options:")
    print("  -v, --verbose    Print detailed debug information")
    print("  --legacy         Output human-readable text (default: TSV)")
    print("  --header         Print TSV header line (for creating new summary file)")
    print("  --na-line        Print NA line for failed sample (requires sample_id, sex, haplotype)")
    print("")
    print("Output formats:")
    print("  Default: Single TSV line (for aggregation into summary file)")
    print("  --legacy: Human-readable text with 'Final structure:' line")
    print("")
    print("Examples:")
    print("  # Create summary file with header")
    print("  python script.py --header > summary.tsv")
    print("")
    print("  # Append sample results")
    print("  python script.py sample.gff SAMPLE1 XX hap1 >> summary.tsv")
    print("  python script.py sample.gff SAMPLE2 XY primary >> summary.tsv")
    print("")
    print("  # Output NA line for failed sample")
    print("  python script.py --na-line SAMPLE3 XX hap1 >> summary.tsv")
    print("")
    print("  # Legacy text output")
    print("  python script.py sample.gff SAMPLE1 XX hap1 --legacy")


if __name__ == "__main__":
    # Handle special flags first
    if len(sys.argv) < 2:
        print_usage()
        sys.exit(1)
    
    # Handle --header flag
    if sys.argv[1] == '--header':
        print(get_tsv_header())
        sys.exit(0)
    
    # Handle --na-line flag
    if sys.argv[1] == '--na-line':
        if len(sys.argv) < 5:
            print("Error: --na-line requires sample_id, sex, and haplotype")
            print("Usage: python script.py --na-line <sample_id> <sex> <haplotype>")
            sys.exit(1)
        sample_id = sys.argv[2]
        sex = sys.argv[3]
        hap_name = sys.argv[4]
        print(format_na_line(sample_id, sex, hap_name))
        sys.exit(0)
    
    # Regular analysis mode - need at least 5 args (script, gff, sample, sex, hap)
    if len(sys.argv) < 5:
        print_usage()
        sys.exit(1)
    
    # Parse positional arguments
    gff_file = sys.argv[1]
    sample_id = sys.argv[2]
    sex = sys.argv[3]
    hap_name = sys.argv[4]
    
    # Validate sex
    if sex not in ['XX', 'XY']:
        print(f"Warning: Unexpected sex value '{sex}' (expected 'XX' or 'XY')")
    
    # Parse optional flags
    verbose = '--verbose' in sys.argv or '-v' in sys.argv
    legacy_output = '--legacy' in sys.argv
    
    # Run analysis
    if verbose:
        print(f"Parsing GFF file: {gff_file}")
        print(f"Sample: {sample_id}, Sex: {sex}, Haplotype: {hap_name}")
    
    try:
        annotations = parse_gff(gff_file, verbose=verbose)
    except FileNotFoundError:
        print(f"Error: GFF file not found: {gff_file}", file=sys.stderr)
        if not legacy_output:
            print(format_na_line(sample_id, sex, hap_name))
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing GFF file: {e}", file=sys.stderr)
        if not legacy_output:
            print(format_na_line(sample_id, sex, hap_name))
        sys.exit(1)
    
    if verbose:
        print(f"Parsed {len(annotations)} OPSIN annotations")
    
    result = determine_array_structure(annotations, verbose=verbose)
    
    # Output
    if legacy_output:
        print_legacy_output(result, sample_id, sex, hap_name)
    else:
        print(format_tsv_line(result, sample_id, sex, hap_name))
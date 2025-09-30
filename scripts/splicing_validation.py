#!/usr/bin/env python3
"""
Splice Variant Validation with Inheritance Awareness
Date: 2025
Description: Validates splice-affecting variants using RNA-seq junction analysis
             with stringent criteria and inheritance pattern consideration.

Requirements: Python 3.6+
No external dependencies required (uses only standard library)
"""

import os
import sys
import argparse
from collections import defaultdict

# Analysis parameters
STRICT_PARAMS = {
    'MIN_READS_PRESENT': 20,
    'MAX_READS_ABSENT': 2,
    'MIN_FOLD_CHANGE': 5
}

HIGH_CONF_PARAMS = {
    'MIN_READS_PRESENT': 50,
    'MAX_READS_ABSENT': 0,
    'MIN_FOLD_CHANGE': 10
}

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Validate splice variants using RNA-seq junctions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Expected file formats:
  Variants file (TSV): sample gene chr position effect inheritance
  Junction files (BED): chr start end name count
  
Sample naming convention:
  Proband: FAMILY.001.junc
  Mother:  FAMILY.002.junc
  Father:  FAMILY.003.junc
  
Example:
  python splicing-validation.py --variants variants.tsv --junctions junction_dir --output results
        """
    )
    parser.add_argument('--variants', required=True, 
                       help='Variant file (TSV format)')
    parser.add_argument('--junctions', required=True, 
                       help='Directory containing junction files')
    parser.add_argument('--output', required=True, 
                       help='Output directory')
    parser.add_argument('--window', type=int, default=1000, 
                       help='Window around variant (bp, default: 1000)')
    parser.add_argument('--mother-suffix', default='.002', 
                       help='Mother sample suffix (default: .002)')
    parser.add_argument('--father-suffix', default='.003', 
                       help='Father sample suffix (default: .003)')
    return parser.parse_args()

def load_variants(variant_file):
    """Load variant information from TSV file
    
    Expected columns: sample gene chr position effect inheritance
    Lines starting with # are treated as comments
    """
    variants = []
    
    if not os.path.exists(variant_file):
        print(f"Error: Variant file not found: {variant_file}", file=sys.stderr)
        sys.exit(1)
    
    try:
        with open(variant_file) as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    print(f"Warning: Line {line_num} has fewer than 6 columns, skipping", 
                          file=sys.stderr)
                    continue
                try:
                    variants.append({
                        'sample': parts[0],
                        'gene': parts[1],
                        'chrom': parts[2],
                        'pos': int(parts[3]),
                        'effect': parts[4],
                        'inheritance': parts[5]
                    })
                except ValueError as e:
                    print(f"Warning: Line {line_num} position not an integer, skipping", 
                          file=sys.stderr)
                    continue
    except Exception as e:
        print(f"Error reading variant file: {e}", file=sys.stderr)
        sys.exit(1)
    
    return variants

def load_junction_counts(junction_file):
    """Load junction counts from BED format file
    
    Expected format: chr start end name count
    Returns dict with junction_key -> count
    """
    junctions = {}
    
    if not os.path.exists(junction_file):
        return junctions
    
    try:
        with open(junction_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    count = int(parts[4])
                    junction_key = f"{chrom}:{start}-{end}"
                    junctions[junction_key] = count
    except Exception as e:
        print(f"Warning: Error reading junction file {junction_file}: {e}", 
              file=sys.stderr)
    
    return junctions

def analyze_variant(variant, junction_dir, window, mother_suffix, father_suffix):
    """Analyze a single variant for splicing changes"""
    
    # Extract family ID from sample name
    family = variant['sample'].split('.')[0]
    
    # Construct filenames for trio
    proband_file = os.path.join(junction_dir, f"{variant['sample']}.junc")
    mother_file = os.path.join(junction_dir, f"{family}{mother_suffix}.junc")
    father_file = os.path.join(junction_dir, f"{family}{father_suffix}.junc")
    
    # Check if files exist
    missing_files = []
    if not os.path.exists(proband_file):
        missing_files.append(proband_file)
    if not os.path.exists(mother_file):
        missing_files.append(mother_file)
    if not os.path.exists(father_file):
        missing_files.append(father_file)
    
    if missing_files:
        print(f"  Warning: Missing junction files: {', '.join(missing_files)}", 
              file=sys.stderr)
        return {
            'variant': variant,
            'aberrant_junctions': [],
            'lost_junctions': [],
            'confidence': 'MISSING_DATA',
            'error': f"Missing files: {', '.join(missing_files)}"
        }
    
    # Load junction data for trio
    proband_junc = load_junction_counts(proband_file)
    mother_junc = load_junction_counts(mother_file)
    father_junc = load_junction_counts(father_file)
    
    # Combine junction data
    all_junctions = set(proband_junc.keys()) | set(mother_junc.keys()) | set(father_junc.keys())
    
    # Filter for junctions near variant
    nearby_junctions = []
    for junction in all_junctions:
        try:
            chrom, coords = junction.split(':')
            start, end = map(int, coords.split('-'))
            
            if chrom == variant['chrom']:
                dist_to_start = abs(start - variant['pos'])
                dist_to_end = abs(end - variant['pos'])
                if min(dist_to_start, dist_to_end) <= window:
                    nearby_junctions.append({
                        'junction': junction,
                        'proband': proband_junc.get(junction, 0),
                        'mother': mother_junc.get(junction, 0),
                        'father': father_junc.get(junction, 0),
                        'distance': min(dist_to_start, dist_to_end)
                    })
        except Exception as e:
            print(f"  Warning: Error parsing junction {junction}: {e}", file=sys.stderr)
            continue
    
    # Analyze based on inheritance pattern
    results = analyze_inheritance_pattern(variant, nearby_junctions)
    return results

def analyze_inheritance_pattern(variant, junctions):
    """Apply appropriate analysis based on inheritance pattern"""
    
    results = {
        'variant': variant,
        'aberrant_junctions': [],
        'lost_junctions': [],
        'confidence': 'LOW'
    }
    
    if variant['inheritance'] == 'de_novo':
        analysis = analyze_de_novo(junctions, STRICT_PARAMS)
    elif variant['inheritance'] in ['maternal', 'paternal']:
        carrier = 'mother' if variant['inheritance'] == 'maternal' else 'father'
        non_carrier = 'father' if carrier == 'mother' else 'mother'
        analysis = analyze_inherited(junctions, carrier, non_carrier, STRICT_PARAMS)
    else:
        print(f"  Warning: Unknown inheritance pattern: {variant['inheritance']}", 
              file=sys.stderr)
        return results
    
    results['aberrant_junctions'] = analysis['aberrant_junctions']
    results['lost_junctions'] = analysis['lost_junctions']
    
    # Assign confidence level
    if len(results['aberrant_junctions']) + len(results['lost_junctions']) > 0:
        if any(j['reads'] >= HIGH_CONF_PARAMS['MIN_READS_PRESENT'] 
               for j in results['aberrant_junctions']):
            results['confidence'] = 'HIGH'
        else:
            results['confidence'] = 'MODERATE'
    
    return results

def analyze_de_novo(junctions, params):
    """Analyze de novo variants
    
    Looks for:
    - Aberrant junctions present in proband but absent in both parents
    - Lost junctions present in both parents but absent in proband
    """
    aberrant = []
    lost = []
    
    for j in junctions:
        parent_max = max(j['mother'], j['father'])
        parent_avg = (j['mother'] + j['father']) / 2
        
        # Check for aberrant junction in proband
        if (j['proband'] >= params['MIN_READS_PRESENT'] and 
            parent_max <= params['MAX_READS_ABSENT']):
            if parent_avg > 0:
                fold = j['proband'] / parent_avg
            else:
                fold = float('inf')
            if fold >= params['MIN_FOLD_CHANGE']:
                aberrant.append({
                    'junction': j['junction'],
                    'reads': j['proband'],
                    'parent_avg': parent_avg,
                    'distance': j['distance']
                })
        
        # Check for lost junction in proband
        if (j['proband'] <= params['MAX_READS_ABSENT'] and 
            parent_avg >= params['MIN_READS_PRESENT']):
            if j['proband'] > 0:
                fold = parent_avg / j['proband']
            else:
                fold = float('inf')
            if fold >= params['MIN_FOLD_CHANGE']:
                lost.append({
                    'junction': j['junction'],
                    'reads': j['proband'],
                    'parent_avg': parent_avg,
                    'distance': j['distance']
                })
    
    return {'aberrant_junctions': aberrant, 'lost_junctions': lost}

def analyze_inherited(junctions, carrier, non_carrier, params):
    """Analyze inherited variants
    
    Looks for:
    - Shared aberrant junctions in proband and carrier parent
    - Lost junctions in both proband and carrier parent
    """
    shared_aberrant = []
    lost = []
    
    for j in junctions:
        carrier_count = j[carrier]
        non_carrier_count = j[non_carrier]
        proband_count = j['proband']
        
        # Shared aberrant junction
        if (proband_count >= params['MIN_READS_PRESENT'] and 
            carrier_count >= params['MIN_READS_PRESENT'] and
            non_carrier_count <= params['MAX_READS_ABSENT']):
            shared_aberrant.append({
                'junction': j['junction'],
                'reads': proband_count,
                'carrier_reads': carrier_count,
                'distance': j['distance']
            })
        
        # Lost in both carriers
        if (proband_count <= params['MAX_READS_ABSENT'] and
            carrier_count <= params['MAX_READS_ABSENT'] and
            non_carrier_count >= params['MIN_READS_PRESENT']):
            lost.append({
                'junction': j['junction'],
                'reads': proband_count,
                'non_carrier_reads': non_carrier_count,
                'distance': j['distance']
            })
    
    return {'aberrant_junctions': shared_aberrant, 'lost_junctions': lost}

def main():
    args = parse_arguments()
    
    # Create output directory
    try:
        os.makedirs(args.output, exist_ok=True)
    except Exception as e:
        print(f"Error creating output directory: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check junction directory exists
    if not os.path.isdir(args.junctions):
        print(f"Error: Junction directory not found: {args.junctions}", file=sys.stderr)
        sys.exit(1)
    
    # Load variants
    variants = load_variants(args.variants)
    print(f"Loaded {len(variants)} variants")
    
    if len(variants) == 0:
        print("No variants to process. Exiting.", file=sys.stderr)
        sys.exit(1)
    
    # Process each variant
    all_results = []
    for variant in variants:
        print(f"\nAnalyzing {variant['sample']} - {variant['gene']} ({variant['effect']})")
        results = analyze_variant(variant, args.junctions, args.window, 
                                 args.mother_suffix, args.father_suffix)
        all_results.append(results)
        
        # Print summary
        if 'error' in results:
            print(f"  Error: {results['error']}")
        else:
            n_aberrant = len(results['aberrant_junctions'])
            n_lost = len(results['lost_junctions'])
            print(f"  Found {n_aberrant} aberrant, {n_lost} lost junctions")
            print(f"  Confidence: {results['confidence']}")
    
    # Write results
    output_file = os.path.join(args.output, 'splice_validation_results.txt')
    try:
        with open(output_file, 'w') as f:
            f.write("Sample\tGene\tEffect\tInheritance\tAberrant\tLost\tConfidence\n")
            for r in all_results:
                v = r['variant']
                f.write(f"{v['sample']}\t{v['gene']}\t{v['effect']}\t{v['inheritance']}\t")
                f.write(f"{len(r['aberrant_junctions'])}\t{len(r['lost_junctions'])}\t")
                f.write(f"{r['confidence']}\n")
        
        print(f"\nResults saved to {output_file}")
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

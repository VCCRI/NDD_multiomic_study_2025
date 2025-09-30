#!/bin/bash

# Monoallelic Expression (MAE) Analysis Pipeline
# Analyzes allele-specific expression in trio families using genomic and RNA-seq data
#
# Requirements:
# - bcftools, samtools, GATK, bedtools
# - VCF file with trio genotypes
# - RNA-seq BAM files for children
# - Reference genome with .dict file

set -e  # Exit on any error

# Configuration - EDIT THESE PATHS FOR YOUR DATA
VCF_FILE="${VCF_FILE:-input/trios.vcf.gz}"
RNA_DIR="${RNA_DIR:-input/rna_bams}"
REF_GENOME="${REF_GENOME:-reference/genome.fa}"
FAMILY_MAP="${FAMILY_MAP:-input/family_mapping.txt}"
OUTPUT_DIR="${OUTPUT_DIR:-output}"

# Create output directories
mkdir -p "${OUTPUT_DIR}/informative_snps"
mkdir -p "${OUTPUT_DIR}/allelic_counts"
mkdir -p "${OUTPUT_DIR}/results"
mkdir -p "temp"

echo "Starting MAE Analysis at $(date)"
echo "Input VCF: ${VCF_FILE}"
echo "RNA directory: ${RNA_DIR}"
echo "Reference: ${REF_GENOME}"
echo "Family mapping: ${FAMILY_MAP}"

# Check required files exist
for file in "${VCF_FILE}" "${REF_GENOME}" "${FAMILY_MAP}"; do
    if [ ! -f "${file}" ]; then
        echo "ERROR: Required file not found: ${file}"
        exit 1
    fi
done

# Check for reference dictionary
DICT_FILE="${REF_GENOME%.*}.dict"
if [ ! -f "${DICT_FILE}" ]; then
    echo "Creating reference dictionary..."
    gatk CreateSequenceDictionary -R "${REF_GENOME}"
fi

# Function to filter VCF for biallelic SNPs
filter_biallelic_snps() {
    local input_vcf=$1
    local output_vcf=$2

    echo "Filtering for biallelic SNPs: $(basename ${input_vcf})"

    bcftools view "${input_vcf}" \
        --min-alleles 2 --max-alleles 2 \
        --types snps \
        --exclude 'INFO/AC=0 | INFO/AN=0' \
        -Oz -o "${output_vcf}.tmp"

    bcftools norm "${output_vcf}.tmp" -d both -Oz -o "${output_vcf}"
    rm -f "${output_vcf}.tmp"

    bcftools index "${output_vcf}"
    gatk IndexFeatureFile -I "${output_vcf}"
}

# Step 1: Identify informative SNPs for each family
echo ""
echo "=== STEP 1: IDENTIFYING INFORMATIVE SNPs ==="

while IFS=$'\t' read -r family child_vcf mother_vcf father_vcf child_bam mother_bam father_bam; do
    # Skip header line
    [[ "$family" == "family_id" ]] && continue

    echo "Processing family ${family}: ${child_vcf}, ${mother_vcf}, ${father_vcf}"

    # Skip if already processed
    if [ -f "${OUTPUT_DIR}/informative_snps/informative_${family}.vcf.gz" ]; then
        echo "Informative SNPs already exist for family ${family}"
        continue
    fi

    # Extract trio samples from main VCF
    bcftools view "${VCF_FILE}" -s "${child_vcf},${mother_vcf},${father_vcf}" \
        -Oz -o "temp/trio_${family}_raw.vcf.gz"

    if [ ! -s "temp/trio_${family}_raw.vcf.gz" ]; then
        echo "WARNING: Failed to extract trio for family ${family}"
        continue
    fi

    # Filter for biallelic SNPs
    filter_biallelic_snps "temp/trio_${family}_raw.vcf.gz" "temp/trio_${family}.vcf.gz"

    # Verify we have 3 samples
    sample_count=$(bcftools query -l "temp/trio_${family}.vcf.gz" | wc -l)
    if [ ${sample_count} -ne 3 ]; then
        echo "WARNING: Expected 3 samples for family ${family}, found ${sample_count}"
        continue
    fi

    # Find informative SNPs (child heterozygous, parents homozygous opposite)
    echo "Finding informative SNPs for family ${family}..."
    bcftools view "temp/trio_${family}.vcf.gz" \
        -i 'GT[0]="0/1" && ((GT[1]="0/0" && GT[2]="1/1") || (GT[1]="1/1" && GT[2]="0/0"))' \
        -Oz -o "temp/informative_${family}_raw.vcf.gz"

    # Clean and sort
    bcftools norm "temp/informative_${family}_raw.vcf.gz" \
        -d both -Oz -o "${OUTPUT_DIR}/informative_snps/informative_${family}.vcf.gz"

    bcftools index "${OUTPUT_DIR}/informative_snps/informative_${family}.vcf.gz"
    gatk IndexFeatureFile -I "${OUTPUT_DIR}/informative_snps/informative_${family}.vcf.gz"

    # Count and report
    informative_count=$(bcftools view -H "${OUTPUT_DIR}/informative_snps/informative_${family}.vcf.gz" | wc -l)
    echo "Family ${family}: ${informative_count} informative SNPs found"

    if [ ${informative_count} -gt 0 ]; then
        # Determine parental origin of reference allele
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' \
            "${OUTPUT_DIR}/informative_snps/informative_${family}.vcf.gz" | \
        awk -v OFS='\t' '{
            chrom = $1; pos = $2; ref = $3; alt = $4;
            child_gt = $5; mother_gt = $6; father_gt = $7;

            # Clean genotypes (remove format fields if present)
            if (index(child_gt, ":") > 0) {
                split(child_gt, c_parts, ":"); child_clean = c_parts[1];
            } else { child_clean = child_gt; }

            if (index(mother_gt, ":") > 0) {
                split(mother_gt, m_parts, ":"); mother_clean = m_parts[1];
            } else { mother_clean = mother_gt; }

            if (index(father_gt, ":") > 0) {
                split(father_gt, f_parts, ":"); father_clean = f_parts[1];
            } else { father_clean = father_gt; }

            # Determine parent of reference allele
            parent_of_ref = "unknown";
            if (mother_clean == "0/0" && father_clean == "1/1") {
                parent_of_ref = "mother";
            }
            else if (mother_clean == "1/1" && father_clean == "0/0") {
                parent_of_ref = "father";
            }

            print chrom, pos, ref, alt, parent_of_ref;
        }' > "${OUTPUT_DIR}/informative_snps/parental_origin_${family}.txt"

        # Report parental origin counts
        mother_count=$(grep -c "mother" "${OUTPUT_DIR}/informative_snps/parental_origin_${family}.txt" || echo 0)
        father_count=$(grep -c "father" "${OUTPUT_DIR}/informative_snps/parental_origin_${family}.txt" || echo 0)
        echo "  Parental origins: ${mother_count} from mother, ${father_count} from father"
    fi

    # Cleanup
    rm -f temp/trio_${family}*.vcf.gz*
    rm -f temp/informative_${family}_raw.vcf.gz*

done < "${FAMILY_MAP}"

# Step 2: RNA-seq allelic counting
echo ""
echo "=== STEP 2: RNA-SEQ ALLELIC COUNTING ==="

while IFS=$'\t' read -r family child_vcf mother_vcf father_vcf child_bam mother_bam father_bam; do
    # Skip header line
    [[ "$family" == "family_id" ]] && continue

    echo "Processing RNA-seq for family ${family}"

    # Check files exist
    if [ ! -f "${RNA_DIR}/${child_bam}" ]; then
        echo "WARNING: RNA BAM not found: ${RNA_DIR}/${child_bam}"
        continue
    fi

    if [ ! -f "${OUTPUT_DIR}/informative_snps/informative_${family}.vcf.gz" ]; then
        echo "WARNING: No informative SNPs for family ${family}"
        continue
    fi

    # Skip if already processed
    if [ -f "${OUTPUT_DIR}/allelic_counts/ase_${family}.txt" ]; then
        echo "ASE results already exist for family ${family}"
        continue
    fi

    # Add read groups to BAM (required by GATK)
    echo "Adding read groups to ${child_bam}..."
    gatk AddOrReplaceReadGroups \
        -I "${RNA_DIR}/${child_bam}" \
        -O "temp/${child_bam}" \
        -RGID "${family}" \
        -RGLB "lib_${family}" \
        -RGPL "ILLUMINA" \
        -RGPU "unit_${family}" \
        -RGSM "${family}_child"

    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to add read groups for family ${family}"
        continue
    fi

    # Run ASEReadCounter with relaxed parameters for RNA-seq
    echo "Running allelic counting for family ${family}..."
    gatk ASEReadCounter \
        -R "${REF_GENOME}" \
        -I "temp/${child_bam}" \
        -V "${OUTPUT_DIR}/informative_snps/informative_${family}.vcf.gz" \
        -O "${OUTPUT_DIR}/allelic_counts/ase_counts_${family}.table" \
        --min-depth 1 \
        --min-mapping-quality 10 \
        --min-base-quality 10

    if [ $? -eq 0 ] && [ -f "${OUTPUT_DIR}/allelic_counts/ase_counts_${family}.table" ]; then
        ase_count=$(tail -n +2 "${OUTPUT_DIR}/allelic_counts/ase_counts_${family}.table" | wc -l)
        echo "ASE counting completed for family ${family}: ${ase_count} positions"

        # Merge with parental origin information
        if [ ${ase_count} -gt 0 ] && [ -f "${OUTPUT_DIR}/informative_snps/parental_origin_${family}.txt" ]; then
            echo "Merging ASE data with parental origin..."

            # Create header
            echo -e "contig\tposition\tvariantID\trefAllele\taltAllele\trefCount\taltCount\ttotalCount\tparentOfRef" > \
                "${OUTPUT_DIR}/allelic_counts/ase_${family}.txt"

            # Sort files for joining
            awk 'NR>1 {print $1":"$2, $0}' "${OUTPUT_DIR}/allelic_counts/ase_counts_${family}.table" | \
            sort -k1,1 > "temp/sorted_ase_${family}.txt"

            awk '{print $1":"$2, $0}' "${OUTPUT_DIR}/informative_snps/parental_origin_${family}.txt" | \
            sort -k1,1 > "temp/sorted_origin_${family}.txt"

            # Join and format
            join -1 1 -2 1 "temp/sorted_ase_${family}.txt" "temp/sorted_origin_${family}.txt" | \
            awk 'BEGIN {OFS="\t"} {print $3,$4,$5,$6,$7,$8,$9,$10,$13}' >> \
                "${OUTPUT_DIR}/allelic_counts/ase_${family}.txt"

            merged_count=$(tail -n +2 "${OUTPUT_DIR}/allelic_counts/ase_${family}.txt" | wc -l)
            echo "Successfully merged ${merged_count} positions for family ${family}"
        fi
    else
        echo "WARNING: ASE counting failed for family ${family}"
    fi

    # Cleanup
    rm -f "temp/${child_bam}"
    rm -f temp/sorted_*_${family}.txt

done < "${FAMILY_MAP}"

# Step 3: Generate summary results
echo ""
echo "=== STEP 3: GENERATING SUMMARY ==="

# Create summary statistics
echo -e "family\tinformative_snps\tmother_origin\tfather_origin\tase_positions\tmonoallelic\tbiallelic\tmaternal_bias\tpaternal_bias" > \
    "${OUTPUT_DIR}/results/summary.txt"

while IFS=$'\t' read -r family child_vcf mother_vcf father_vcf child_bam mother_bam father_bam; do
    # Skip header line
    [[ "$family" == "family_id" ]] && continue

    # Initialize counters
    total_informative=0
    mother_origin=0
    father_origin=0
    ase_positions=0
    monoallelic=0
    biallelic=0
    maternal_bias=0
    paternal_bias=0

    # Count informative SNPs
    if [ -f "${OUTPUT_DIR}/informative_snps/informative_${family}.vcf.gz" ]; then
        total_informative=$(bcftools view -H "${OUTPUT_DIR}/informative_snps/informative_${family}.vcf.gz" | wc -l)
    fi

    # Count parental origins
    if [ -f "${OUTPUT_DIR}/informative_snps/parental_origin_${family}.txt" ]; then
        mother_origin=$(grep -c "mother" "${OUTPUT_DIR}/informative_snps/parental_origin_${family}.txt" 2>/dev/null || echo 0)
        father_origin=$(grep -c "father" "${OUTPUT_DIR}/informative_snps/parental_origin_${family}.txt" 2>/dev/null || echo 0)
    fi

    # Analyze ASE patterns
    if [ -f "${OUTPUT_DIR}/allelic_counts/ase_${family}.txt" ]; then
        ase_positions=$(tail -n +2 "${OUTPUT_DIR}/allelic_counts/ase_${family}.txt" | wc -l)

        # Calculate expression patterns (>90% = monoallelic, 10-90% = biallelic)
        mae_stats=$(awk 'NR>1 && $8>0 {
            ref_freq = $6/$8; alt_freq = $7/$8;

            if (ref_freq > 0.9) {
                mono++;
                if ($9=="mother") mat_bias++; else if ($9=="father") pat_bias++;
            }
            else if (alt_freq > 0.9) {
                mono++;
                if ($9=="father") mat_bias++; else if ($9=="mother") pat_bias++;
            }
            else if (ref_freq >= 0.1 && alt_freq >= 0.1) {
                bi++;
            }
        } END {
            print (mono+0), (bi+0), (mat_bias+0), (pat_bias+0);
        }' "${OUTPUT_DIR}/allelic_counts/ase_${family}.txt")

        read monoallelic biallelic maternal_bias paternal_bias <<< "$mae_stats"
    fi

    # Write summary
    echo -e "${family}\t${total_informative}\t${mother_origin}\t${father_origin}\t${ase_positions}\t${monoallelic}\t${biallelic}\t${maternal_bias}\t${paternal_bias}" >> \
        "${OUTPUT_DIR}/results/summary.txt"

done < "${FAMILY_MAP}"

# Copy key results files
cp "${FAMILY_MAP}" "${OUTPUT_DIR}/results/"
cp "${OUTPUT_DIR}"/allelic_counts/ase_*.txt "${OUTPUT_DIR}/results/" 2>/dev/null || true

# Create README
cat > "${OUTPUT_DIR}/results/README.txt" << EOF
MAE Analysis Results
====================

This directory contains results from monoallelic expression (MAE) analysis of trio families.

Files:
- summary.txt: Summary statistics for all families
- ase_*.txt: Allelic expression data for each family
- family_mapping.txt: Input family-to-sample mapping

Analysis:
- Informative SNPs: Child heterozygous, parents homozygous opposite
- Monoallelic: >90% expression from one allele
- Biallelic: 10-90% expression from both alleles
- Parental bias: Which parent's allele is preferentially expressed

Generated: $(date)
EOF

# Cleanup
rm -rf temp

echo ""
echo "=== ANALYSIS COMPLETE ==="
echo "Results saved in: ${OUTPUT_DIR}/results/"
echo "Summary: ${OUTPUT_DIR}/results/summary.txt"
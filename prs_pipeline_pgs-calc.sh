#!/bin/bash
# Polygenic Risk Score Calculation Pipeline
# Description: Calculates PRS using pgsc_calc for specified traits

#PBS -N prs_calculation
#PBS -l walltime=10:00:00
#PBS -l mem=100GB
#PBS -l ncpus=4

# Set error handling
set -euo pipefail

# Configuration
WORK_DIR="${PBS_O_WORKDIR:-$(pwd)}"
OUTPUT_DIR="${WORK_DIR}/prs_results"
TEMP_DIR="${WORK_DIR}/temp"

# Create directories
mkdir -p ${OUTPUT_DIR} ${TEMP_DIR}

# Load required modules
module load nextflow
module load java

echo "Starting PRS calculation pipeline at $(date)"

# Function to run PRS calculation
calculate_prs() {
    local TRAIT_ID=$1
    local TRAIT_NAME=$2
    local INPUT_VCF=$3
    local OUTPUT_PREFIX="${OUTPUT_DIR}/${TRAIT_NAME}"
    
    echo "Calculating PRS for ${TRAIT_NAME} (${TRAIT_ID})"
    
    # Run pgsc_calc
    nextflow run pgscatalog/pgsc_calc \
        -profile conda \
        --input ${INPUT_VCF} \
        --target_build GRCh38 \
        --pgs_id ${TRAIT_ID} \
        --run_ancestry true \
        --output ${OUTPUT_PREFIX} \
        -work-dir ${TEMP_DIR}/work_${TRAIT_NAME} \
        -resume
    
    if [ $? -eq 0 ]; then
        echo "Successfully calculated PRS for ${TRAIT_NAME}"
    else
        echo "Warning: PRS calculation failed for ${TRAIT_NAME}"
    fi
}

# Main analysis
main() {
    # Check input arguments
    if [ $# -lt 1 ]; then
        echo "Usage: $0 <input_vcf> [samplesheet]"
        exit 1
    fi
    
    INPUT_VCF=$1
    SAMPLESHEET=${2:-""}
    
    # Verify input file exists
    if [ ! -f "${INPUT_VCF}" ]; then
        echo "Error: Input VCF not found: ${INPUT_VCF}"
        exit 1
    fi
    
    # Define traits to analyze
    # Format: PGS_ID,TRAIT_NAME
    TRAITS=(
        "PGS003753,ADHD"
        "PGS001021,Anxiety" 
        "PGS000327,Autism"
        "PGS001829,Depression"
    )
    
    # Calculate PRS for each trait
    for TRAIT_PAIR in "${TRAITS[@]}"; do
        IFS=',' read -r PGS_ID TRAIT_NAME <<< "${TRAIT_PAIR}"
        calculate_prs ${PGS_ID} ${TRAIT_NAME} ${INPUT_VCF}
    done
    
    # Combine results
    echo "Combining PRS results..."
    python3 - <<EOF
import pandas as pd
import os
import glob

output_dir = "${OUTPUT_DIR}"
combined_scores = None

# Find all PRS result files
score_files = glob.glob(f"{output_dir}/*/aggregated_scores.txt")

for score_file in score_files:
    trait_name = os.path.basename(os.path.dirname(score_file))
    
    # Read scores
    scores = pd.read_csv(score_file, sep='\t')
    scores['Trait'] = trait_name
    
    # Combine
    if combined_scores is None:
        combined_scores = scores
    else:
        combined_scores = pd.concat([combined_scores, scores], ignore_index=True)

# Save combined results
if combined_scores is not None:
    combined_scores.to_csv(f"{output_dir}/combined_prs_scores.csv", index=False)
    print(f"Combined {len(score_files)} trait scores")
else:
    print("No PRS scores found to combine")
EOF
    
    # Clean up temporary files
    rm -rf ${TEMP_DIR}/work_*
    
    echo "PRS calculation complete at $(date)"
    echo "Results saved to: ${OUTPUT_DIR}"
}

# Run main analysis
main "$@"

#!/bin/bash

# analyze_mutations_and_apobec3.sh
# Usage: ./analyze_mutations_and_apobec3.sh <variants_dir> <stats_dir> <output_dir>

if [ $# -ne 3 ]; then
    echo "Usage: $0 <variants_dir> <stats_dir> <output_dir>"
    echo "Example: $0 /path/to/variants/ivar /path/to/samtools_stats /path/to/output"
    exit 1
fi

VARIANTS_DIR="$1"
STATS_DIR="$2"
OUTPUT_DIR="$3"

mkdir -p "$OUTPUT_DIR"

# Function to process a single sample
process_sample() {
    local variant_file="$1"
    local stats_file="$2"
    local sample_name="$3"
    
    # Get reads mapped
    local reads=$(grep "reads mapped:" "$stats_file" | awk '{print $4}')
    
    # Calculate mutations - only count SNPs (where ALT is a single base)
    local all_muts=$(awk -F'\t' '$4!~/^-/ && length($4)==1 && $4~/^[ATGC]$/ {count++} END{print (count>0?count:0)}' "$variant_file")
    local major_muts=$(awk -F'\t' '$4!~/^-/ && length($4)==1 && $4~/^[ATGC]$/ && $11>0.5 {count++} END{print (count>0?count:0)}' "$variant_file")
    local major_apobec=$(awk -F'\t' '$4!~/^-/ && length($4)==1 && $4~/^[ATGC]$/ && $11>0.5 && ($3=="G" && $4=="A" || $3=="C" && $4=="T") {count++} END{print (count>0?count:0)}' "$variant_file")
    local minor_muts=$(awk -F'\t' '$4!~/^-/ && length($4)==1 && $4~/^[ATGC]$/ && $11<=0.5 {count++} END{print (count>0?count:0)}' "$variant_file")
    local minor_apobec=$(awk -F'\t' '$4!~/^-/ && length($4)==1 && $4~/^[ATGC]$/ && $11<=0.5 && ($3=="G" && $4=="A" || $3=="C" && $4=="T") {count++} END{print (count>0?count:0)}' "$variant_file")
    
    # Calculate percentages
    local major_apobec_perc=$([ "$major_muts" -gt 0 ] && awk -v ma=$major_apobec -v mm=$major_muts 'BEGIN{printf "%.1f", (ma/mm)*100}' || echo "0.0")
    local minor_apobec_perc=$([ "$minor_muts" -gt 0 ] && awk -v ma=$minor_apobec -v mm=$minor_muts 'BEGIN{printf "%.1f", (ma/mm)*100}' || echo "0.0")
    
    # Calculate APM
    local apm=$(awk -v ma=$minor_apobec -v reads=$reads 'BEGIN{printf "%.1f", (ma/reads)*1000000}')
    
    # Print results in tab-separated format
    printf "%s\t%s\t%s\t%s\t%s\t%.1f\t%s\t%s\t%.1f\t%.1f\n" \
        "$sample_name" "$reads" "$all_muts" "$major_muts" "$major_apobec" \
        "$major_apobec_perc" "$minor_muts" "$minor_apobec" "$minor_apobec_perc" "$apm"
}

# Create and write header to summary file
SUMMARY_FILE="$OUTPUT_DIR/mutations_and_apobec3_summary.txt"
{
    echo "Analysis of mutations (SNPs only, excluding insertions and deletions)"
    echo "================================"
    echo
    echo "NOTE: This analysis focuses on single nucleotide polymorphisms (SNPs) only."
    echo "Results do not include insertions, deletions, or other complex variants."
    echo "Major mutations: frequency > 50%"
    echo "Minor mutations: frequency <= 50%"
    echo "APM: APOBEC3 mutations Per Million MPXV reads"
    echo
    echo "Summary Table - Mutations by Sample"
    echo "------------------------------------------"
    echo
    printf "%-40s %12s %10s %12s %20s %12s %20s %12s\n" \
        "Sample" "MPXV Reads" "All SNPs" "Major SNPs" "Major APOBEC3" "Minor SNPs" "Minor APOBEC3" "APM"
    printf "%s\n" "$(printf '%.0s-' {1..140})"
} > "$SUMMARY_FILE"

# Process all samples and collect results
total_reads=0
total_muts=0
total_major=0
total_major_apobec=0
total_minor=0
total_minor_apobec=0

for variant_file in "$VARIANTS_DIR"/*.tsv; do
    base_name=$(basename "$variant_file" .tsv)
    stats_file="$STATS_DIR/${base_name}.sorted.bam.stats"
    
    if [ ! -f "$stats_file" ]; then
        echo "Warning: Stats file not found for $base_name"
        continue
    fi
    
    result=$(process_sample "$variant_file" "$stats_file" "$base_name")
    echo "$result" >> "$SUMMARY_FILE"
    
    # Update totals
    read -r _ reads all_muts major_muts major_apobec _ minor_muts minor_apobec _ apm <<< "$result"
    total_reads=$((total_reads + reads))
    total_muts=$((total_muts + all_muts))
    total_major=$((total_major + major_muts))
    total_major_apobec=$((total_major_apobec + major_apobec))
    total_minor=$((total_minor + minor_muts))
    total_minor_apobec=$((total_minor_apobec + minor_apobec))
done

# Add separator and totals
{
    printf "%s\n" "$(printf '%.0s-' {1..140})"
    total_major_perc=$([ "$total_major" -gt 0 ] && awk -v ma=$total_major_apobec -v mm=$total_major 'BEGIN{printf "%.1f", (ma/mm)*100}' || echo "0.0")
    total_minor_perc=$([ "$total_minor" -gt 0 ] && awk -v ma=$total_minor_apobec -v mm=$total_minor 'BEGIN{printf "%.1f", (ma/mm)*100}' || echo "0.0")
    total_apm=$(awk -v ma=$total_minor_apobec -v reads=$total_reads 'BEGIN{printf "%.1f", (ma/reads)*1000000}')
    printf "%-40s %12s %10s %12s %12s (%4.1f%%) %12s %12s (%4.1f%%) %12.1f\n" \
        "TOTAL" "$total_reads" "$total_muts" "$total_major" "$total_major_apobec" "$total_major_perc" \
        "$total_minor" "$total_minor_apobec" "$total_minor_perc" "$total_apm"
    echo
    echo "Analysis Summary:"
    echo "------------------------"
    echo "Total samples analyzed: $(ls "$VARIANTS_DIR"/*.tsv | wc -l)"
    echo "Total mutations (all samples): $total_muts"
    echo "Total MPXV reads (all samples): $total_reads"
    echo "Total major mutations: $total_major"
    echo "Total major APOBEC3 mutations: $total_major_apobec ($total_major_perc%)"
    echo "Total minor mutations: $total_minor"
    echo "Total minor APOBEC3 mutations: $total_minor_apobec ($total_minor_perc%)"
    echo "Total APM: $total_apm"
} >> "$SUMMARY_FILE"

echo "Analysis complete. Results written to: $SUMMARY_FILE"

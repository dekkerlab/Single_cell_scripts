#!/bin/bash
#BSUB -J merge_dedup_BX_summary_UU
#BSUB -q long
#BSUB -n 10
#BSUB -W 30:00
#BSUB -R "rusage[mem=46000]"
#BSUB -oo merge_dedup_BX_summary_UU.out
#BSUB -eo merge_dedup_BX_summary_UU.err

input="droplethic_merged.pairs.gz"
dedup_tmp="body.dedup.txt"
sorted_output="droplethic_merged_dedup_sorted.pairs.gz"
summary_output="BX_pairwise_summary_UU_ranked.txt"

# === Step 1: Merge all .pairs.gz files ===
echo "[1/4] Removing headers and merging all .pairs.gz files from ./fastq/"
zcat ./fastq/*.pairs.gz | grep -v '^#' | gzip > "$input"

# === Step 2: Deduplicate pairs ===
echo "[2/4] Deduplicating to reduce data size..."
zcat "$input" | awk '{
    key = $2" "$3" "$4" "$5" "$6" "$7" "$11
    if (!seen[key]) {
        seen[key]=1
        print
    }
}' > "$dedup_tmp"

# === Step 3: Sort deduplicated pairs ===
echo "[3/4] Sorting deduplicated pairs by chrom1, pos1, chrom2, pos2..."
sort -k11,11 -k2,2 -k3,3n -k4,4 -k5,5n --parallel=12 --buffer-size=4G "$dedup_tmp" | gzip > "$sorted_output"

# === Step 4: Summarize BX barcodes for UU pairs ===
echo "[4/4] Summarizing UU pair counts by BX barcode..."
zcat "$sorted_output" | awk '
$8=="UU" {
    BX = $11
    pair = $2":"$3"-"$4":"$5
    total[BX]++
    uniq_pair[BX":"pair] = 1
}
END {
    for (key in uniq_pair) {
        split(key, arr, ":")
        BX = arr[1]
        unique_count[BX]++
    }

    for (bx in total) {
        print bx, total[bx], (bx in unique_count ? unique_count[bx] : 0)
    }
}' | sort -k2,2nr > "$summary_output"

# Optional cleanup
# rm "$dedup_tmp"

echo "✅ All steps completed!"
echo "• Merged: $input"
echo "• Deduplicated + Sorted: $sorted_output"
echo "• BX Summary (UU only): $summary_output"


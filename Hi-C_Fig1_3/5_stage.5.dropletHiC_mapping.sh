#!/bin/bash
#BSUB -J droplethic_pipeline_DEa
#BSUB -q long
#BSUB -n 20
#BSUB -W 24:00
#BSUB -R "rusage[mem=20000]"     
#BSUB -R "span[hosts=1]"
#BSUB -oo pipeline_a.out
#BSUB -eo pipeline_a.err

# === Setup ===
source /home/xiangru.huo-umw/conda/ENTER/etc/profile.d/conda.sh
conda activate droplethic-tools

echo "[1/6] Give the files..."
R1=DE-SI_3A_A1_23322GLT4_S1_L004_R1_001.fastq.gz #change the data name
R3=DE-SI_3A_A1_23322GLT4_S1_L004_R3_001.fastq.gz #change the data name
R2=DE-SI_3A_A1_23322GLT4_S1_L004_R2_001.fastq.gz #change the data name

echo "[2/6] Mapping with BWA-MEM and converting to BAM..."
bwa mem -t 20 /home/xiangru.huo-umw/reference/hg38/hg38.fa \
    "$R3" "$R1" | \
    samtools view -@ 8 -b | \
    samtools sort -n -@ 12 -o sorted_a.bam

echo "[3/6] Extracting read name -> barcode from R2..."
zcat "$R2" | \
awk 'NR % 4 == 1 {read=substr($1, 2)} NR % 4 == 2 {print read "\t" substr($0, 1, 16)}' \
> readname_to_bx_a.txt

echo "[4/6] Adding BX tags to BAM using Python..."
python <<EOF
import pysam

# Load barcode map
bx = {}
with open("readname_to_bx_a.txt") as f:
    for line in f:
        name, code = line.strip().split()
        bx[name] = code

inbam = pysam.AlignmentFile("sorted_a.bam", "rb")
outbam = pysam.AlignmentFile("with_bx_a.bam", "wb", header=inbam.header)

for read in inbam:
    tag = bx.get(read.query_name.split()[0], "0")
    read.set_tag("BX", tag, value_type='Z')
    outbam.write(read)

inbam.close()
outbam.close()
EOF

echo "[5/6] Running pairtools parse and sort..."
pairtools parse --add-columns mapq,BX --drop-sam \
	-c /home/xiangru.huo-umw/reference/hg38/hg38.chrom.sizes \
	with_bx_a.bam | \
pairtools sort --nproc 20 -o droplethic_output_a.pairs.gz

echo "[6/6] Cleaning up intermediate files..."
#rm readname_to_bx_a.txt sorted_a.bam with_bx_a.bam

echo "Pipeline completed successfully!"


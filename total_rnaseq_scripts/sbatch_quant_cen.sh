#!/usr/bin/env bash
#SBATCH --job-name=cen_quant
#SBATCH --partition=nice
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

module load bedtools2/2.31.0-xmfvklh

samplename=$1

echo "${samplename}"

bedtools intersect -c -sorted \
    -a /home/ss18818/projects/trs/censat_hor_with_strand_R_numbered_chrsorted.bed \
    -b /home/ss18818/projects/trs/star_out_filtered/${samplename}_99_filtered_sorted.bam \
    -g /home/ss18818/projects/trs/t2tv2_ncbi.chrom.sizes \
> /home/ss18818/projects/trs/cen_quant_bed/${samplename}_99_filtered_sorted.bam.bed

echo "Quantified flag 99 BAM"

bedtools intersect -c -sorted \
    -a /home/ss18818/projects/trs/censat_hor_with_strand_R_numbered_chrsorted.bed \
    -b /home/ss18818/projects/trs/star_out_filtered/${samplename}_83_filtered_sorted.bam \
    -g /home/ss18818/projects/trs/t2tv2_ncbi.chrom.sizes \
> /home/ss18818/projects/trs/cen_quant_bed/${samplename}_83_filtered_sorted.bam.bed

echo "Quantified flag 83 BAM"

bedtools intersect -c -sorted \
    -a /home/ss18818/projects/trs/censat_hor_with_strand_R_numbered_chrsorted.bed \
    -b /home/ss18818/projects/trs/star_out_filtered/${samplename}_99_83_filtered_sorted.bam \
    -g /home/ss18818/projects/trs/t2tv2_ncbi.chrom.sizes \
> /home/ss18818/projects/trs/cen_quant_bed/${samplename}_99_83_filtered_sorted.bam.bed

echo "Quantified the merged BAM"


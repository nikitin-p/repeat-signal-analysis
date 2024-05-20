#!/usr/bin/env bash
#SBATCH --job-name=bam_to_bigwig
#SBATCH --partition=nice
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

module load py-deeptools/3.5.3-2ap3qbb  

samplename=$1

echo "${samplename}"

bamCoverage \
    --binSize 1 \
    --numberOfProcessors 8 \
    --normalizeUsing BPM \
    -b /home/ss18818/projects/trs/star_out_filtered/${samplename}_99_filtered_sorted.bam \
    -o /home/ss18818/projects/trs/total_rnaseq_bigwig/${samplename}_99_filtered_sorted.reverse.bigWig

echo "Generated bigWig for 99."

bamCoverage \
    --binSize 1 \
    --numberOfProcessors 8 \
    --normalizeUsing BPM \
    -b /home/ss18818/projects/trs/star_out_filtered/${samplename}_83_filtered_sorted.bam \
    -o /home/ss18818/projects/trs/total_rnaseq_bigwig/${samplename}_83_filtered_sorted.forward.bigWig

echo "Generated bigWig for 83."


#!/usr/bin/env bash
#SBATCH --job-name=filter_and_merge
#SBATCH --partition=nice
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G

samplename=$1

module load samtools/1.19.2-tp4qbgg

echo "${samplename}"

samtools view \
    -f 99 \
    -@ 8 \
    -b \
    -o star_out_filtered/${samplename}_99_filtered.bam \
    star_out/${samplename}Aligned.sortedByCoord.out.bam

echo "Filtered for flag 99"

samtools sort \
    -@ 8 \
    -O bam \
    -o star_out_filtered/${samplename}_99_filtered_sorted.bam \
    star_out_filtered/${samplename}_99_filtered.bam

echo "Sorted flag 99 BAM"

samtools index \
    -@ 8 \
    -o star_out_filtered/${samplename}_99_filtered_sorted.bai \
    star_out_filtered/${samplename}_99_filtered_sorted.bam

echo "Indexed flag 99 BAM"

samtools view \
    -f 83 \
    -@ 8 \
    -b \
    -o star_out_filtered/${samplename}_83_filtered.bam \
    star_out/${samplename}Aligned.sortedByCoord.out.bam \

echo "Filtered for flag 83"

samtools sort \
    -@ 8 \
    -O bam \
    -o star_out_filtered/${samplename}_83_filtered_sorted.bam \
    star_out_filtered/${samplename}_83_filtered.bam

echo "Sorted flag 83 BAM"

samtools index \
    -@ 8 \
    -o star_out_filtered/${samplename}_83_filtered_sorted.bai \
    star_out_filtered/${samplename}_83_filtered_sorted.bam

echo "Indexed flag 83 BAM"

samtools merge \
    -@ 8 \
    -o star_out_filtered/${samplename}_99_83_filtered.bam \
    star_out_filtered/${samplename}_99_filtered_sorted.bam \
    star_out_filtered/${samplename}_83_filtered_sorted.bam
 
echo "Merged flag 99 and flag 83 BAMs"

samtools sort \
    -@ 8 \
    -O bam \
    -o star_out_filtered/${samplename}_99_83_filtered_sorted.bam \
    star_out_filtered/${samplename}_99_83_filtered.bam \

echo "Sorted the merged BAM"

samtools index \
    -@ 8 \
    -o star_out_filtered/${samplename}_99_83_filtered_sorted.bai \
    star_out_filtered/${samplename}_99_83_filtered_sorted.bam

echo "Indexed the merged BAM"


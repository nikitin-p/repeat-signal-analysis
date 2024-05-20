#!/usr/bin/env bash
#SBATCH --job-name=filter_and_merge_se
#SBATCH --partition=nice
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G

samplename=$1

module load samtools/1.19.2-tp4qbgg

echo "${samplename}"

samtools view \
    -F 4 \
    -@ 8 \
    -b \
    -o star_out_filtered/${samplename}_mapped.bam \
    star_out/${samplename}Aligned.sortedByCoord.out.bam

echo "Filtered for mapped reads"

samtools view \
    -f 16 \
    -@ 8 \
    -b \
    -o star_out_filtered/${samplename}_mapped_f16.bam \
    star_out_filtered/${samplename}_mapped.bam \

echo "Filtered for flag 16"

samtools sort \
    -@ 8 \
    -O bam \
    -o star_out_filtered/${samplename}_mapped_f16_sorted.bam \
    star_out_filtered/${samplename}_mapped_f16.bam 

echo "Sorted flag 16 BAM"

samtools index \
    -@ 8 \
    -o star_out_filtered/${samplename}_mapped_f16_sorted.bai \
    star_out_filtered/${samplename}_mapped_f16_sorted.bam 

echo "Indexed flag 16 BAM"

samtools view \
    -F 16 \
    -@ 8 \
    -b \
    -o star_out_filtered/${samplename}_mapped_F16.bam \
    star_out_filtered/${samplename}_mapped.bam \

echo "Filtered against flag 16"

samtools sort \
    -@ 8 \
    -O bam \
    -o star_out_filtered/${samplename}_mapped_F16_sorted.bam \
    star_out_filtered/${samplename}_mapped_F16.bam 

echo "Sorted flag non-16 BAM"

samtools index \
    -@ 8 \
    -o star_out_filtered/${samplename}_mapped_F16_sorted.bai \
    star_out_filtered/${samplename}_mapped_F16_sorted.bam 

echo "Indexed flag non-16 BAM"

samtools merge \
    -@ 8 \
    -o star_out_filtered/${samplename}_f16_F16_merged.bam \
    star_out_filtered/${samplename}_mapped_f16_sorted.bam \
    star_out_filtered/${samplename}_mapped_F16_sorted.bam 
 
echo "Merged flag 16 and non-flag 16 BAMs"

samtools sort \
    -@ 8 \
    -O bam \
    -o star_out_filtered/${samplename}_f16_F16_merged_sorted.bam \
    star_out_filtered/${samplename}_f16_F16_merged.bam 

echo "Sorted the merged BAM"

samtools index \
    -@ 8 \
    -o star_out_filtered/${samplename}_f16_F16_merged_sorted.bai \
    star_out_filtered/${samplename}_f16_F16_merged_sorted.bam 

echo "Indexed the merged BAM"


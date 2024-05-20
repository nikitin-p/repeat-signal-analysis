#!/usr/bin/env bash
#SBATCH --job-name=airway_iso1_star
#SBATCH --partition=nice
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=36G

module load star/2.7.11a-pgsk3s4

samplename=$1

STAR \
    --runThreadN 8 \
    --genomeDir /home/ss18818/projects/trs/star_index_t2t \
    --readFilesIn /home/ss18818/projects/trs/total_rnaseq/${samplename}_R1.fastq.gz \
    --readFilesCommand gunzip -c \
    --outTmpDir /home/ss18818/mnt/scratch/ss18818/${samplename}_tmp \
    --outSAMtype BAM SortedByCoordinate \
    --winAnchorMultimapNmax 500 \
    --outFilterMultimapNmax 100 \
    --outSAMmultNmax 1 \
    --outFilterMatchNminOverLread 0.9 \
    --outFileNamePrefix /home/ss18818/projects/trs/star_out/${samplename}


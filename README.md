# repeat-signal-analysis
Nextflow pipeline for analysis of short-read data (CAGE and total RNA-seq) with T2T-CHM13v2.0 human reference genome.

Description of files:

* `download_reads.txt` Download all the FASTQ libraries of total RNA-seq and CAGE.

* `analysis.txt` All commands used for mappping reads and further analysis.

* `hor_strand.Rmd` BED of centromere alpha-satellite HOR repeat arrays with strand-specificity.

* `total_rna_centromere_heatmap.Rmd` Heatmaps of total RNA expression counts in centromere HOR arrays.

* `cage_centromere_heatmap.Rmd` Heatmaps of CAGE expression counts in centromere HOR arrays.

* `total_rnaseq_scripts` Scripts for total-rna seq mapping and analysis

    * `download_fastq.sh` Download reads
    
    * `.yaml` and `.csv` files contain parameters for nfcore/rnaseq run

    * `t2tv2_ncbi.chrom.sizes` file with T2T-CHM13v2.0 chromosomes sizes for https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009914755.4/

    * `generate_t2t_chrom_sizes.txt` script for chromosomez sizes caclulation

    * `tool_versions.txt` tools versions

    * `sbatch_quant_cen.sh` and `sbatch_quant_cen_se.sh` sbatch scripts for quantifocation of total RNA-seq

sbatch_quant_cen.sh и sbatch_quant_cen_se.sh
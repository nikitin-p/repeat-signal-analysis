# repeat-signal-analysis
Nextflow pipeline for analysis of short-read data (CAGE and total RNA-seq) with T2T-CHM13v2.0 human reference genome.

Description of files:

* `download_reads.txt` Download all the FASTQ libraries of total RNA-seq and CAGE.

* `analysis.txt` All commands used for mappping reads and further analysis.

* `hor_strand.Rmd` BED of centromere alpha-satellite HOR repeat arrays with strand-specificity.

* `total_rna_centromere_heatmap.Rmd` Heatmaps of total RNA expression counts in centromere HOR arrays.

* `cage_centromere_heatmap.Rmd` Heatmaps of CAGE expression counts in centromere HOR arrays.

# Data analysis of the south Wales Geochemical Gradient (SGG) 16S rRNA gene dataset

This repository contains all R and bash scripts utilised to process and analyse the SGG 16S rRNA gene dataset.
Resulting data is under preparation for submission under the title "Hydrogeological controls of biogeochemically relevant bacterial consortia in subsurface mine waters", aimed at **ISMEJ**.

## Folder organization and content:

 - *00_Data:*

	Metadata; if cloned, can contain raw FASTQ generated for this study.

 - *01_QC_Filtering:*

	`trim_16S_reads.sh` trims and filters Illumina adapter sequences with [`cutadapt`](https://github.com/marcelm/cutadapt)

 - *02_SV_Processing:*
	
	`DADA2_SGG_BDataPip.R` uses [`DADA2`](https://github.com/benjjneb/dada2) to infer SVs (sequence variants) from trimmed, quality-filtered FASTQs.

 - *03_Reports:*

	`SGG-Manuscript-Figures_Data_Pond.Rmd`
	
	`SGG-Manuscript-Figures_Data_Pond_BW.Rmd`

	`sgg_geochem.Rmd`

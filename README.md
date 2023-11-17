# edna_citizen

Here, we use environmental DNA based sequencing to explore fish diversity in the Oslofjord and compare this data with other citizen science (Artsdatabanken) and the beach seine survey. This repository includes scripts to run cutadapt (https://github.com/marcelm/cutadapt), dada2 (https://github.com/benjjneb/dada2) and statistical analyses for obtaining fish diversity data as well as data comparisons and visualization.

## pull repository

```
cd path to repositories
git clone github.com/alper1976/edna_citizen.git
```

## Authors and acknowledgment
Scripts were written by Lone Kvalheim, Eivind Stensrud, Olli Hyvärinen and Alexander Eiler.

## License
This Code is subject to the terms of the MIT License. 

## Project status
Results from this project have been submitted to a peer-reviewed scientific journal.

## Folders and code
The analyses code is divided into multiple folders "data", "seq_data_processing" and "stats" representing the code to analyze raw sequencing data and perform statistical analysis, respectively.

### data
The "data" folder contains the output data from the sequence data analysis. In addition, it contains mined data from Artsdatabanken and data from the beach seine survey.

environmental_data.csv - data on water properties
hal.data2.csv
hal.data.csv
historic.data_big.csv
historic.data.csv
phyloseq_mifish_lca_right.rds
phyloseq_mifish_lca.rds
species_table_n.csv - species tables from the eDNA analysis
species_table.csv - species tables from the eDNA analysis
strandnot 2010-2022.csv - data from the beach seine survey

### seq_data_processing
Here you can find the code to obtain the Amplicon sequence variants (ASV) tables from the raw sequencing data. You can also find code for the taxonomic analyses using an in-house fish database (ScandiFish) and LCA analysis.

### stats
Here you can find the code for running statistical analyses and visualizations.

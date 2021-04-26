# Research Compendium for the VivID Project (Genetic Portion)
### The Epidemiology of Plasmodium vivax Among Adults in the Democratic Republic of the Congo
#### PMCID: 

```
├── PopGenome_Analysis  
├── README.md
├── eda_batch_effects # Exploratory 
├── hardfilters_variants
├── lists
├── raw_variants
├── regions  
├── scrape_pubseqs
├── wgs_pe_improved_ViVIDSmpls
├── wgs_pe_improved_global
├── wgs_qc_improved_ViVIDSmpls
├── wgs_qc_improved_global
└── wgs_se_improved_global

```

The `PopGenome_Analysis` directory contains the various analysis scripts, which follow a conventional numbering system for order in which they are to be run. 

VCFs will need to be created by first downloading the publicly available sequences with the snakemake modules in `scrape_pubseqs`. Sequences will then need to be aligned with `wgs_pe_improved*` snakemake modules with respect to global and VivID samples. QC will need to be performed in the same manner. Then, variants can be called with the snakemake modules in `raw_variants` followed by application of filtering with the `hardfilters_variants` module. 
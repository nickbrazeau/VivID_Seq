# README

# Next-Generation Sequences
## Step 0: Download Publicly Available Sequences
### Paired End Reads
1. From the `vivid_seq_public_NGS.xlsx` file scrape ENA to convert SRS to SRR and then create PE and SE master tab.txt files for download. This is all done with the `01-scrape_ena_NGS.R` script.  
2. Run the scripts within `scrape_pubseqs/snakescrape_ENA_PE`

### Single End Reads
1. From above -- scraped ENA and made SE master tab.txt file.
2. Run the scripts within `scrape_pubseqs/snakescrape_ENA_SE`   

# Sanger Sequences
## Step 0: Download Publicly Available Sequences
1. From the `vivid_seq_public_sangerseq.xlsx` primary taken from [Rodrigues et al. 2018](https://www.nature.com/articles/s41598-018-19554-0) download sanger fastas from accession numbers through NCBI (using `rentrez`). This process is wrapped in the `03-scrape_ncbi_sanger.R` script
2. Resulting downloaded multiple-sequence fasta is written to: `TODO`

## Step 1: MARS and MSA

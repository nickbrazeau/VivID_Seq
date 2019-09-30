
# Download Next-Generation Sequences
## Step 0: Download Publicly Available Sequences (NGS)
### Paired End Reads
1. From the `vivid_seq_public_NGS.xlsx` file scrape ENA to convert SRS to SRR and then create PE and SE master tab.txt files for download. This is all done with the `01-scrape_ena_NGS.R` script.  
2. Run the scripts within `scrape_pubseqs/snakescrape_ENA_PE`

### Single End Reads
1. From above -- scraped ENA and made SE master tab.txt file.
2. Run the scripts within `scrape_pubseqs/snakescrape_ENA_SE`   

## Step 0: Download Publicly Available Sequences (Sanger)
1. From the `vivid_seq_public_sangerseq.xlsx` primary taken from [Rodrigues et al. 2018](https://www.nature.com/articles/s41598-018-19554-0) download sanger fastas from accession numbers through NCBI (using `rentrez`). This process is wrapped in the `03-scrape_ncbi_sanger.R` script
2. Resulting downloaded multiple-sequence fasta is written to: `TODO`

# Polarize mtDNA
## Step 1: MARS and MSA
_TO DO_

# Align Sequences
## Step 1.1: WGS PE Global
_TO DO_ -- `wgs_pe_improved_global`

## Step 1.2: WGS SE Global
_TO DO_ -- `wgs_se_improved_global`

## Step 1.3: WGS SE VivID
_TO DO_ -- `wgs_pe_improved_ViVIDSmpls`

# Quality Control
## Step 2.1: Global QC
_TO DO_ -- `wgs_qc_improved_global`

## Step 2.2: VivID QC
_TO DO_ -- `wgs_qc_improved_ViVIDSmpls`

# Variant Calling and Recalibration
## Step 3: Find _Pv_ Lab Strain Gold Standard Variants
_TO DO_ -- `labstrains_variants`

## Step 4: Variant Filter _Pv_ Calls
_TO DO_-- `filtered_variants`

# Cell hashing dataset (Stoeckius et al)
The raw **cell hashing** datasets can be downloaded from [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1603-1) or GEO database.
The RData file used in our analysis can be downloaded from [google drive](https://drive.google.com/file/d/12y0fW_Y9OdhBLns_2gpjo2Xq25c4qnGY/view?usp=sharing). It contains both the count matrix and label for each barcode.

# 10X PBMC dataset
The raw **10x 10k PBMC** dataset (cellranger 3.0.0) was downloaded from https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3. Detailed information can be found in [10x_PBMC_dataset_page](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3).

# 10x 900 Neuron dataset
The raw 10x **Neuron 900** dataset (cellranger 2.1.0) was downloaded from https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/neurons_900. Detailed information can be found in [10x_Neuron_dataset_page](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/neurons_900).

# Kang et al sampleC dataset
The raw bam file was downloaded from [SRA website](https://sra-pub-src-1.s3.amazonaws.com/SRR5398237/C.merged.bam.1). More information can be found in [SRR5398237](https://www.ncbi.nlm.nih.gov/sra/SRX2693024[accn]). We used "bamtofastq" tool to generate appropriate formatted fastq files and then applied cellranger 2.0 to get the raw count matrix, which can be downloaded from [google drive](https://drive.google.com/file/d/1Z1Vxzpu17kWwZGo6f2BMKo9eLjofmdrk/view?usp=sharing).

#!/bin/bash
# Cell hashing dataset (Stoeckius et al.)
fileId=12y0fW_Y9OdhBLns_2gpjo2Xq25c4qnGY
fileName=Cell_hashing_raw.RData
curl -sc /tmp/cookie "https://drive.google.com/uc?export=download&id=${fileId}" > /dev/null
code="$(awk '/_warning_/ {print $NF}' /tmp/cookie)"
curl -Lb /tmp/cookie "https://drive.google.com/uc?export=download&confirm=${code}&id=${fileId}" -o ./data/${fileName}

# Multiplexed PBMC dataset (Kang et al.)
wget -O ./data/Multiplexed_PBMC_raw.h5 "https://drive.google.com/uc?export=download&id=1Z1Vxzpu17kWwZGo6f2BMKo9eLjofmdrk"

# Two additional 10x datasets
wget -O ./data/PBMC_10k_raw.h5 "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_raw_feature_bc_matrix.h5"
wget -O ./data/Neurons_900_raw.h5 "http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neurons_900/neurons_900_raw_gene_bc_matrices_h5.h5"

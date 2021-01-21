#!/bin/bash
# Cell hashing dataset (Stoeckius et al.)
fileName=Cell_hashing_raw.RData
fileId=12y0fW_Y9OdhBLns_2gpjo2Xq25c4qnGY
curl -sc /tmp/cookie "https://drive.google.com/uc?export=download&id=${fileId}" > /dev/null
code="$(awk '/_warning_/ {print $NF}' /tmp/cookie)"
curl -Lb /tmp/cookie "https://drive.google.com/uc?export=download&confirm=${code}&id=${fileId}" -o ./data/${fileName}
# Multiplexed PBMC dataset (Kang et al.)
wget -O ./data/Multiplexed_PBMC_raw.h5 "https://drive.google.com/uc?export=download&id=1Z1Vxzpu17kWwZGo6f2BMKo9eLjofmdrk"
# Two additional 10x datasets
wget -O ./data/PBMC_10k_raw.h5 "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_raw_feature_bc_matrix.h5"
wget -O ./data/Neurons_900_raw.h5 "http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neurons_900/neurons_900_raw_gene_bc_matrices_h5.h5"

# Other files to reproduce the results in manuscript
fileName=pbmc8k_retained.rds
fileId=1YiYD0yt8dbkK0AC6jm9g2o85tNRsiryZ
curl -sc /tmp/cookie "https://drive.google.com/uc?export=download&id=${fileId}" > /dev/null
code="$(awk '/_warning_/ {print $NF}' /tmp/cookie)"
curl -Lb /tmp/cookie "https://drive.google.com/uc?export=download&confirm=${code}&id=${fileId}" -o ./data/${fileName}
wget -O ./data/neuron900_retained.rds "https://drive.google.com/uc?export=download&id=1cwz4KDW4KdMziYRKVTieyu4mqfSjsYL-"
wget -O ./data/BlueYellowColormaps_V1.RData "https://drive.google.com/uc?export=download&id=1WHv5il6L26hAii961wKfvf3AZFZenpx9"

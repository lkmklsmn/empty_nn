#!/bin/bash
# Raw datasets
## Cell hashing dataset (Stoeckius et al.)
function download_large_googlefiles {
  curl -sc /tmp/cookie "https://drive.google.com/uc?export=download&id=${fileId}" > /dev/null
  code="$(awk '/_warning_/ {print $NF}' /tmp/cookie)"
  curl -Lb /tmp/cookie "https://drive.google.com/uc?export=download&confirm=${code}&id=${fileId}" -o ./data/${fileName}
}
fileName=cell_hashing_raw.RData
fileId=12y0fW_Y9OdhBLns_2gpjo2Xq25c4qnGY
download_large_googlefiles

## Multiplexed PBMC dataset (Kang et al.)
wget -O ./data/multiplexed_PBMC_raw.h5 "https://drive.google.com/uc?export=download&id=1Z1Vxzpu17kWwZGo6f2BMKo9eLjofmdrk"

## Two additional 10x datasets
wget -O ./data/pbmc_8k_raw.h5 "https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_raw_gene_bc_matrices_h5.h5"
wget -O ./data/neurons_900_raw.h5 "http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neurons_900/neurons_900_raw_gene_bc_matrices_h5.h5"

# Other files to reproduce the results in manuscript
## reproduce Figure 2
wget -O ./data/cell_hashing_results.RData "https://drive.google.com/uc?export=download&id=1AXUzQEp1OJmZgLTwTHf1szqD-9nE092v"
fileName=reference_pbmc_3k.rds
fileId=1qsVgHJ6TmorW4lxFbD9x-iiBNKZY-5b4
download_large_googlefiles
fileName=cell_hashing_original_paper.rds
fileId=1oviPPff0-GV6t4p-eaFx5gF8fbU0-m-u
download_large_googlefiles
fileName=cell_hashing_retained.rds
fileId=114emrk75iiF9Eu4QHzk_0-gYljiVBHMa
download_large_googlefiles

## reproduce Figure 3
wget -O ./data/multiplexed_PBMC_results.RData "https://drive.google.com/file/d/1EqpH_0qmXTPcwV29NYPcWhA0uuTIx3xb/view?usp=sharing"
wget -O ./data/multiplexed_PBMC_demuxlet.best "https://drive.google.com/uc?export=download&id=1fVvTZE8GE3424Dp1xCqethe_8qtv5b8n"

## reproduce Figure 4
wget -O ./data/mouse_nuclei_retained.rds "https://drive.google.com/file/d/19UhPPCy_f-wxDD-7n6DVJkpiLWCcIXPD/view?usp=sharing"
wget -O ./data/mouse_nuclei_2k_raw.rds "https://drive.google.com/file/d/1Vh5iB_zHFhDSvGO7KG50zzLZwofYElia/view?usp=sharing"
fileName=unspliced.csv
fileId=15t-4EqT7RDFnjkBJiW6yzEykpAcak4pW
download_large_googlefiles

fileName=spliced.csv
fileId=1hih3YOSummstHyVMRx9pcgZp2Bc3N-AH
download_large_googlefiles

## reproduce Figure 5
fileName=pbmc_8k_retained.rds
fileId=1YiYD0yt8dbkK0AC6jm9g2o85tNRsiryZ
download_large_googlefiles
wget -O ./data/neuron900_retained.rds "https://drive.google.com/uc?export=download&id=1cwz4KDW4KdMziYRKVTieyu4mqfSjsYL-"
wget -O ./data/BlueYellowColormaps_V1.RData "https://drive.google.com/uc?export=download&id=1WHv5il6L26hAii961wKfvf3AZFZenpx9"

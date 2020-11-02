if (dir.exists("./data")){
       next
} else {
        dir.create("./data")
}

website1 <- "http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neurons_900/neurons_900_raw_gene_bc_matrices_h5.h5"
website2 <- "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_raw_feature_bc_matrix.h5"
website3 <- "https://drive.google.com/uc?export=download&id=1Z1Vxzpu17kWwZGo6f2BMKo9eLjofmdrk"
download.file(website1,destfile="./data/pbmc_10k_raw.h5")
download.file(website2,destfile="./data/neurons_900_raw.h5")
download.file(website3,destfile="./data/kang.sampleC.h5")

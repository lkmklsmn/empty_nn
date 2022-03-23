runSeurat <-
function(counts,gene.use=FALSE,resolution){
        require(Seurat)
        tmp <- CreateSeuratObject(counts = counts)
        tmp <- tmp[-grep("^RPS", rownames(tmp)),]
        tmp <- tmp[-grep("^RPL", rownames(tmp)),]
        tmp <- NormalizeData(tmp)
        tmp <- FindVariableFeatures(tmp)
        if (!gene.use){
                gene.use <- VariableFeatures(tmp)
        } else {
                gene.use <- gene.use[gene.use %in% rownames(tmp)]
        }
        tmp <- ScaleData(tmp, features = gene.use)
        tmp <- RunPCA(tmp,features=gene.use)
        tmp <- FindNeighbors(tmp, dims = 1:10)
        tmp <- FindClusters(tmp, resolution = resolution)
        tmp <- RunTSNE(tmp, dims = 1:10, check_duplicates = FALSE)
        return(tmp)
}

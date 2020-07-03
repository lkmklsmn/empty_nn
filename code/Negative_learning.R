library(caret)
library(keras)

# Negative Learning ####
neg_learning <- function(counts,neg.threshold=100,n_splits=5,
                         iteration=3,batch_size=16,epoch=10,verbose=TRUE,training_verbose=0){
        if (nrow(counts) < ncol(counts)){
                stop(paste0("Please transpose counts matrix before running neg_learning\n",
                            "  rows are barcodes, columns are genes"))
        }
        # Calculate total reads
        n_counts <- Matrix::rowSums(counts)
        names(n_counts) <- rownames(counts)
        n_counts <- n_counts[n_counts>=10]
        
        # Negative samples & Unlabeled samples
        negative <- names(n_counts[which(n_counts <= neg.threshold)])
        if (verbose){
                print(paste0("negtive samples ",length(negative))) 
        }
        negative <- sample(negative,10000)
        unlabel <- names(n_counts[which(n_counts > neg.threshold)])
        if (verbose){
                print(paste0("unlabeled samples ",length(unlabel)))
                print(paste0("unlabeled samples were split into ",n_splits," pieces")) 
        }
        
        # Keep 2k most frequent genes in negative samples
        gene.use <- names(tail(sort(Matrix::colSums(counts[negative,])), 2000))
        counts_2k <- counts[names(n_counts),gene.use]
        
        # Normalization
        if (verbose){
                print("data normalization")
        }
        norm.counts.2k <- sweep(counts_2k, 1, n_counts, '/')
        
        # Create a dataframe to store predicted p-value
        cv_p <- data.frame(matrix(0, ncol = iteration, nrow = length(unlabel)))
        rownames(cv_p) <- unlabel
        colnames(cv_p) <- paste0("iteration",seq(1,ncol(cv_p)))
        
        if (verbose){
                print("start training")
        }
        # Negative-Unlabled Learning
        for (i in seq(iteration)){
                if (verbose){
                        print(paste0("iteration ",i))
                }
                # Split unlabled samples into {n_splits} splits
                df <- data.frame(total_counts = n_counts[unlabel],label=1)
                require(caret)
                df$folds<- createFolds(factor(df$label), k = n_splits, list = FALSE)
                
                for (j in seq(1,n_splits)){
                        # Training and Testing data
                        if (verbose){
                                print(paste0("training split ",j))
                        }
                        train <- c(negative,rownames(df[df$folds==j,]))
                        test <- rownames(df[df$folds!=j,])
                        x_train <- data.matrix(norm.counts.2k[train,])
                        x_test <- data.matrix(norm.counts.2k[test,])
                        y_train <- c(rep(0, length(negative)),rep(1,nrow(df[df$folds==j,])))
                        random_indices <- sample(1:length(y_train))
                        x_train <- x_train[random_indices,]
                        y_train <- y_train[random_indices]
                        y_train <- to_categorical(y_train)
                        
                        # NN
                        model_neg <- neg_create_model(input=dim(x_train)[2])
                        model_neg %>% fit(
                                x_train, y_train,
                                batch_size = batch_size,
                                epochs = epoch,
                                verbose = training_verbose,
                                validation_split = 0.2
                        )
                        y_test_pred <- model_neg %>% predict(x_test)
                        cv_p[test,i] <- y_test_pred[,2]
                }
        }
        # return a boolean vector showing cells or not, predictions, and total counts
        cv_p$mean.crossval <- apply(cv_p,1,mean)
        nn_bcs <- rownames(cv_p[cv_p$mean.crossval>0.5,])
        nn.keep <- rownames(counts) %in% nn_bcs
        return(list("nn.keep"=nn.keep,"prediction"=cv_p))
}

neg_create_model <- function(input){
        model <- keras_model_sequential()
        model %>% 
                layer_dense(units = 256, activation = 'relu', input_shape = input) %>% 
                layer_dropout(rate = 0.4) %>% 
                layer_dense(units = 128, activation = 'relu') %>%
                layer_dropout(rate = 0.3) %>%
                layer_dense(units = 2, activation = 'softmax')
        model %>% compile(
                loss = 'binary_crossentropy',
                optimizer = optimizer_rmsprop(),
                metrics = c('accuracy')
        )
        model
}

runSeurat <- function(counts,gene.use,resolution){
        tmp <- CreateSeuratObject(counts = counts)
        tmp <- tmp[-grep("^RPS", rownames(tmp)),]
        tmp <- tmp[-grep("^RPL", rownames(tmp)),]
        tmp <- NormalizeData(tmp)
        gene.use <- gene.use[gene.use %in% rownames(tmp)]
        tmp <- ScaleData(tmp, features = gene.use)
        tmp <- RunPCA(tmp,features=gene.use)
        tmp <- FindNeighbors(tmp, dims = 1:10)
        tmp <- FindClusters(tmp, resolution = resolution)
        tmp <- RunTSNE(tmp, dims = 1:10, check_duplicates = FALSE)
        tmp <- RunUMAP(object = tmp, dims = 1:10, min.dist = 0.75)
        return(tmp)
}


library(keras)

# Curri vs Non_curri learning ####
Empty_NN <- function(counts,expected=NA,threshold=NA,n.quantile=20,k=1,
                     curri.learning=TRUE,batch_size=16,epoch=10,verbose=TRUE){
        if (nrow(counts) < ncol(counts)){
                stop(paste0("Please transpose counts matrix before running EmptyNN\n",
                            "  rows are barcodes, columns are genes"))
        }
        # Calculate total reads ####
        n_counts <- Matrix::rowSums(counts)
        names(n_counts) <- rownames(counts)
        
        if (is.na(expected) & is.na(threshold)){
                stop("Please input expected no. of cells / threshold")
        }
        # Calculate cutoff of Cell Ranger 2.0 ####
        if (is.na(threshold)){
                o <- order(n_counts, decreasing = TRUE)
                top <- n_counts[head(o, n = expected)]
                threshold = round(quantile(top, 0.99)*0.1,0)
        } else {threshold = threshold}
        
        if (verbose){
                print(paste0("threshold is ",threshold))
        }
        assign("threshold",threshold, envir = .GlobalEnv) 
        
        # Remove cells with less than 10 reads & Keep 5k most frequent genes ####
        counts_5k <- counts[which(n_counts > 10),names(tail(sort(Matrix::colSums(counts)), 5000))]
        
        # Split barcodes into quantiles by total counts ####
        above <- which(n_counts >= threshold)
        below <- which(n_counts < threshold & n_counts >10)
        
        split_into_quantiles <- function(tmp,quantile=n.quantile){
                qs <- seq(0, 1, length = quantile+1)
                #split(names(tmp), cut(tmp, breaks=c(quantile(tmp, probs = qs))))
                split(names(tmp), .bincode(tmp, breaks=c(quantile(tmp, probs = qs))))
        }
        
        top_splits <- split_into_quantiles(n_counts[above])
        top_splits <- top_splits[rev(1:length(top_splits))]
        bottom_splits <- split_into_quantiles(n_counts[below])
        if (length(top_splits) >= length(bottom_splits)){
                min.q <- length(bottom_splits)
                top_splits <- split_into_quantiles(n_counts[above],quantile=min.q)
                top_splits <- top_splits[rev(1:length(top_splits))]
        } else {
                min.q <- length(top_splits)
                bottom_splits <- split_into_quantiles(n_counts[below],quantile=min.q)
        }
        
        min_samples <- min(unlist(lapply(c(top_splits,bottom_splits),length)))
        if (batch_size > min_samples){
                stop(paste0("Batch size is larger than samples in the quantile.\n",
                            "  Suggested maximum batch size is ",min_samples))
        }
        
        # Training and Testing
        n <- length(top_splits)
        test <- c(unlist(top_splits[c((n-k+1):n)]),unlist(bottom_splits[c((n-k+1):n)]))
        train  <- setdiff(rownames(counts_5k), test)
        x_train <- data.matrix(counts_5k[train,])
        if (verbose){
                print(paste0("training samples ",dim(x_train)[1]))
        }
        x_test <- data.matrix(counts_5k[test,])
        if (verbose){
                print(paste0("testing samples ",dim(x_test)[1]))
        }
        y_train <- rep(0, nrow(x_train))
        names(y_train) <- rownames(x_train)
        y_train[intersect(names(n_counts[above]), rownames(x_train))] <- 1
        
        random_indices <- sample(1:length(y_train))
        x_train <- x_train[random_indices,]
        y_train <- y_train[random_indices]
        
        
        y_train <- to_categorical(y_train)
        if (verbose){
                print(paste0("splitting into ",length(top_splits)," quantiles"))
                print(paste0("predicting bcs with total counts from ",min(n_counts[test])," to ",max(n_counts[test])))
        }
        
        if (curri.learning){
                print("Start curri learning")
                model_curri <- create_model(input=dim(x_train)[2])
                curri_accs <- lapply(1:(n-k), function(x){
                        if (verbose){
                                print(paste0("training quantile",x))
                        }
                        training_accs <- unlist(lapply(1:epoch, function(k){
                                top_barcodes <- sample(top_splits[[x]], batch_size)
                                bottom_barcodes <- sample(bottom_splits[[x]], batch_size)
                                batch_indices <- c(match(top_barcodes, rownames(x_train)), match(bottom_barcodes, rownames(x_train)))
                                y <- c(rep(1, batch_size), rep(0, batch_size))
                                r <- sample(1:length(batch_indices))
                                batch_indices <- batch_indices[r]
                                y <- y[r]
                                stats <- train_on_batch(model_curri, x_train[batch_indices,], 
                                                        to_categorical(y), class_weight = NULL, sample_weight = NULL)
                                stats[[2]]
                        }))
                        training_acc <- mean(training_accs)
                        if (verbose){
                                print(paste0("training accuracy: ",training_acc))
                        }
                })
                y_test_pred <- model_curri %>% predict(x_test)
                y_test_pred <- y_test_pred[,2]
                names(y_test_pred) <- test
        }
        else {
                # Non-curri
                if (verbose){
                        print("Start non-curri learning")
                }
                model_non_curri <- create_model(input=dim(x_train)[2])
                model_non_curri %>% fit(
                        x_train, y_train,
                        batch_size = batch_size*2,
                        epochs = epoch,
                        #steps_per_epoch = 10,
                        verbose = 2,
                        validation_split = 0.1
                )
                y_test_pred <- model_non_curri %>% predict(x_test)
                y_test_pred <- y_test_pred[,2]
                names(y_test_pred) <- test
        }
        max <- max(n_counts[test])
        nn_bcs <- c(names(y_test_pred[y_test_pred>0.5]), names(n_counts[n_counts>max]))
        nn.keep <- rownames(counts) %in% nn_bcs
        return(list("nn.keep"=nn.keep,"prediction"=y_test_pred))
}

create_model <- function(input){
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
        tmp <- RunTSNE(tmp, dims = 1:10,check_duplicates = FALSE)
        tmp <- RunUMAP(object = tmp, dims = 1:10, min.dist = 0.75)
        return(tmp)
}





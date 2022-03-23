emptynn <-
function(counts,threshold=100,k=5,iteration=5,
                    batch_size=16,epoch=10,verbose=TRUE,training_verbose=0){
        require(keras)
        require(Matrix)
        if (nrow(counts) < ncol(counts)){
                stop(paste0("Please transpose counts matrix before running EmptyNN\n",
                            "  rows are barcodes, columns are genes"))
        }
        # Calculate total reads
        n_counts <- Matrix::rowSums(counts)
        names(n_counts) <- rownames(counts)
        n_counts <- n_counts[n_counts>=10]

        # Negative samples & Unlabeled samples
        negative <- names(n_counts[which(n_counts <= threshold)])
        if (verbose){
                print(paste0("there are ",length(negative)," in P set"))
        }
        if (length(negative)>=10000){
                negative <- sample(negative,10000)
        } else {
                negative <- negative
        }

        unlabel <- names(n_counts[which(n_counts > threshold)])

        if (verbose){
                print(paste0("there are ",length(unlabel)," in U set"))
                print(paste0("Samples in U set were split into ",k," folds"))
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
                # Split unlabled samples into {k} folds
                df <- data.frame(total_counts = n_counts[unlabel],label=1)
                require(caret)
                df$folds<- createFolds(factor(df$label), k = k, list = FALSE)

                for (j in seq(1,k)){
                        # Training and Testing data
                        if (verbose){
                                print(paste0("training fold ",j))
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
        # return a boolean vector showing cells or not & predictions
        cv_p$mean.crossval <- apply(cv_p,1,mean)
        nn_bcs <- rownames(cv_p[cv_p$mean.crossval>0.5,])
        nn.keep <- rownames(counts) %in% nn_bcs
        return(list("nn.keep"=nn.keep,"prediction"=cv_p))
}

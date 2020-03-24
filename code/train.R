create_model <- function(){
        model <- keras_model_sequential()
        model %>% 
                layer_dense(units = 256, activation = 'relu', input_shape = c(5000)) %>% 
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

# Curri vs Non_curri learning ####
Empty_NN <- function(counts,threshold=200,curri.learning=TRUE,batch_size=32,epoch=10){
        # Calculate total reads ####
        n_counts <- Matrix::rowSums(counts)
        names(n_counts) <- rownames(counts)
        
        # Remove cells with less than 10 reads ####
        counts <- counts[which(n_counts > 10), ]
        n_counts <- n_counts[which(n_counts > 10)]
        
        # Keep 5k most frequent genes ####
        counts_5k <- counts[, names(tail(sort(Matrix::colSums(counts)), 5000))]
        
        # Split barcodes into quantiles by total counts ####
        above <- which(n_counts >= threshold)
        below <- which(n_counts < threshold)
        split_into_quantiles <- function(tmp){
                qs <- seq(0, 1, length = 21)
                split(names(tmp), cut(tmp, breaks=c(quantile(tmp, probs = qs))))
        }
        top_splits <- split_into_quantiles(n_counts[above])
        top_splits <- top_splits[rev(1:length(top_splits))]
        bottom_splits <- split_into_quantiles(n_counts[below])
        
        # Training and Testing
        test <- c(unlist(top_splits[c(17:20)]),unlist(bottom_splits[c(17:20)]))
        train  <- setdiff(rownames(counts_5k), test)
        x_train <- data.matrix(counts_5k[train,])
        print(paste0("training sample ",dim(x_train)[1]))
        x_test <- data.matrix(counts_5k[test,])
        print(paste0("testing sample ",dim(x_test)[1]))
        y_train <- rep(0, nrow(x_train))
        names(y_train) <- rownames(x_train)
        y_train[intersect(names(n_counts[above]), rownames(x_train))] <- 1
        y_train <- to_categorical(y_train)
        y_test <- global[test,2]
        print(paste0("predicting bcs with total counts from ",min(n_counts[test])," to ",max(n_counts[test])))

        if (curri.learning){
                print("Start curri learning")
                model_curri <- create_model()
                curri_accs <- lapply(1:16, function(x){
                        print(paste0("training quantile",x))
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
                        print(paste0("training accuracy: ",training_acc))
                })
                y_test_pred <- model_curri %>% predict(x_test)
                y_test_pred <- y_test_pred[,2]
                names(y_test_pred) <- test
                tl <- unclass(table(y_test_pred > 0.5, y_test))
                testing_acc <- sum(diag(tl))/sum(tl)
                print(paste0("testing accuracy: ",testing_acc))
        }
        else {
                # Non-curri 
                print("Start non-curri learning")
                model_non_curri <- create_model()
                model_non_curri %>% fit(
                        x_train, y_train,
                        batch_size = batch_size*2,
                        epochs = epoch,
                        #steps_per_epoch = 10,
                        verbose = 2,
                        validation_split = 0.1
                )
                y_test_pred <- model_non_curri %>% predict(x_test)
                tl <- unclass(table(y_test_pred[,2] > 0.5, y_test))
                testing_acc <- sum(diag(tl))/sum(tl) 
        }
        return(list(testing_acc,y_test_pred))
}
neg_create_model <-
function(input){
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

  ## Initial Setup
  
#  Loading all required libraries


library(keras)
library(caret)


## Setting the data up

#First we have to recover the treated data.


PosUTR <- scan(file="PosUTR",what="character")
PosUTR <- PosUTR[-which(duplicated(PosUTR))]

NegUTR <- scan(file="NegUTR",what="character")
NegUTR <- NegUTR[-which(duplicated(NegUTR))]


# Now we face an issue with our sequences. They are 1250 characters long strings, but R sees them as indivisible elements. We have to transform them before working with them, turning each 1250 characters long sequence into 1250 sequences of only one character.

print("hard")

posNum <- length(PosUTR) # Number of positive sequences
SeqUTR <- append(PosUTR,NegUTR)

n2 <- nchar(SeqUTR[1]) -1 # Number of divisions to make

secuenciasOH <- sapply(1 + 1*0:n2, function(z) substr(SeqUTR, z, z))  # Obtaining the split characters
df <- data.frame(secuenciasOH, "Value" = 0) # Saving them into a dataframe
indx <- sapply(df, is.factor) 
df[indx] <- lapply(df[indx], function(x) as.character(x)) # Factor --> char conversion

df[1:posNum,]$Value <- 1 # We have the value of the positive sequences to positive


# Since we are going to train five different models, we will divide our data into five different dataframes. Each of them will have the same number of positive and negative sequences. Since we dont
# have many positive sequences, we will train every model with every positive sequences. This way, every model will be trained with the same set of positive sequences, but a different set of negative sequences.

print("triangulo")

start = posNum + 1
end = start + posNum - 1
df.a <- df[c(1:posNum,start:end),]
df.b <- df[c(1:posNum,(start + posNum):(end + posNum)),]
df.c <- df[c(1:posNum,(start + (posNum*2)):(end + (posNum*2))),]
df.d <- df[c(1:posNum,(start + (posNum*3)):(end + (posNum*3))),]
df.e <- df[c(1:posNum,(start + (posNum*4)):(end + (posNum*4))),]

dfs <- list(df.a, df.b, df.c, df.d, df.e)


# Now we have to divide our data into Training and Test for our model training.

print("particiones")

output <- c("Value")
trains <- list()
tests <- list()

for (i in 1:5){
  partition <- createDataPartition(dfs[[i]][[output]],
                                     p = 0.8,
                                     list = FALSE,
                                     times = 1)
  trains <- append(trains, list(dfs[[i]][partition,]))
  tests <- append(tests, list(dfs[[i]][-partition,]))
}



# Once divided, we have to adapt the data to the format the model expects.

print("todatamatrix")

train.x <- list()
train.y <- list()

test.x <- list()
test.y <- list()

for (i in 1:5){

  train.x <- append(train.x, list(data.matrix(trains[[i]][,1:1250])))
  train.y <- append(train.y, list(data.matrix(trains[[i]][,1251])))

  test.x <- append(test.x, list(data.matrix(tests[[i]][,1:1250])))
  test.y <- append(test.y, list(data.matrix(tests[[i]][,1251])))

  train.x[[i]] <- array_reshape(train.x[[i]], c(nrow(trains[[i]]),1250,1))
  test.x[[i]] <- array_reshape(test.x[[i]], c(nrow(tests[[i]]),1250,1))
 }


# Now it is time to build our model. It follows the specifications mentioned in the thesis report. 

print("model")

  batch_size <- 125
  epochs <- 125 
  input_shape <- c(1250,1)
  learn_rate = 0.0001
  
  modelConvolu <- keras_model_sequential()
  modelConvolu %>% 
    layer_conv_1d(filters = 48, kernel_size = 75, activation = "relu", input_shape = input_shape)%>%
    layer_flatten() %>%
    layer_dense(units = 50, activation = 'relu') %>% 
    layer_dense(units = 1, activation="sigmoid") %>%
    
    summary(model)
  
  modelConvolu %>% compile(
    loss = loss_binary_crossentropy,
    optimizer = optimizer_nadam(lr = learn_rate),
    metrics = c('accuracy')
  )


# And now we train that model with the data obtained. After each training epoch, we will obtain certain statistics related to the sensitivity and specificity of the model, since they are not metrics that Keras can return from the model.

print("train")

for (i in 1:5){
  print(i)
  sens <- NULL
  spec <- NULL
  history <- NULL
 
  
   modelConvolu <- keras_model_sequential()
  modelConvolu %>% 
    layer_conv_1d(filters = 48, kernel_size = 75, activation = "relu", input_shape = input_shape)%>%
    layer_flatten() %>%
    layer_dense(units = 50, activation = 'relu') %>% 
    layer_dense(units = 1, activation="sigmoid")
  
  modelConvolu %>% compile(
    loss = loss_binary_crossentropy,
    optimizer = optimizer_nadam(lr = learn_rate),
    metrics = c('accuracy')
  )
  
    for (epoch in 1:epochs){
    historial <- modelConvolu %>% fit(
      x = train.x[[i]],
      y = train.y[[i]],
      epochs = 1,  
      batch_size = batch_size, 
      validation_data = list(test.x[[i]], test.y[[i]]),
      verbose = 2)
    
    history <- c(history, historial)
    
    ## Prediccion 
    pred <- modelConvolu %>% predict(test.x[[i]], batch_size = batch_size)
    prediction.y = round(pred)
    # Confusion matrix
    CM = table(prediction.y, test.y[[i]])
    if (length(CM) > 2){
      sens <- c(sens,sensitivity(CM))
      spec <- c(spec,specificity(CM))
    }
    else {  # This happens in the weird cases the model only predicts positive or negative
      if (sum(grepl(0,prediction.y)) > 0){
        sens <- c(sens, 0)
        spec <- c(spec, 1)
      }
      else{
        sens <- c(sens, 0)
        spec <- c(spec, 1)
      }
    }
  }
}


# Lastly, we will plot the graphs showing the results of the training of this model.


nums <- NULL
nums <- (1:(epochs))*2
valaccs <- NULL
accs <- NULL
valoss <- NULL
loss <- NULL
for (i in nums){ 
  valaccs <- c(valaccs,history[i]$metrics$val_acc) 
  accs <- c(accs,history[i]$metrics$acc) 
  loss <- c(loss,history[i]$metrics$loss) 
  valoss <- c(valoss,history[i]$metrics$val_loss) 
  }

plot(valoss, ylim = c(min(valoss,loss), max(valoss,loss)), cex = 0, ylab="loss", xlab="epochs")
lines(valoss, col="green")
lines(loss, col="blue")
legend("bottomleft", legend=c("training","validacion"), col=c("blue","green"), lty = 1)
dev.print(png, 'loss.png')

plot(valaccs, ylim = c(min(valaccs,accs), max(valaccs,accs)), cex = 0, ylab = "accuracy", xlab = "epochs")
lines(valaccs, col="green")
lines(accs, col="blue")
legend("topleft", legend=c("training","validacion"), col=c("blue","green"), lty = 1)
dev.print(png, 'acc.png')

plot(sens, ylim = c(min(sens, spec, sens + spec - 1), max(sens, spec, sens+spec-1)), cex = 0)
lines(sens,col="red")
lines(spec,col="blue")
lines((sens + spec) -1, col="orange")
legend("bottomright", legend=c("Specificity","Sensitivity", "Youden"), col=c("blue","red", "orange"), lty = 1)
dev.print(png, 'youden.png')



mean_acc = sum(tail(valaccs,n=7))/7 # Estimation of accuracy obtained
mean_acc
resultadoacc <- file("ResultAccuracy")
writeLines(toString(mean_acc), con=resultadoacc)




## Loading the libraries

library(biomaRt)
library(SetupSequences)
library(keras)
library(caret)

## Obtaining the sequences

library(tidyverse)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(GenomicRanges)

normalizeDF <- function(dataframe){
 mn <- min(dataframe)
 mx <- max(dataframe)
 
 return((dataframe - mn) / (mx - mn))
}

ml_table = readRDS(file = "../Data/ml_data_rf.rds") %>% 
  dplyr::select(class) %>% 
  rownames_to_column("id") %>% 
  tidyr::separate(id, c("gene", "seqnames", "start", "end", "strand", "tmp"), sep=":", remove=FALSE)


# converting to Granges object
ml_table_gr = makeGRangesFromDataFrame(ml_table, keep.extra.columns = TRUE)

# adding "chr" in front of seqnames
newStyle <- mapSeqlevels(seqlevels(ml_table_gr), "UCSC")
ml_table_gr <- renameSeqlevels(ml_table_gr, newStyle)

# Extract sequences using the package BSgenome
data = data.frame(
  id = ml_table$id,
  class = ml_table$class,
  sequence = BSgenome::getSeq(Hsapiens, ml_table_gr, as.character = TRUE) #as.character=TRUE does not work?
)
data$sequence <- as.character(data$sequence)
table(data$class)


# Next step is grouping all the non 3'UTR sequences into a group, so that we can compare 3'UTR sequences with non 3'UTR sequences regardless of the type.


old_data <- data # Needed for the valdation tests

data <- data[data$class == "ICE" | data$class == "UTR",]
data$detailclass = data$class
data$class = ifelse(data$class %in% "UTR", "UTR", "Non_3_UTR")  %>% as.factor() %>% as.integer()
data$class <- data$class - 1 #0 = No3UTR 1 = 3UTR
table(data$class)

# Validation test
test1 = table(data$class)["0"] == nrow(data) - table(old_data$class)["UTR"]; if (test1) print ("Correct 1") else print ("ERROR 1")
test2 = table(data$class)["1"] == table(old_data$class)["UTR"]; if (test2) print ("Correct 2") else print ("ERROR 2")



# We will separate these sequences into the "Positive" group (3'UTR) and the "Negative" group (non-3'UTR).
                                                                                             

PositiveSequences = data[data$class == 1,]
NegativeSequences = data[data$class == 0,]


# And save these sequences to a file.


write.csv(PositiveSequences,file="../Data/PositiveUTRSequences_untreated_sid", row.names = FALSE, quote = FALSE)
write.csv(NegativeSequences,file="../Data/NegativeUTRSequences_untreated_sid", row.names = FALSE, quote = FALSE)


## Treating the sequences

# First, we recover the sequences we obtained with ensembl.

PosSequences <- PositiveSequences
NegSequences <- NegativeSequences


# Next we start with the treatment of sequences. We have to:

# - Reverse them
# - Make sure all sequences are 250 characters long
# - One-hot encode them

# For that, we need some functions from the custom package developed for this project.

# We save the sequences for later use
PosSequences$rawseq <- PosSequences$sequence
NegSequences$rawseq <- NegSequences$sequence

# Example of sequence: AACCGT
PosSequences$sequence[1]  # Example of a positive sequence
NegSequences$sequence[1]  # Example of a negative sequence


old_seq = PosSequences$sequence[1]

# strReverse:  AACCGT --> TGCCAA
PosSequences$sequence <- strReverse(PosSequences$sequence)
NegSequences$sequence <- strReverse(NegSequences$sequence)

PosSequences$sequence[1]  # The positive sequence, once reversed
NegSequences$sequence[1]  # The negative sequence, once reversed

# Validation
doublereverse = strReverse(PosSequences$sequence[1])
test1 = doublereverse == old_seq; if (test1) print ("Correct 1") else print ("ERROR 1")

# padding_sequences: Inserts padding characters ("X") in sequences shorter than 250 characters, and trims sequences longer than 250 characters
PosSequences$sequence <- padding_sequences(PosSequences$sequence)
NegSequences$sequence <- padding_sequences(NegSequences$sequence)

PosSequences$sequence[1]  # The positive sequence, with some characters added so it can be 250 characters long
NegSequences$sequence[1]  # The negative sequence, trimmed so it becomes 250 characters long

# Validation
test1 = all(sapply(list(min(nchar(PosSequences$sequence)),max(nchar(PosSequences$sequence)), mean(nchar(PosSequences$sequence))), function(x) x == 250))
test2 = all(sapply(list(min(nchar(NegSequences$sequence)),max(nchar(NegSequences$sequence)), mean(nchar(NegSequences$sequence))), function(x) x == 250))
if (test1) print ("Correct 1") else print ("ERROR 1")
if (test2) print ("Correct 2") else print ("ERROR 2")

# to_onehot:   TGCCAA --> 00010010010010001000    (A = 1000, C = 0100, G = 0010, T = 0001)
PosSequences$sequence <- to_onehot(PosSequences$sequence)
NegSequences$sequence <- to_onehot(NegSequences$sequence)

PosSequences$sequence[1]  # The positive sequence, reversed, padded and one-hot encoded
NegSequences$sequence[1]  # The negative sequence, reversed, trimmed and one-hot encoded

# Validation
test1 = all(sapply(list(min(nchar(PosSequences$sequence)),max(nchar(PosSequences$sequence)), mean(nchar(PosSequences$sequence))), function(x) x == 1250))
test2 = all(sapply(list(min(nchar(NegSequences$sequence)),max(nchar(NegSequences$sequence)), mean(nchar(NegSequences$sequence))), function(x) x == 1250))
if (test1) print ("Correct 1") else print ("ERROR 1")
if (test2) print ("Correct 2") else print ("ERROR 2")

# Final check: Since some of the sequences might be corrupt, we will delete the ones that are corrupt, if they exist
# We can know if a sequence is corrupt by looking for the W character. When one-hot encoding, we encoded everything that was not an A, T, C or G with a "W"

if (any(grepl("W",PosSequences$sequence))){
print("Found 'W'")
PosSequences <- PosSequences[-which(grepl("W",PosSequences$sequence)),]
}

if (any(grepl("W",NegSequences$sequence))){
print("Found 'W'")
NegSequences <- NegSequences[-which(grepl("W",NegSequences$sequence)),]
}

# Once treated, we can save them in a file for later use.


fileCon<-file("../Data/Pos3UTR")
write.table(PosSequences$sequence,file = fileCon, quote=FALSE, row.names = FALSE, col.names = FALSE)
fileCon<-file("../Data/Neg3UTR")
write.table(NegSequences$sequence,file = fileCon, quote=FALSE, row.names = FALSE, col.names = FALSE)


## Training the network
# First we have to recover the treated data.


PosUTR <- scan(file="../Data/Pos3UTR",what="character")
length(PosUTR)
PosSequences_nodup <- PosSequences[-which(duplicated(PosUTR)),]
PosUTR <- PosUTR[-which(duplicated(PosUTR))]
length(PosUTR)

NegUTR <- scan(file="../Data/Neg3UTR",what="character")
length(NegUTR)
NegSequences_nodup <- NegSequences[-which(duplicated(NegUTR)),]
NegUTR <- NegUTR[-which(duplicated(NegUTR))]
length(NegUTR)

# Validation test
test1 = all(NegSequences_nodup$sequence == NegUTR) 
test2 = all(PosSequences_nodup$sequence == PosUTR)
if (test1) print ("Correct 1") else print ("ERROR 1")
if (test2) print ("Correct 2") else print ("ERROR 2")

# Now we face an issue with our sequences. They are 1250 characters long strings, but R sees them as indivisible elements. We have to transform them before working with them, turning each 1250 characters long sequence into 1250 sequences of only one character.


posNum <- length(PosUTR) # Number of positive sequences
SeqUTR <- append(PosUTR,NegUTR)

old_seq = SeqUTR[1] # For validation tests

n2 <- nchar(SeqUTR[1]) -1 # Number of divisions to make

secuenciasOH <- sapply(1 + 1*0:n2, function(z) substr(SeqUTR, z, z))  # Obtaining the split characters
df <- data.frame(secuenciasOH, "Value" = 0) # Saving them into a dataframe
indx <- sapply(df, is.factor) 
df[indx] <- lapply(df[indx], function(x) as.character(x)) # Factor --> char conversion

df[1:posNum,]$Value <- 1 # We have the value of the positive sequences to positive

table(df$Value)

# Validation tests
test1 = paste(unlist(secuenciasOH[1,]),sep="", collapse="") == old_seq; if (test1) print ("Correct 1") else print ("ERROR 1")
test2 = table(df$Value)["0"] == length(NegUTR); if (test2) print ("Correct 2") else print ("ERROR 2")
test3 = table(df$Value)["1"] == length(PosUTR); if (test3) print ("Correct 3") else print ("ERROR 3")
test4 = all(length(PosUTR) + length(NegUTR) == length(SeqUTR), length(SeqUTR) == nrow(df)); if (test4) print ("Correct 4") else print ("ERROR 4")

# Since we are going to train five different models, we will divide our data into five different dataframes. Each of them will have the same number of positive and negative sequences. Since we don't have many positive sequences, we will train every model with every positive sequences. This way, every model will be trained with the same set of positive sequences, but a different set of negative sequences.

# Select all the negative sequences
set.seed(1234)
chosenSeqs <- sample(nrow(df) - posNum, size = posNum*5, replace = FALSE, prob = NULL) # Selected rows from 0 to max-posNum
chosenSeqs <- chosenSeqs + posNum

sets <- list()

sets[[1]] <- 1:posNum
sets[[2]] <- (posNum+1):(posNum*2)
sets[[3]] <- ((posNum*2)+1):(posNum*3)
sets[[4]] <- ((posNum*3)+1):(posNum*4)
sets[[5]] <- ((posNum*4)+1):(posNum*5)

df.a <- df[c(1:posNum,chosenSeqs[sets[[1]]]),]
df.b <- df[c(1:posNum,chosenSeqs[sets[[2]]]),]
df.c <- df[c(1:posNum,chosenSeqs[sets[[3]]]),]
df.d <- df[c(1:posNum,chosenSeqs[sets[[4]]]),]
df.e <- df[c(1:posNum,chosenSeqs[sets[[5]]]),]

dfs <- list(df.a, df.b, df.c, df.d, df.e) # 

totalSeq = (nrow(df.a) + nrow(df.b) + nrow(df.c) + nrow(df.d) + nrow(df.e)) 
length(SeqUTR) - totalSeq # Number of sequences that do not appear in the dataframes
test1 = all(sapply(list(table(df.a$Value),table(df.b$Value), table(df.c$Value), table(df.d$Value)), function(x) x == table(df.e$Value)))
test2 = table(df.a$Value)["0"] == posNum
test3 = table(df.a$Value)["1"] == posNum
if (test1) print ("Correct 1") else print ("ERROR 1")
if (test2) print ("Correct 2") else print ("ERROR 2")
if (test3) print ("Correct 3") else print ("ERROR 3")


start = posNum + 1
end = start + posNum - 1
df.a <- df[c(1:posNum,start:end),]
df.b <- df[c(1:posNum,(start + posNum):(end + posNum)),]
df.c <- df[c(1:posNum,(start + (posNum*2)):(end + (posNum*2))),]
df.d <- df[c(1:posNum,(start + (posNum*3)):(end + (posNum*3))),]
df.e <- df[c(1:posNum,(start + (posNum*4)):(end + (posNum*4))),]

dfs <- list(df.a, df.b, df.c, df.d, df.e) # 

# Validation tests
totalSeq = (nrow(df.a) + nrow(df.b) + nrow(df.c) + nrow(df.d) + nrow(df.e)) 
length(SeqUTR) - totalSeq # Number of sequences that do not appear in the dataframes
test1 = all(sapply(list(table(df.a$Value),table(df.b$Value), table(df.c$Value), table(df.d$Value)), function(x) x == table(df.e$Value)))
test2 = table(df.a$Value)["0"] == posNum
test3 = table(df.a$Value)["1"] == length(SeqUTR) - table(df.a$Value)["0"]
if (test1) print ("Correct 1") else print ("ERROR 1")
if (test1) print ("Correct 2") else print ("ERROR 2")
if (test1) print ("Correct 3") else print ("ERROR 3")

# Now we have to divide our data into Training and Test for our model training.


output <- c("Value")
trains <- list()
tests <- list()
usedForTrainPos <- list()
usedForTrainNeg <- list()

set.seed(1234)

for (i in 1:5){
partition <- createDataPartition(dfs[[i]][[output]],
                                p = 0.8,
                                list = FALSE,
                                times = 1)
trains <- append(trains, list(dfs[[i]][partition,])) # 28210
tests <- append(tests, list(dfs[[i]][-partition,])) # 7852
# Sum of both: 35262

partPos <- partition[partition <= posNum]
partNeg <- partition[partition > posNum]
partNeg <- partNeg - posNum

usedForTrainPos <- append(usedForTrainPos,list(unlist(partPos)))
usedForTrainNeg <- append(usedForTrainNeg,list(unlist(chosenSeqs[sets[[i]][partNeg]])))
}

# Validation test

total = 0

for (i in 1:5) {
total = total + nrow(trains[[i]]) + nrow(tests[[i]])
}

test1 = total == totalSeq; if (test1) print ("Correct 1") else print ("ERROR 1")
test2 = sum(duplicated(usedForTrainPos)) == 0; if (test2) print ("Correct 2") else print ("ERROR 2")
test3 = sum(duplicated(usedForTrainNeg)) == 0; if (test3) print ("Correct 3") else print ("ERROR 3")

#usedForTrainPos <- usedForTrainPos[-which(duplicated(usedForTrainPos))] # Only the positive genes could be duplicated

# Checking the type distribution of the train sets

for (i in 1:5){
negtr <- usedForTrainNeg[[i]]
print(min(negtr))
print(max(negtr))
print(table(NegSequences_nodup$detailclass[negtr]))
}

# Once divided, we have to adapt the data to the format the model expects.

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

# We have to prepare the dataframe that will contain the sequences and the output of the models.

# Next step is making a dataframe of every control gene and predicting positive or negative
genesN <- NegSequences_nodup
genesP <- PosSequences_nodup

##### We have to make some transformations to these secuences

n2 <- nchar(genesN$sequence[1]) - 1

secuencesOH <- sapply(1 + 1*0:n2, function(z) substr(genesN$sequence, z, z))  
dfn <- data.frame(secuencesOH, class=genesN$detailclass, id=genesN$id) 
indx <- sapply(dfn, is.factor) 
dfn[indx] <- lapply(dfn[indx], function(x) as.character(x))  #

secuencesOH <- sapply(1 + 1*0:n2, function(z) substr(genesP$sequence, z, z))  
dfp <- data.frame(secuencesOH, class=genesP$detailclass, id=genesP$id) 
indx <- sapply(dfp, is.factor) 
dfp[indx] <- lapply(dfp[indx], function(x) as.character(x))  

# We have got a dataframe with the ensembl_id, sequence and hgnc symbol in genesP and genesN, and a dataframe with the split sequence and hgnc symbol in dfp and dfn

finalSequencesN <- data.matrix(dfn[,1:1250])
finalSequencesN <- array_reshape(finalSequencesN,c(nrow(finalSequencesN),1250,1))
finalSequencesP <- data.matrix(dfp[,1:1250])
finalSequencesP <- array_reshape(finalSequencesP,c(nrow(finalSequencesP),1250,1))

genesP$model1 <- 0
genesP$model2 <- 0
genesP$model3 <- 0
genesP$model4 <- 0
genesP$model5 <- 0

genesN$model1 <- 0
genesN$model2 <- 0
genesN$model3 <- 0
genesN$model4 <- 0
genesN$model5 <- 0

# Now it is time to build our model. It follows the specifications mentioned in the thesis report. 

batch_size <- 125
epochs <- 25 
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

config = tensorflow::tf$ConfigProto(gpu_options = list(allow_growth = TRUE))
sess = tensorflow::tf$Session(config = config)
keras::k_set_session(session = sess)

#with(tf$device("/device:GPU:0"), {

models <- list()

for (i in 1:5){

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
}
## Prediction 
pred <- modelConvolu %>% predict(finalSequencesP, batch_size = batch_size)
genesP[,5+i] = pred

pred <- modelConvolu %>% predict(finalSequencesN, batch_size = batch_size)
genesN[,5+i] = pred

models <- append(models, modelConvolu)
}

#})

# First, we have to obtain the data needed, for example the weights ¿or the biases?. For that, we need to access its layers 

lrp <- function(model, dataframe, indexSeq){
  print(paste0("Sequence number ",indexSeq))
  layers <- model$layers
  input <- as.matrix(as.numeric(dataframe[indexSeq,1:1250])) # posNum+10
  
  ## STEP 1: CREATE DATA STRUCTURES
  activs <- list()
  weights <- list()
  biases <- list()
  
  for (i in 1:length(layers)){
    lay <- layers[[i]]
    layer_name <- lay$name
    
    intermediate_layer_model <- keras_model(inputs = modelConvolu$input,
                                            outputs = get_layer(modelConvolu, layer_name)$output)
    intermediate_output <- predict(intermediate_layer_model, array_reshape(input, c(1,1250,1)))
    
    if (length(lay$weights) > 0){ # If not a flatten layer
      activs[[i]] <- intermediate_output
      weights[[i]] <- lay$get_weights()[[1]]
      biases[[i]] <- lay$get_weights()[[2]]
    }
  }
  
  ## STEP 2: OUTPUT LAYER
  
  Rk <- list()
  Rk[["output"]] <- as.numeric(activs[[length(layers)]])
  
  if (Rk[["output"]] < 0.5) # If the output is lesser than 0.5 it means it is a negative sequence, and we have to compute the relevance towards the negative class
    Rk[["output"]] <- 1 - Rk[["output"]]
  print(Rk[["output"]])
  ## STEP 3: DENSE LAYER
  
  pos_weights <- weights[[4]]
  for (i in 1:length(weights[[4]]))
    pos_weights[i] <- max(weights[[4]][i],0) # Not using the apply family functions so it does not alter the shape
  
  ##  if (min(pos_weights) == 0) print("CORRECT 1") else print ("ERROR 1")
  
  denom = 0
  
  for (i in 1:length(activs[[3]])) { # Calculating the denominator
    value <- activs[[3]][i] * pos_weights[i]
    denom <- denom + value
  }
  
  for (i in 1:length(activs[[3]])) { # Calculating each numerator and, with it, the relevance of each neuron of the dense layer
    numer <- activs[[3]][i] * pos_weights[i]
    if (denom != 0) {
      division <- numer/denom
    } else {
      division <- 0
    }
    Rk[["dense"]][i] <- division * Rk[["output"]][1]
  }
  
  ## STEP 4: CONVOLUTIONAL LAYER
  
  pos_weights <- weights[[3]]
  for (i in 1:length(weights[[3]]))
    pos_weights[i] <- max(weights[[3]][i],0) # Not using the apply family functions so it does not alter the shape
  
  ##  if (min(pos_weights) == 0) print("CORRECT 1") else print ("ERROR 1")
  
  for (i in 1:length(activs[[1]]))
    Rk[["conv"]][i] = 0
  
  for (k in 1:length(Rk[["dense"]])) { # For each neuron in the layer analyzed before we have to repeat the main cycle
    denom <- 0
    
    # pos_weights[j,k]
    #  individual_weights = length(pos_weights) / length(Rk[["dense"]])
    
    for (j in 1:length(activs[[1]])) { # For each neuron in the actual layer, calculate denominator
      value <- activs[[1]][j] * pos_weights[j,k]
      denom <- denom + value
    }
    
    for (j in 1:length(activs[[1]])) {# For each neuron in the actual layer, calculate numerator and relevance j<-k
      numer <- activs[[1]][j] * pos_weights[j,k]
      if (denom != 0) {
        division <- numer/denom
      } else {
        division <- 0
      }
      Rk[["conv"]][j] <- Rk[["conv"]][j] + (division * Rk[["dense"]][k])
    }
  }
  
  ## STEP 5: INPUT LAYER
  
  pos_weights <- weights[[1]]
  
  for (i in 1:length(weights[[1]]))
    pos_weights[i] <- max(weights[[1]][i],0) # Not using the apply family functions so it does not alter the shape
  
  ##  if (min(pos_weights) == 0) print("CORRECT 1") else print ("ERROR 1")
  
  
  used <- list()
  for (i in 1:1250){
    Rk[["input"]][i] = 0
    used[["input"]][i] = 0
  }
  
  for (k in 1:length(activs[[1]])){ # For each one of the convolution cells
    #  if (Rk[["conv"]][k] == 0)
    #   next
    ## First we have to obtain the list of input nodes that contribute to this value. We know that there are as many nodes as the size of the filters:
    num_neurons <- length(weights[[1]][,,1])
    
    ## Now we have to know which are these neurons. For that reason, we have to obtain the position of this cell relative to the convolution filter
    # Since a%%b returns a value in [0,b-1] and we want a value in [1,b], we have to do it this way:
    cells_per_kernel <- length(activs[[1]][,,1])
    k_relative <- ((k-1)%%cells_per_kernel)+1
    # And now we obtain the list of input values that affected this particular cell
    j_neurons <- k_relative:(k_relative + num_neurons - 1)
    
    ## To obtain the weights, we have to know to which kernel this cell belongs to
    kernel <- ceiling(k/cells_per_kernel)
    
    ## And now we can obtain said weights
    kernel_weights <- pos_weights[,,kernel]
    
    ## Now we can calculate the denominator for this filter, following the formula
    denom <- sum(input[j_neurons] * kernel_weights)
    
    ## And, for each j neuron, we calculate its numerator, and add its value to its relevance
    for (j in 1:length(j_neurons)){
      numer <- input[j_neurons[j]] * kernel_weights[j] * Rk[["conv"]][k]
      if (denom != 0) {
        division <- numer/denom
      } else {
        division <- 0
      }
      Rk[["input"]][j_neurons[j]] = Rk[["input"]][j_neurons[j]] + division
      used[["input"]][j_neurons[j]] = used[["input"]][j_neurons[j]] + 1
    }
  }
  
  ## STEP 6: RETURNING RESULT
  
 # print(sum(Rk[["input"]]) == sum(Rk[["output"]]))
 # print(paste("Input", sum(Rk[["input"]])))
 # print(paste("Output", sum(Rk[["output"]])))
 # print(paste("Resta", sum(Rk[["output"]]) - sum(Rk[["input"]])))
  
  
  return(Rk[["input"]] / used[["input"]])
  #return(Rk[["input"]])
  
}

model <- models[[1]]

relevancesP <- genesP
relevancesN <- genesN

posS <- df[df$Value == 1,] # Quizá cambiar al formato del BasicExample, con el which
negS <- df[df$Value == 0,]

for (i in 1:250){
  relevancesP[[paste0("relevance",i)]] <- 0
  relevancesN[[paste0("relevance",i)]] <- 0
}

num_seq <- 500

for (i in 1:num_seq){
  print(paste0("Vamos con la secuencia positiva ", i))
  relevances <- lrp(modelConvolu, posS, i)
  
  for (j in 1:250){
    k = 1 + (5*(j - 1))
    relevancesP[i,paste0("relevance", j)] <- relevances[k] + relevances[k + 1] + relevances[k + 2] + relevances[k + 3] + relevances[k + 4]
  }
  
  write.csv(relevancesP[1:i,],"../Results/Relevances/posRelevances", row.names = FALSE)
}

for (i in 1:num_seq){
  print(paste0("Vamos con la secuencia negativa ", i))
  relevances <- lrp(modelConvolu, negS, i)
  
  for (j in 1:250){
    k = 1 + (5*(j - 1))
    relevancesN[i,paste0("relevance", j)] <- relevances[k] + relevances[k + 1] + relevances[k + 2] + relevances[k + 3] + relevances[k + 4]
  }
  
  write.csv(relevancesN[1:i,],"../Results/Relevances/negRelevances", row.names = FALSE)
}




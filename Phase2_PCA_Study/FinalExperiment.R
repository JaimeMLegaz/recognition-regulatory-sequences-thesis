# Requisitos: DATOS/ExperimentoUTR/SeqUTR
#             Results/experimentos

# Devuelve: history.rds, sens.rds, spec.rds

library(keras)
library(biomaRt)
library(caret)
library(beepr)
library(pROC)
library(SetupSequences)

set.seed(1337)

## Obtención de las secuencias
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

genesP <- scan(file = "DATOS/ExperimentoUTR/GenesPositivos", what = "raw") # Contains only the positive genes
genesN <- scan(file = "DATOS/ExperimentoUTR/nuevosGenes", what = "raw") # Contains ALL genes (later we will remove the positives)

genesN <- genesN[-which(genesN %in% genesP)] # We remove the negative genes present in the positive list

secuenciasPos <- getBM(attributes = c("ensembl_gene_id","5utr", "hgnc_symbol"),
                       filters = "hgnc_symbol",
                       values = genesP,
                       mart = ensembl)
secuenciasNeg<- getBM(attributes = c("ensembl_gene_id","5utr","hgnc_symbol"),
                      filters = "hgnc_symbol",
                      values = genesN,
                      mart = ensembl)

# We will now remove all sequences that ensembl couldn't give us, indicated by the "Sequence unavailable" string. And again, just in case, will remove the genes from the negatives list that also appear in the positives list
secuenciasPos <- secuenciasPos[!(secuenciasPos$`5utr` == "Sequence unavailable"),] # 3657
secuenciasNeg <- secuenciasNeg[!(secuenciasNeg$`5utr` == "Sequence unavailable"),] # 3094

# We remove positive genes from the negative list (just to be 100% sure)
secuenciasNeg <- secuenciasNeg[-which(secuenciasNeg$ensembl_gene_id %in% secuenciasPos$ensembl_gene_id),]

# After making sure twice that there are no positive genes in the negative list, we will now check how many sequences appear in both lists:

pos <- secuenciasPos
neg <- secuenciasNeg

shared <- data.frame(positive=integer(), negative=integer())

pos$shared <- 0

for (i in 1:nrow(pos)){
  print(i)
  for (j in 1:nrow(neg)){
    if (pos[i,]$`5utr` == neg[j,]$`5utr`){
      pos[i,]$shared <- pos[i,]$shared + 1
      shared[nrow(shared)+1,] <- list(i,j)
    }
  }
}

# Esto borrarlo en un futuro, es solo para copia de seguridad
write.csv(pos,file="DATOS/dataframePos", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.csv(neg,file="DATOS/dataframeNeg", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.csv(shared,file="DATOS/dataframeShared", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Now we remove the negative sequences that appear in the positive genes list

copiaSegNeg <- secuenciasNeg # COPIA DE SEGURIDADD, BORRAR DESPUES

negShared_nonDup <- shared$negative
negShared_nonDup <- negShared_nonDup[-which(duplicated(negShared_nonDup))]

secuenciasNeg <- secuenciasNeg[-negShared_nonDup,]

# We now fix the sequences

secuenciasPos$'5utr' <- strReverse(secuenciasPos$'5utr')
secuenciasNeg$'5utr' <- strReverse(secuenciasNeg$'5utr')

secuenciasPos$'5utr' <- padding_sequences(secuenciasPos$'5utr')
secuenciasNeg$'5utr' <- padding_sequences(secuenciasNeg$'5utr')

# Hasta aquí parece que todo va bien, todas las secuencias de UTR tienen 250 caracteres

secuenciasPos$'5utr' <- to_onehot(secuenciasPos$'5utr')
secuenciasNeg$'5utr' <- to_onehot(secuenciasNeg$'5utr')

seqPos <- secuenciasPos$'5utr'
seqNeg <- secuenciasNeg$'5utr'

# We now save the data in a file

seqUTR <- append(seqPos,seqNeg) # 3094 pos, 3094 neg, en ese orden

fileCon<-file("DATOS/ExperimentoUTR/SeqUTR_2")
write.table(seqUTR,file = fileCon, quote=FALSE, row.names = FALSE, col.names = FALSE)
write.csv(secuenciasNeg,file = "DATOS/ExperimentoUTR/dataframeSeqNeg", quote=FALSE, row.names = FALSE)

####  ---  ####


# Vamos a obtener las secuencias de los promotores de nuestros genes positivos y negativos que tenemos guardadas en un archivo.
seqUTR <- scan(file = "DATOS/ExperimentoUTR/SeqUTR",what="character") 

# Vamos a modificar estas secuencias para que tengan la forma que buscamos:

n2 <- nchar(seqUTR[1]) -1#/ 4 -1

secuenciasOH <- sapply(1 + 1*0:n2, function(z) substr(seqUTR, z, z))  # Obtenemos las cadenas divididas en subcadenas de longitud 4 (codificación one-hot)
df <- data.frame(secuenciasOH, "Valor" = 0) # Para checkear que es correcto, df[1,]
indx <- sapply(df, is.factor) 
df[indx] <- lapply(df[indx], function(x) as.character(x))  # Convertimos todos los factores en enteros

# Ahora mismo la columna "Valor" marca que todas son secuencias negativas Sin embargo, de la 1 a la 3657 son positivas Vamos a cambiarlo:
df[1:3657,]$Valor <- 1


# Ahora tenemos que dividir todo ese dataframe en 5 dataframes distintos. Los cinco comparten el set de positivos
start = 3658
end = start + 3656
df.a <- df[c(1:3657,start:end),]
df.b <- df[c(1:3657,(start + 3657):(end + 3657)),]
df.c <- df[c(1:3657,(start + (3657*2)):(end + (3657*2))),]
df.d <- df[c(1:3657,(start + (3657*3)):(end + (3657*3))),]
df.e <- df[c(1:3657,(start + (3657*4)):(end + (3657*4))),]

dfs <- list(df.a, df.b, df.c, df.d, df.e)

# Ahora tendremos que dividir nuestros datos en Training y Test, PARA CADA UNO DE LOS CINCO DATAFRAMES

prom.var.salida <- c("Valor")
trains <- list()
tests <- list()



for (i in 1:5){
  partition <- createDataPartition(dfs[[i]][[prom.var.salida]],
                                   p = 0.8,
                                   list = FALSE,
                                   times = 1)
  trains <- append(trains, list(dfs[[i]][partition,]))
  tests <- append(tests, list(dfs[[i]][-partition,]))
}

# Now we have "trains" containing the 5 training sets, and "tests" containing the 5 test sets


# Una vez divididos, tendremos que adaptar su estructura a la que quiere el modelo de Convolución
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

# Now we have two lists containing all blocks of training and test, for our 5 experiments

# Next step is making a dataframe of every control gene and predicting positive or negative
genesN <- secuenciasNeg
genesP <- secuenciasPos
 ##### We have to make some transformations to these secuences

n2 <- nchar(genesN$`5utr`[1]) - 1

secuenciasOH <- sapply(1 + 1*0:n2, function(z) substr(genesN$`5utr`, z, z))  # Obtenemos las cadenas divididas en subcadenas de longitud 4 (codificación one-hot)
dfn <- data.frame(secuenciasOH, gene=genesN$hgnc_symbol) # Para checkear que es correcto, df[1,]
indx <- sapply(dfn, is.factor) 
dfn[indx] <- lapply(dfn[indx], function(x) as.character(x))  # Convertimos todos los factores en enteros

secuenciasOH <- sapply(1 + 1*0:n2, function(z) substr(genesP$`5utr`, z, z))  # Obtenemos las cadenas divididas en subcadenas de longitud 4 (codificación one-hot)
dfp <- data.frame(secuenciasOH, gene=genesP$hgnc_symbol) # Para checkear que es correcto, df[1,]
indx <- sapply(dfp, is.factor) 
dfp[indx] <- lapply(dfp[indx], function(x) as.character(x))  # Convertimos todos los factores en enteros

# We have got a dataframe with the ensembl_id, sequence and hgnc symbol in genes, and a dataframe with the split sequence and hgnc symbol in df

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


for (i in 1:5){
  prueba <- i
  
  modelConvolu <- NULL
  history <- NULL
  sens <- NULL
  spec <- NULL
  
  history = NULL
  
  batch_size <- 125 # Arbitrario
  epochs <- 125 
  input_shape <- c(1250,1)
  learn_rate = 0.0001
  
  modelConvolu <- keras_model_sequential()
  modelConvolu %>% 
    layer_conv_1d(filters = 192, kernel_size = 75, activation = "relu", input_shape = input_shape)%>%
    layer_max_pooling_1d() %>%
    layer_dropout(0.25) %>%
    layer_flatten() %>%
    layer_dense(units = 50, activation = 'relu') %>% 
    layer_dense(units = 1, activation="sigmoid") %>% # 1 nodo y sigmpoide clasifica 0 o 1
    
    summary(model)
  
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
  ## Prediccion 
  pred <- modelConvolu %>% predict(finalSequencesP, batch_size = batch_size)
  genesP[,3+i] = pred
  
  pred <- modelConvolu %>% predict(finalSequencesN, batch_size = batch_size)
  genesN[,3+i] = pred
}

genes$totalSum <- genes$model1 + genes$model2 + genes$model3 + genes$model4 + genes$model5 

write.csv(genes[,c(3,4,5,6,7,8,9)],file="Results/finalEvaluation", row.names = FALSE, col.names = FALSE, quote = FALSE)
genesSpecials <- genes[which(genes$totalSum >= 4),]
write.csv(genesSpecials[,c(3,4,5,6,7,8,9)],file="Results/finalEvaluation>4", row.names = FALSE, col.names = FALSE, quote = FALSE)

genesNorep <- genesSpecials[-which(duplicated(genesSpecials$hgnc_symbol)),]

start = 3658
end = start + 3656
genes.a <- genes[1:3657,]
genes.b <- genes[(start):(end),]
genes.c <- genes[(start + (3657)):(end + (3657)),]
genes.d <- genes[(start + (3657*2)):(end + (3657*2)),]
genes.e <- genes[(start + (3657*3)):(end + (3657*3)),]

genesSp <- NULL
genesSp <- c(genesSp,genes.a[which(genes.a$model1 == 1),3])
genesSp <- c(genesSp,genes.b[which(genes.b$model2 == 1),3])
genesSp <- c(genesSp,genes.c[which(genes.c$model3 == 1),3])
genesSp <- c(genesSp,genes.d[which(genes.d$model4 == 1),3])
genesSp <- c(genesSp,genes.e[which(genes.e$model5 == 1),3])
write.csv(genesSp,file="Results/finalEvaluationSpecials", row.names = FALSE, col.names = FALSE, quote = FALSE)

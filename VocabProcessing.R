vocabMarginalize <- function(vocabList, membList, fileName, wordNames, nWords, sortByMemb){
nDays <- length(vocabList)

vocabDist <- vocabList[[1]]
membVect <- membList[[1]] 
for(i in 2:nDays){
  
  vocabDist <- vocabDist + vocabList[[i]]
  membVect <-  membVect+membList[[i]]  
}

vocabWordList <- list()
vocabProbList <- list()
for(topic in 1: dim(vocabDist)[[1]]){
  
  topicVocabOrder <- order(vocabDist[topic,],decreasing = T)
  
  vocabWordList[[topic]] <- wordNames[topicVocabOrder[1:nWords]]
  vocabProbList[[topic]] <- vocabDist[topic, topicVocabOrder[1:nWords]]/ sum(vocabDist[topic,])
}
  
vocabWordList <- do.call("rbind", vocabWordList)
vocabProbList <- do.call("rbind", vocabProbList)
vocabWordList <- cbind(vocabWordList, as.vector(membVect))


vocabWordList <- as.data.frame(vocabWordList)
names(vocabWordList) <- 1:dim(vocabDist)[[1]]

vocabProbList <- as.data.frame(vocabProbList)
names(vocabProbList) <- 1:dim(vocabDist)[[1]]

if(sortByMemb){
  vocabWordList <- vocabWordList[,order(membVect, decreasing = T)]
  vocabProbList <- vocabProbList[,order(membVect, decreasing = T)]
}
write.csv(vocabWordList, paste(fileName,"Words.csv",sep = ""))
write.csv(vocabProbList, paste(fileName,"Probs.csv",sep = ""))
return(list(vocabWordList, vocabProbList, membVect, vocabDist))
}
test<-  vocabMarginalize(vocab[[4]], vocab[[5]], fileName = "Test",wordNames = wordList[,2], nWords = , sortByMemb = T)

fileList <-list.files("./DataRuns/")
namesFiles <- strsplit(fileList,".", fixed = T)
namesFiles <- do.call("rbind", namesFiles)
namesFiles <- namesFiles[,1]
for(i in 1:length(fileList)){
  print(fileList[i])
  load(paste("./DataRuns/", fileList[i],sep = ""))
  
  vocabMarginalize(vocab[[4]], vocab[[5]], namesFiles[i], wordNames, 20, T)
  
}


vocabSpecificity <- function(vocabList, membList, fileName, wordNames, nWords, lowerLim = 100){
  procText <- vocabMarginalize(vocabList, membList, fileName, wordNames, nWords, F)
  
  vocabDist <- procText[[4]]
  vocabProb <- vocabDist
  membDist <<- procText[[3]]
  wordCount <- colSums(vocabDist)
  wordMarg <<- colSums(vocabDist)/sum(colSums(vocabDist))
  
  for(i in 1:nrow(vocabDist)){
    vocabProb[i,] <- vocabDist[i,]/sum(vocabDist[i,])
  }
  
  membProb <- rowSums(vocabDist)/sum(rowSums(vocabDist))
  
  for(i in 1:nrow(vocabDist)){
    for(j in 1:ncol(vocabDist)){
      vocabProb[i,j] <- exp(log(vocabProb[i,j]) + log(membProb[i])- log(wordMarg[j]))
    }
  }
  
  filter <- which(wordCount >= lowerLim)
  
  filterSet <- wordNames[filter]
  filterVocabProb <- vocabProb[,filter]
  
  
  topicNames <- matrix("", nrow= nrow(vocabProb), ncol = nWords )
  topicProbs <- matrix(0, nrow= nrow(vocabProb), ncol = nWords )
  for(i in 1:nrow(filterVocabProb)){
    wordOrder <- order(filterVocabProb[i,],decreasing = T)
    words <- filterSet[wordOrder][1:nWords]
    probs <- filterVocabProb[i,wordOrder][1:nWords]
    
    topicNames[i,] <- words
    topicProbs[i,] <- probs
  }
  
  return(list(topicNames, topicProbs))
}

specList <- array(0, dim = c(366, 22,8000))
specList[1,,] <- matrix(0, 22, 8000)
vocabSpecificityPerWord <- function(vocabList, nWords, nTopics){
  specList <- array(0, dim = c( nTopics,nWords,366))
  for(day in 1:366){
    print(day)
  vocabDist <- vocabList[[day]]
  vocabProb <- vocabDist
  wordCount <- colSums(vocabDist)
  wordMarg <<- colSums(vocabDist)/sum(colSums(vocabDist))
  
  for(i in 1:nrow(vocabDist)){
    vocabProb[i,] <- vocabDist[i,]/sum(vocabDist[i,])
  }
  
  membProb <- rowSums(vocabDist)/sum(rowSums(vocabDist))
  print(dim(vocabDist))
  for(i in 1:nrow(vocabDist)){
    for(j in 1:ncol(vocabDist)){
      vocabProb[i,j] <- exp(log(vocabProb[i,j]) + log(membProb[i])- log(wordMarg[j]))
    }
  }
  print(dim(vocabProb))
  specList[,,day] <- as.matrix(vocabProb)
  
  }
 return(specList)
}

vocabFreqPerWord <- function(vocabList, nWords, nTopics){
  specList <- array(0, dim = c( nTopics,nWords,366))
  for(day in 1:366){
    print(day)
    vocabDist <- vocabList[[day]]
   
    specList[,,day] <- as.matrix(vocabDist)
    
  }
  return(specList)
}

vocabPropPerWord <- function(vocabList, nWords, nTopics){
  specList <- array(0, dim = c( nTopics,nWords,366))
  for(day in 1:366){
    print(day)
    vocabDist <- vocabList[[day]]
    
    for(i in 1:nrow(vocabDist)){
      vocabDist[i,] <- vocabDist[i,]/sum(vocabDist[i,])
    }
    
    specList[,,day] <- as.matrix(vocabDist)
    
  }
  return(specList)
}

vocabWeightedPropPerWord <- function(vocabList,specificityWeights, nWords, nTopics){
  specList <- array(0, dim = c( nTopics,nWords,366))
  for(day in 1:366){
    print(day)
    vocabDist <- vocabList[[day]]*specificityWeights[,,day]
    
    for(i in 1:nrow(vocabDist)){
      vocabDist[i,] <- vocabDist[i,]/sum(vocabDist[i,])
    }
    
    specList[,,day] <- as.matrix(vocabDist)
    
  }
  return(specList)
}

wordNames <- as.vector(read.csv("../Duke Project/wordList.csv")[,2])
load("readout1000.Rdata")

test <- vocabSpecificityPerWord(vocab[[4]], length(wordNames), 22)
testSmooth <- vocabSpecificityPerWord(vocab[[4]], length(wordNames), 22)
test2 <- vocabFreqPerWord(vocab[[4]], length(wordNames), 22)
vocabProp <- vocabPropPerWord(vocab[[4]], length(wordNames), 22)

testSmooth[is.nan(testSmooth)] <- 0

testWeightedProportions <- vocabWeightedPropPerWord(vocab[[2]], testSmooth,length(wordNames), 22 )
write.csv(test[[1]], "HighProbabilityWords.csv")
vocabList <- vocab[[4]]
write.csv(t(test[[1]][c(19,21),]), "voc.csv")
Run16Assign = do.call("c",vocab[[1]])

write.csv(Run16Assign, "16TopicAssign.csv")

table(blockAssign)
top.Stuff<- list()
for(i in 1:16){
  
  top.Stuff[[i]] <- matrix("a", 20, 12)
}
ZimProbs <- testSmooth[9,7963,]
TrayvonProbs <- testSmooth[9,7344,]
TrayvonMartinProbs <-testSmooth[9,7345,]
LanzaProbs <- testSmooth[9,3946,]
test2[,3946,]

ZimProps <- testWeightedProportions[9,7963,]
trayvonProp <- testWeightedProportions[9,7344,]
trayvonMProp <- testWeightedProportions[9,7345,]

LanzaProbs <- testWeightedProportions[9,6219,]
SandiHookProps <- testWeightedProportions[9,3946,]


ZimProps <- vocabProp[9,7963,]
trayvonProp <- vocabProp[9,7344,]
LanzaProbs <- vocabProp[9,6219,]
SandiHookProps <- vocabProp[9,3946,]
ts.plot(cbind(ZimProps, trayvonProp,trayvonMProp))
ts.plot(cbind(LanzaProbs, SandiHookProps))
which(wordNames == "sandi.hook")
for(i in 1:12){
  for(j in 1:16){
    vectWords <- names(vocabList[[i]][j,order(vocabList[[i]][j,], decreasing = T)])[1:20]
    top.Stuff[[j]][,i] <-vectWords 
    
  }
  
  
  
  
}

for(i in 1:16){
  
  write.csv(top.Stuff[[i]],paste("Topic", i, "overYear.csv", sep = ""))
  
}
networkParams

  
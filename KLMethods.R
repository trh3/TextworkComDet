#######Number of Topic Assessor

vocabDecompKL <- function(vocab, wordDataLength, nDays, nTopics){
    mat <- vocab[[4]][[1]]
    for(i in 2:nDays){
      mat <- mat+ vocab[[4]][[1]]
    }
    
    assignVect <- do.call("c", vocab[[1]])
    d<- svd(mat)[[1]]
    m2 <- matrix(0, nrow = length(wordDataLength), ncol = nTopics)
    indices <- cbind(1:length(wordDataLength), assignVect)
    m2[indices] <- 1
    m2 <- wordDataLength %*% m2 
    m2 <- m2/sum(m2)
    
    return(kl.dist(d, m2)[[3]])
  }
  


# vocabDecompKLDaily <- function(vocab, wordDataLength,assignVect, nDays, nTopics){
#   
#   
#   
#   d<- svd(vocab)[[1]]
#   
#   m2 <- matrix(0, nrow = length(wordDataLength), ncol = nTopics)
#   
#   indices <- cbind(1:length(wordDataLength), assignVect)
#   
#   m2[indices] <- 1
#   
#   m2 <- wordDataLength %*% m2 
#   m2 <- m2/sum(m2)
#   
#   return(kl.dist(d, m2)[[3]])
#   
# }



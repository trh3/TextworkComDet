

####The Real Text Functions

library(Rcpp)
library(RcppArmadillo)
library(gtools)

sourceCpp("TopicAnalysisFunction.cpp")
sourceCpp("TopicAnalysisFunctionLongInit2.cpp")
sourceCpp("TopicAnalysisFunctionLong.cpp")
sourceCpp("TopicAnalysisFunctionTrue.cpp")
sourceCpp("DomainTopicLikelihoods.cpp")

fullTopicAnalyze <- function(full.document.matrix, target.K, alpha.prior, beta.prior,lag.window = 1 , iter = 10, control){
  data.list <- list()  
  print("Initializing Data Structures")
  for(i in 1:control[1]){
    data.list[[i]] <- full.document.matrix[which(full.document.matrix[,"date"] == i),3:dim(full.document.matrix)[[2]]]
  }
  
  
  topic.assign.list <- list()
  topic.word.matrix.list <- list()
  memb.num.list <- list()
  
  print("Fitting Day 1 Documents")
  dayone.init <- post.to.topic.initialization(data.list[[1]], target.K, alpha.prior, beta.prior, dayone.iter)
  print(sum(dayone.init[[3]]))
  topic.assign.list[[1]] <- dayone.init[[1]]
  topic.word.matrix.list[[1]] <- dayone.init[[2]]
  memb.num.list[[1]] <- dayone.init[[3]]
  daily.word.matrix.list <- list()
  
  daily.word.matrix.list[[1]] <- dayone.init[[2]]
  daily.memb.num <- list()
  daily.memb.num[[1]] <- dayone.init[[3]]
  
  print(sum(memb.num.list[[1]]))
  print("Sweeping Through Days")
  empty.mat <- matrix(0, dim(topic.word.matrix.list[[1]])[[1]], dim(topic.word.matrix.list[[1]])[[2]])
  empty.memb <- rep(0, dim(topic.word.matrix.list[[1]])[[1]])
  for(i in 2:control[1]){
    print(i)
    prev.day.memb.num <- memb.num.list[[i-1]]
    prev.day.topics <- topic.word.matrix.list[[i-1]]
    
    if(((i - lag.window) <= 0) | (lag.window == 0)) {
      day.info <- post.to.topic.longitudinal.init(data.list[[i]], prev.day.topics, empty.mat, prev.day.memb.num, empty.memb, alpha.prior, beta.prior)
      
    }else{
      day.info <- post.to.topic.longitudinal.init(data.list[[i]], prev.day.topics, daily.word.matrix.list[[i-lag.window]], prev.day.memb.num, daily.memb.num[[i-lag.window]], alpha.prior, beta.prior)
      
    }
    
    
    topic.assign.list[[i]] <- day.info[[1]]
    topic.word.matrix.list[[i]] <- day.info[[2]]
    memb.num.list[[i]] <- day.info[[3]]
    daily.word.matrix.list[[i]] <- day.info[[4]]
    daily.memb.num[[i]] <- day.info[[5]] 
  }
  
  return(list(topic.assign.list, topic.word.matrix.list, memb.num.list, daily.word.matrix.list, daily.memb.num))
  
}



post.assign.full.initialization <- function(full.document.matrix, target.K, alpha.prior, beta.prior,lag.window = 1 ,dayone.iter = 10, control){
  data.list <- list()  
  print("Initializing Data Structures")
  for(i in 1:control[1]){
    data.list[[i]] <- full.document.matrix[which(full.document.matrix[,"date"] == i),3:dim(full.document.matrix)[[2]]]
  }
  
  
  topic.assign.list <- list()
  topic.word.matrix.list <- list()
  memb.num.list <- list()
  
  print("Fitting Day 1 Documents")
  dayone.init <- post.to.topic.initialization(data.list[[1]], target.K, alpha.prior, beta.prior, dayone.iter)
  print(sum(dayone.init[[3]]))
  topic.assign.list[[1]] <- dayone.init[[1]]
  topic.word.matrix.list[[1]] <- dayone.init[[2]]
  memb.num.list[[1]] <- dayone.init[[3]]
  daily.word.matrix.list <- list()
  
  daily.word.matrix.list[[1]] <- dayone.init[[2]]
  daily.memb.num <- list()
  daily.memb.num[[1]] <- dayone.init[[3]]
  
  print(sum(memb.num.list[[1]]))
  print("Sweeping Through Days")
  empty.mat <- matrix(0, dim(topic.word.matrix.list[[1]])[[1]], dim(topic.word.matrix.list[[1]])[[2]])
  empty.memb <- rep(0, dim(topic.word.matrix.list[[1]])[[1]])
  for(i in 2:control[1]){
    print(i)
    prev.day.memb.num <- memb.num.list[[i-1]]
    prev.day.topics <- topic.word.matrix.list[[i-1]]
    
    if(((i - lag.window) <= 0) | (lag.window == 0)) {
      day.info <- post.to.topic.longitudinal.init(data.list[[i]], prev.day.topics, empty.mat, prev.day.memb.num, empty.memb, alpha.prior, beta.prior)
      
    }else{
      day.info <- post.to.topic.longitudinal.init(data.list[[i]], prev.day.topics, daily.word.matrix.list[[i-lag.window]], prev.day.memb.num, daily.memb.num[[i-lag.window]], alpha.prior, beta.prior)
      
    }

    
    topic.assign.list[[i]] <- day.info[[1]]
    topic.word.matrix.list[[i]] <- day.info[[2]]
    memb.num.list[[i]] <- day.info[[3]]
    daily.word.matrix.list[[i]] <- day.info[[4]]
    daily.memb.num[[i]] <- day.info[[5]] 
  }
  
  return(list(topic.assign.list, topic.word.matrix.list, memb.num.list, daily.word.matrix.list, daily.memb.num))
  
}

fullTopicAnalyze <- function(full.document.matrix, target.K, alpha.prior, beta.prior,lag.window = 1 , iter = 10,  nDays){
  data.list <- list()  
  print("Initializing Data Structures")
  for(i in 1:nDays){
    data.list[[i]] <- full.document.matrix[which(full.document.matrix[,"date"] == i),3:dim(full.document.matrix)[[2]]]
  }
  
  
  topic.assign.list <- list()
  topic.word.matrix.list <- list()
  memb.num.list <- list()
  print("Fitting Day 1 Documents")
  dayone.init <- post.to.topic.initialization(data.list[[1]], target.K, alpha.prior, beta.prior, iter)
  topic.assign.list[[1]] <- dayone.init[[1]]
  topic.word.matrix.list[[1]] <- dayone.init[[2]]
  memb.num.list[[1]] <- dayone.init[[3]]
  daily.word.matrix.list <- list()
  
  daily.word.matrix.list[[1]] <- dayone.init[[2]]
  daily.memb.num <- list()
  daily.memb.num[[1]] <- dayone.init[[3]]
  
  print(sum(memb.num.list[[1]]))
  print("Sweeping Through Days")
  empty.mat <- matrix(0, dim(topic.word.matrix.list[[1]])[[1]], dim(topic.word.matrix.list[[1]])[[2]])
  empty.memb <- rep(0, dim(topic.word.matrix.list[[1]])[[1]])
  for(i in 2:nDays){
    print(i)
    prev.day.memb.num <- memb.num.list[[i-1]]
    prev.day.topics <- topic.word.matrix.list[[i-1]]
    
    if(((i - lag.window) <= 0) | (lag.window == 0)) {
      day.info <- post.to.topic.longitudinal.init(data.list[[i]], prev.day.topics, empty.mat, prev.day.memb.num, empty.memb, alpha.prior, beta.prior, iter = iter)
      
    }else{
      day.info <- post.to.topic.longitudinal.init(data.list[[i]], prev.day.topics, daily.word.matrix.list[[i-lag.window]], prev.day.memb.num, daily.memb.num[[i-lag.window]], alpha.prior, beta.prior, iter = iter)
      
    }
    
    
    topic.assign.list[[i]] <- day.info[[1]]
    topic.word.matrix.list[[i]] <- day.info[[2]]
    memb.num.list[[i]] <- day.info[[3]]
    daily.word.matrix.list[[i]] <- day.info[[4]]
    daily.memb.num[[i]] <- day.info[[5]] 
  }
  
  return(list(topic.assign.list, topic.word.matrix.list, memb.num.list, daily.word.matrix.list, daily.memb.num))
  
}



top.word.vects <- function(topic.word.matrix, word.names, target.topics, target.prob = .99){
  for(k in target.topics){
    quantiles <- quantile(topic.word.matrix[k,], probs = target.prob)
    print(quantiles)
    print(paste("Topic:", k))
    print(word.names[which((topic.word.matrix[k,] >= quantiles) & (topic.word.matrix[k,] > 0)) ])
  }
}
post.to.topic.initialization <- function(dayone.document.matrix, target.K, alpha.prior, beta.prior, iter = 10){
  n.docs <- dim(dayone.document.matrix)[[1]]
  n.words <- dim(dayone.document.matrix)[[2]]
  ##Initialize Random document assignment
  topic.word.matrix <- matrix(0, target.K, n.words)
  topic.assign <- vector(length = n.docs)
  memb.num <- vector(length = target.K)
  print("Assigning Random Topics")
  for(i in 1:n.docs){
    topic.assign[i] <- sample(1:target.K, 1)
    topic.word.matrix[topic.assign[i],] <- topic.word.matrix[topic.assign[i],] + unlist(dayone.document.matrix[i,])
    memb.num[topic.assign[i]] <- memb.num[topic.assign[i]] +1
  }
  
  ##Begin the iterations:
  things <- likeli(docs = as.matrix(dayone.document.matrix), wordDist = topic.word.matrix, assign = topic.assign,membNum = memb.num,
                    priors = c(alpha.prior, beta.prior, target.K, iter))
   topic.assign <- as.vector(things[[1]])
   topic.word.matrix <-things[[2]]
   memb.num <- as.vector(things[[3]])
   
   return(list(topic.assign, topic.word.matrix, memb.num))
}


post.to.topic.longitudinal.init <- function(document.matrix, topic.word.matrix,lag.matrix, memb.num,lag.memb.num, alpha.prior, beta.prior, iter = 10){
  n.docs <- dim(document.matrix)[[1]]
  n.words <- dim(topic.word.matrix)[[2]]
  target.K <- dim(topic.word.matrix)[[1]]
  original.word.topic.matrix <- topic.word.matrix
  original.memb.num <- memb.num
  ##Initialize Random document assignment
  
  
  output <- likeliLongInit(as.matrix(document.matrix), topic.word.matrix, lag.matrix, memb.num,lag.memb.num, c(alpha.prior, beta.prior, target.K, iter))
  
  
  ##Begin the iterations:
  return(list(as.vector(output[[1]]), output[[2]], as.vector(output[[3]]), output[[4]], as.vector(output[[5]])));
}
post.to.topic.assign.dayone <- function(document.matrix, domain.likelihoods, word.dist,domain.vect, assign, membnum,alpha.prior, beta.prior, iter = 10){
 
  things <- likeliTrue(as.matrix(document.matrix), word.dist, 
                       assign,  membnum, domains = domain.vect,
                       domain.likelihoods,
                   priors = c(alpha.prior, beta.prior))
  
  return(things)
}
post.to.topic.longitudinal <- function(document.matrix, topic.word.matrix, domain.likelihoods,assign, domains ,lag.matrix, memb.num,lag.memb.num, alpha.prior, beta.prior){
  n.docs <- dim(document.matrix)[[1]]
  n.words <- dim(topic.word.matrix)[[2]]
  target.K <- dim(topic.word.matrix)[[1]]
  original.word.topic.matrix <- topic.word.matrix
  original.memb.num <- memb.num
  ##Initialize Random document assignment
  priors <- c(alpha.prior, beta.prior)
  
 
  output <- likeliLong(as.matrix(document.matrix), topic.word.matrix, lag.matrix, memb.num,lag.memb.num,assign, domains, domain.likelihoods, priors)
  
  
  ##Begin the iterations:
  return(list(as.vector(output[[1]]), output[[2]], as.vector(output[[3]]), output[[4]], as.vector(output[[5]])));
}

InitializeData <- function(full.document.matrix, control){
  
  data.list <- list();
  domain.list <- list();
  
  print("Dividing Up Data");
  for(i in 1:control[1]){
    print(i)
    data.list[[i]] <- full.document.matrix[which(full.document.matrix[,"date"] == i),3:dim(full.document.matrix)[[2]]]
    domain.list[[i]] <- full.document.matrix[which(full.document.matrix[,"date"] == i),2]
  }
  
  return(list(data.list, domain.list))
  
}
post.assign.full <- function(data.list, domain.list,word.dist.list, assign.list, membNum.list, domainLikelihoods,  alpha.prior, beta.prior,lag.window = 1 , control){

    
  print("Fitting Day 1 Documents")
  dayone.init <- post.to.topic.assign.dayone(data.list[[1]], domainLikelihoods[,,1], word.dist.list[[1]], domain.list[[1]], assign.list[[1]], membNum.list[[1]],  alpha.prior, beta.prior, dayone.iter)

  assign.list[[1]] <- dayone.init[[1]]
  word.dist.list[[1]] <- dayone.init[[2]]
  membNum.list[[1]] <- dayone.init[[3]]
  daily.word.matrix.list <- list()
  
  daily.word.matrix.list[[1]] <- dayone.init[[2]]
  daily.memb.num <- list()
  daily.memb.num[[1]] <- dayone.init[[3]]
  
  empty.mat <- matrix(0, dim(word.dist.list[[1]])[[1]], dim(word.dist.list[[1]])[[2]])
  empty.memb <- rep(0, dim(word.dist.list[[1]])[[1]])
  print("Sweeping Through Days")
  for(i in 2:control[1]){
    print(i)
    prev.day.memb.num <- membNum.list[[i-1]]
    prev.day.topics <- word.dist.list[[i-1]]
    
    if(((i - lag.window) <= 0) | (lag.window == 0)) {
      day.info <- post.to.topic.longitudinal(data.list[[i]], prev.day.topics, domainLikelihoods[,,i], assign.list[[i]], domain.list[[i]],
                                             empty.mat, prev.day.memb.num, empty.memb, alpha.prior, beta.prior)
    }else{
      day.info <- post.to.topic.longitudinal(data.list[[i]], prev.day.topics, domainLikelihoods[,,i], assign.list[[i]], domain.list[[i]], daily.word.matrix.list[[i-lag.window]],
                                                  prev.day.memb.num, daily.memb.num[[i-lag.window]], alpha.prior, beta.prior)
    }
    
    
    assign.list[[i]] <- day.info[[1]]
    word.dist.list[[i]] <- day.info[[2]]
    membNum.list[[i]] <- day.info[[3]]
    daily.word.matrix.list[[i]] <- day.info[[4]]
    daily.memb.num[[i]] <- day.info[[5]] 
  }
  
  return(list(assign.list, word.dist.list, membNum.list, daily.word.matrix.list, daily.memb.num))
  
}


####Legacy
####Text Processing Script

####Gets starting values for first days topic assignments
# post.to.topic.initializationlegacy <- function(dayone.document.matrix, target.K, alpha.prior, beta.prior, iter = 10){
#   n.docs <- dim(dayone.document.matrix)[[1]]
#   n.words <- dim(dayone.document.matrix)[[2]]
#   print(n.docs)
#   print(n.words)
#   ##Initialize Random document assignment
#   topic.word.matrix <- matrix(0, target.K, n.words)
#   print(topic.word.matrix)
#   topic.assign <- vector(length = n.docs)
#   memb.num <- vector(length = target.K)
#   for(i in 1:n.docs){
#     topic.assign[i] <- sample(1:target.K, 1)
#     topic.word.matrix[topic.assign[i],] <- topic.word.matrix[topic.assign[i],] + unlist(dayone.document.matrix[i,])
#     memb.num[topic.assign[i]] <- memb.num[topic.assign[i]] +1
#   }
#   
#   ##Begin the iterations:
#   
#   for(i in 1:iter){
#     for(j in 1:n.docs){
#       print(paste("doc: ", j))
#       curr.assign <- topic.assign[j]
#       topic.word.matrix[curr.assign,]<- topic.word.matrix[curr.assign,] - unlist(dayone.document.matrix[j,])
#       memb.num[curr.assign] <- memb.num[curr.assign]-1
#       
#       front <- (memb.num + alpha.prior)/(n.docs-1 + target.K*alpha.prior)
#       topic.prob <- vector(length = target.K) 
#       for(k in 1:target.K){
#         denom <- sum(log(sum(topic.word.matrix[k,]) +n.words*beta.prior+(1:sum(dayone.document.matrix[j,])-1)))
#         num <- 0
#         for(w in 1:n.words){
#           if(dayone.document.matrix[j,w] > 0){
#             num <- num+sum(log(topic.word.matrix[k,w]+beta.prior+1:dayone.document.matrix[j,w]- 1))
#           }
#         }
#         
#         topic.prob[k] <- num-denom
#       }
#       
#       
#       topic.prob <- log(front) + topic.prob
#       
#       topic.prob <- exp(topic.prob - max(topic.prob))
#       topic.prob <- topic.prob/sum(topic.prob)
#       topic.assign[j] <- sample(1:target.K, 1, replace = T, topic.prob)
#       
#       topic.word.matrix[topic.assign[j],] <- topic.word.matrix[topic.assign[j],] + unlist(dayone.document.matrix[j,])
#       
#       memb.num[topic.assign[j]] <- memb.num[topic.assign[j]]+1
#     }
#     
#   }
#   return(list(topic.assign, topic.word.matrix, memb.num))
# }

###Sweeps the topic assignment out to the next day
# post.to.topic.longitudinal.init <- function(document.matrix, topic.word.matrix, memb.num, alpha.prior, beta.prior){
#   n.docs <- dim(document.matrix)[[1]]
#   n.words <- dim(topic.word.matrix)[[2]]
#   target.K <- dim(topic.word.matrix)[[1]]
#   original.word.topic.matrix <- topic.word.matrix
#   original.memb.num <- memb.num
#   ##Initialize Random document assignment
#   
#   
#   topic.assign <- vector(length = n.docs)
#   for(i in 1:n.docs){
#     print(paste("Doc:", i))
#     front <- (memb.num + alpha.prior)/(n.docs-1 + target.K*alpha.prior)
#     topic.prob <- vector(length = target.K) 
#     for(k in 1:target.K){
#       denom <- sum(log(sum(topic.word.matrix[k,]) +n.words*beta.prior+(1:sum(dayone.document.matrix[i,])-1)))
#       num <- 0
#       for(w in 1:n.words){
#         if(dayone.document.matrix[i,w] > 0){
#           num <- num+sum(log(topic.word.matrix[k,w]+beta.prior+1:dayone.document.matrix[i,w]- 1))
#         }
#       }
#       
#       topic.prob[k] <- num-denom
#     }
#     
#     
#     topic.prob <- log(front) + topic.prob
#     
#     topic.prob <- exp(topic.prob - max(topic.prob))
#     topic.prob <- topic.prob/sum(topic.prob)
#     topic.assign[i] <- sample(1:target.K, 1, replace = T, topic.prob)
#     
#     topic.word.matrix[topic.assign[i],] <- topic.word.matrix[topic.assign[i],] + unlist(dayone.document.matrix[i,])
#     
#     memb.num[topic.assign[i]] <- memb.num[topic.assign[i]]+1
#   }
#   
#   ##Begin the iterations:
#   return(list(topic.assign, topic.word.matrix - original.word.topic.matrix, memb.num - original.memb.num))
# }

###Once we have starting values, this function uses the domain specific information to help clarify assignments
# post.to.topic.assign.full.info.dayone <- function(dayone.document.matrix, topic.word.matrix,topic.assign,
#                                                   domain.list, domain.interest.matrix, domain.n.posts, memb.num, alpha.prior, beta.prior, chi.prior){
#   n.docs <- dim(dayone.document.matrix)[[1]]
#   n.words <- dim(dayone.document.matrix)[[2]]
#   
#   target.K <- dim(topic.word.matrix)[[1]]
#   print(target.K)
#   for(j in 1:n.docs){
#     print(paste("doc: ", j))
#     curr.assign <- topic.assign[j]
#     topic.word.matrix[curr.assign,]<- topic.word.matrix[curr.assign,] - unlist(dayone.document.matrix[j,])
#     memb.num[curr.assign] <- memb.num[curr.assign]-1
#     
#     front <- (memb.num + alpha.prior)/(n.docs-1 + target.K*alpha.prior)* (domain.interest.matrix[,domain.list[j]]*domain.n.posts[,domain.list[j]] + chi.prior)/(domain.n.posts[domain.list[j]] - 1 + target.K*chi.prior)
#     topic.prob <- vector(length = target.K) 
#     for(k in 1:target.K){
#       denom <- sum(log(sum(topic.word.matrix[k,]) +n.words*beta.prior+(1:sum(dayone.document.matrix[j,])-1)))
#       num <- 0
#       for(w in 1:n.words){
#         if(dayone.document.matrix[j,w] > 0){
#           num <- num+sum(log(topic.word.matrix[k,w]+beta.prior+1:dayone.document.matrix[j,w]- 1))
#         }
#       }
#       
#       topic.prob[k] <- num-denom
#     }
#     
#     
#     topic.prob <- log(front) + topic.prob
#     
#     topic.prob <- exp(topic.prob - max(topic.prob))
#     topic.prob <- topic.prob/sum(topic.prob)
#     topic.assign[j] <- sample(1:target.K, 1, replace = T, topic.prob)
#     
#     topic.word.matrix[topic.assign[j],] <- topic.word.matrix[topic.assign[j],] + unlist(dayone.document.matrix[j,])
#     
#     memb.num[topic.assign[j]] <- memb.num[topic.assign[j]]+1
#   }
#   
#   
#   return(list(topic.assign, topic.word.matrix, memb.num))
#   
#   
# }
###LEGACY
# post.to.topic.assign.full.info<- function(document.matrix, topic.word.matrix,
#                                           domain.list, domain.interest.matrix, domain.n.posts, memb.num, alpha.prior, beta.prior, chi.prior){
#   n.docs <- dim(dayone.document.matrix)[[1]]
#   n.words <- dim(word.dist)[[2]]
#   target.K <- dim(topic.word.matrix)[[1]]
#   original.word.dist <- topic.word.matrix
#   original.memb.num <- memb.num
#   
#   topic.assign <- vector(length = n.docs)
#   for(j in 1:n.docs){
#     print(paste("doc: ", j))
#     
#     
#     front <- (memb.num + alpha.prior)/(n.docs-1 + target.K*alpha.prior)* (domain.interest.matrix[,domain.list[j]]*domain.n.posts[,domain.list[j]] + chi.prior)/(domain.n.posts[domain.list[j]] - 1 + target.K*chi.prior)
#     topic.prob <- vector(length = target.K) 
#     for(k in 1:target.K){
#       denom <- sum(log(sum(topic.word.matrix[k,]) +n.words*beta.prior+(1:sum(document.matrix[j,])-1)))
#       num <- 0
#       for(w in 1:n.words){
#         if(dayone.document.matrix[j,w] > 0){
#           num <- num+sum(log(topic.word.matrix[k,w]+beta.prior+1:document.matrix[j,w]- 1))
#         }
#       }
#       
#       topic.prob[k] <- num-denom
#     }
#     
#     
#     topic.prob <- log(front) + topic.prob
#     
#     topic.prob <- exp(topic.prob - max(topic.prob))
#     topic.prob <- topic.prob/sum(topic.prob)
#     topic.assign[j] <- sample(1:target.K, 1, replace = T, topic.prob)
#     
#     topic.word.matrix[topic.assign[j],] <- topic.word.matrix[topic.assign[j],] + unlist(document.matrix[j,])
#     
#     memb.num[topic.assign[j]] <- memb.num[topic.assign[j]]+1
#   }
#   
#   return(list(topic.assign, topic.word.matrix - original.topic.word.matrix, memb.num - original.memb.num))
#   
#   
# }



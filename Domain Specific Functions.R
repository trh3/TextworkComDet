#####Domain Interest and PostRate Functions

####Given Activation and Given Topic Assignment and Topic Interests, return topic specific activity, and interest parameters
sourceCpp("PostRegressionProcess.cpp")


###Initial Block Setup Process:


PostRegressionProcess <- function(assignmentVector, domainVector, dayVector, domainInterestMatrix, activationMatrix, 
                                    ndays, ndomains, ntopics){
  
  control = c(ndays, ndomains,ntopics)
  
  dataset = postRegressionData(assignmentVector, domainVector, dayVector, domainInterestMatrix, activationMatrix, control)
  
  return(dataset);
};

postDataLikelihood <- function(data, topic.excite.parameters, domain.post.rates){
  
  rate <- domain.post.rates[data[,4]]*data[,2] + data[,3]*topic.excite.parameters
  return(sum(dpois(data[,1], rate, log = T)))
  
  
}
post.to.topic.regression <- function(data, n.topics, topic.excite.params, domain.post.rates,
                                     regparampriorsmean, regparampriorscale, prop.scale){

  top.exc.para <- topic.excite.params


  print("Sampling Post To Topic Regression Parameters")
  
  for(i in 1:n.topics)  {
    topic.data <- data[which(data[,6] == i),]
    if(dim(topic.data)[[1]] > 0){
      
    prop.top.exc.para <- rtruncnorm(1,a=0, mean = top.exc.para[i], sd =  prop.scale)
    
    curr.rate <- postDataLikelihood(topic.data, top.exc.para[i], domain.post.rates)
    prop.rate.exc <- postDataLikelihood(topic.data, prop.top.exc.para, domain.post.rates)
      curr.exc.prior <- log(dtruncnorm(top.exc.para[i], a = 0,  mean = regparampriorsmean, sd = regparampriorscale))
    prop.exc.prior <- log(dtruncnorm(prop.top.exc.para, a= 0, mean = regparampriorsmean, sd = regparampriorscale))

    exc.ratio <- (log(dtruncnorm(top.exc.para[i], a= 0,  mean = prop.top.exc.para, sd = prop.scale))-
                    log(dtruncnorm(prop.top.exc.para, a =0,  mean = top.exc.para[i], sd = prop.scale)))

    
    exc.a <-  (prop.rate.exc - curr.rate) + (prop.exc.prior - curr.exc.prior) + exc.ratio
    
    if(log(runif(1,0,1)) < exc.a){ top.exc.para[i] <- prop.top.exc.para}


  }
  }
  return(top.exc.para)
  
}

domain.base.rate.regression <- function(data, n.domains, topic.excite.params, domain.post.rates,
                                     regparampriorsmean, regparampriorscale, prop.scale){
  require(truncnorm)
  top.exc.para <- topic.excite.params

  
  print("Sampling Domain Base Rates")  
  for(i in 1:n.domains){

    topic.data <- data[which(data[,4] == i),]
   
    domain.post.ratesProp <- rtruncnorm(1, a = 0,mean = domain.post.rates[i], sd =  prop.scale)
 

    curr.rate <- sum(dpois(topic.data[,1],(topic.data[,2]*domain.post.rates[i] +topic.data[,3]*top.exc.para[topic.data[,6]]), log = T), na.rm = T)
    prop.rate.exc <-sum(dpois(topic.data[,1], (topic.data[,2]*domain.post.ratesProp +topic.data[,3]*top.exc.para[topic.data[,6]]),log = T), na.rm = T)
    curr.exc.prior <- log(dtruncnorm(domain.post.rates[i], a = 0, mean = regparampriorsmean, sd = regparampriorscale))
    prop.exc.prior <- log(dtruncnorm(domain.post.ratesProp, a = 0, mean = regparampriorsmean, sd = regparampriorscale))
    
    exc.ratio <- (log(dtruncnorm(domain.post.rates[i], a = 0,  mean = domain.post.ratesProp, sd = prop.scale))-
                    log(dtruncnorm(domain.post.ratesProp, a = 0,  mean = domain.post.rates[i], sd = prop.scale)))
    
    
    exc.a <-  (prop.rate.exc - curr.rate) + (prop.exc.prior - curr.exc.prior) + exc.ratio
    
    if(log(runif(1,0,1)) < exc.a){domain.post.rates[i] <- domain.post.ratesProp}
    
    ###Topic Interest Parameter
    
  }
  return(domain.post.rates)
  
}

block.domain.eval <- function(domain.topic.interest, block.topic.interest){

likelihood <- log(ddirichlet(domain.topic.interest, block.topic.interest))

}

domain.specific.topic.interest.sampler <- function(processed.data, domain.topic.interest.matrix, topic.excite.params,
                                                   domain.post.rates, chi.prior, domain.interest.prior, block.interest.matrix, block.membership,blockIntConcentration,
                                                  n.domains){
  
    
  domain.topic.interest.matrix.copy <- domain.topic.interest.matrix
  print("Sampling Domain Topic Interests")
  for(i in 1:n.domains){

    domain.subset <- processed.data[which(processed.data[,4] == i),]
    possible.topics <- unique(domain.subset[,6])
    
    total.posts <- sum(domain.subset[,1])
    toggle = T
    
    prop.topic.interest.vector <- rdirichlet(1,domain.topic.interest.matrix[,i]*total.posts + chi.prior)
    prop.likelihood <- sum(dpois(x = domain.subset[,1], lambda = (prop.topic.interest.vector[domain.subset[,6]]*domain.post.rates[i] + domain.subset[,7]*prop.topic.interest.vector[domain.subset[,6]]*topic.excite.params[domain.subset[,6]]), log = T), na.rm =T) + 
	  block.domain.eval(prop.topic.interest.vector, block.interest.matrix[which(block.interest.matrix[,1] == block.membership[i]),-1]*blockIntConcentration+1)
  
    curr.likelihood <-sum(dpois(x = domain.subset[,1], lambda = (domain.subset[,2]*domain.post.rates[i] + domain.subset[,3]*topic.excite.params[domain.subset[,6]]), log = T),na.rm =T) + 
	  block.domain.eval(domain.topic.interest.matrix[,i], block.interest.matrix[which(block.interest.matrix[,1] == block.membership[i]),-1]*blockIntConcentration+1)
  
    
    ratio <- log(ddirichlet(domain.topic.interest.matrix[,i],prop.topic.interest.vector*total.posts + chi.prior ))- log(ddirichlet(prop.topic.interest.vector, domain.topic.interest.matrix[,i]*total.posts + chi.prior))

    a <- (prop.likelihood - curr.likelihood)  + ratio
    if(log(runif(1,0,1)) < a){
      domain.topic.interest.matrix.copy[,i] <- prop.topic.interest.vector
      
    }
    
    
  }
  return(domain.topic.interest.matrix.copy)
  
}

block.sums.calculator <- function(domain.topic.interest.matrix, block.memb,n.blocks, n.topics){
  block.average <- matrix(0, n.topics, n.blocks)
  for(i in 1:n.blocks){
    if(is.null(dim(domain.topic.interest.matrix[,which(block.memb == i)]))){
      block.average[,i] <- domain.topic.interest.matrix[,which(block.memb == i)]
    }else{
      block.average[,i] <-rowSums(domain.topic.interest.matrix[,which(block.memb == i)])
    }
  }
  return(block.average)
}

blockLikelihoodCalculator <- function(domainTopicInterest, blockSums,  blockMemb, nBlocks, nTopics, nDomains){
  require(gtools)
  domainBlockLikelihoods <- matrix(0, nDomains, nBlocks)
  
  for(i in 1:length(blockMemb)){
    for(j in 1:nBlocks){
      if(j == blockMemb[i]){
        domainBlockLikelihoods[i,j] = log(ddirichlet(domainTopicInterest[,i], alpha = blockSums[,j]))
      }else{
        domainBlockLikelihoods[i,j] = log(ddirichlet(domainTopicInterest[,i], alpha = blockSums[,j]))
        
      }
      
    }
  }
  return(domainBlockLikelihoods)  
}

SameElements <- function(a, b) return(identical(sort(a), sort(b)))

blockLabeller <- function(topicInterestList){
  
  topicInterestList <- lapply(X = topicInterestList, FUN = sort)
  
  uniqueBlocks <- unique(topicInterestList) 
  
  blockMemb <- vector(length = length(topicInterestList))
  membNum <- vector(length = length(uniqueBlocks))  
  keys <- lapply(uniqueBlocks, interestHash)
  keys <- unlist(keys)
  for(i in 1:length(topicInterestList)){
    
    for(j in 1:length(uniqueBlocks)){
      
      if(identical(topicInterestList[[i]], uniqueBlocks[[j]])){
        break
      }

    }
    blockMemb[i] = keys[j]
    membNum[j] = membNum[j] + 1
  }
  membNum <- data.table(keys = keys, membNum = membNum)
  setkey(membNum, "keys")
return(list(topicInterestList, uniqueBlocks, blockMemb, membNum))    
}

BlockGenerator  <- function(nTopics){
  comb3 <- combinations(n = nTopics, r = 3)
  comb2 <- combinations(n = nTopics, r = 2)
  comb1 <- combinations(n = nTopics, r = 1)
  comb3 <- split(comb3, row(comb3))
  comb2 <- split(comb2, row(comb2))
  comb1 <- split(comb1, row(comb1))
  allcomb <- c(comb3, comb2, comb1)
  allcomb[[length(allcomb) +1]] <- c(1:nTopics)
  
  blockMat <- matrix(0, nrow = length(allcomb), ncol = nTopics+1)
  for( i in 1: length(allcomb)){
    blockMat[i, allcomb[[i]] + 1] = 1
  }

  keys <- unlist(lapply(allcomb, interestHash))
  keys <- as.numeric(keys)
  blockMat[,1] = keys  
  membNum <- data.table(keys = keys, membNum = rep(0, length(keys)))
  setkey(membNum, keys)
  return(list(blockMat, membNum))
}

initialMembershipAlloc <- function(blockMembNumbers, blockMembVector){
  
  for(i in blockMembVector){
    blockMembNumbers[list(i), membNum:=blockMembNumbers[list(i), membNum] + 1]
  }
  return(blockMembNumbers)
}

  


processedDataToDomainTopicCount <- function(processedData){
  data<-data.table(processedData)
  dataSum <- data[,sum(V1), by = c("V4", "V6")]
  setkey(dataSum, V4)
  return(dataSum)
}


blockLikelihoodAssign <- function(domainTopicCounts, domainTopicInterests, blockMembNumbers, blockMembVector,
                                  fullLinkData, currentNetworkParams, controlList, topicInterestDyads, priors){
  nBlocksPrior = priors[["nBlocksPrior"]]
  nMembPrior = priors[["nMembPrior"]]
  for(i in 1:controlList[["nDomains"]]){
  print(paste("Assigning Block for Domain", i))
   currDomain  <- domainTopicCounts[list(i)]
   topicsInterest <- which(currDomain[,V1 > 0])
   if(length(topicsInterest) >= 3){
   
     comb3 <- combinations(n = length(topicsInterest) ,v = topicsInterest,  r = 3)
     comb2 <- combinations(n = length(topicsInterest) ,v = topicsInterest,  r = 2)
     comb1 <- combinations(n = length(topicsInterest) ,v = topicsInterest,  r = 1)
     
     comb3 <- split(comb3, row(comb3))
     comb2 <- split(comb2, row(comb2))
     comb1 <- split(comb1, row(comb1))
     
     allcomb <- c(comb3, comb2, comb1)
   }else{
     if(length(topicsInterest) == 2){
     comb2 <- combinations(n = length(topicsInterest) ,v = topicsInterest,  r = 2)
     comb1 <- combinations(n = length(topicsInterest) ,v = topicsInterest,  r = 1)
     
     comb2 <- list(comb2)
     comb1 <- split(comb1, row(comb1))
     
     allcomb <- c(comb2, comb1)
     }else{
       comb1 <- combinations(n = length(topicsInterest) ,v = topicsInterest,  r = 1)
       comb1 <- list(comb1)
       allcomb <- comb1
     }
   }
   blockMembNumbers[list(blockMembVector[i]), membNum := blockMembNumbers[list(blockMembVector[i]), membNum] - 1]
   currBlocks <- sum(blockMembNumbers[,membNum > 0])
   print(currBlocks)
   allcomb[[length(allcomb)+1]]= 1:nTopics
   ddirch <- list()
   keys <- unlist(lapply(allcomb, interestHash))
   for(j in 1:length(allcomb)){
     comp <- rep(1, nTopics)
     comp[allcomb[[j]]] <- 100
     ddirch[[keys[j]]] <- log(ddirichlet(domainTopicInterests[,i], comp))
     if(blockMembNumbers[list(keys[[j]]),membNum] == 0){
       ddirch[[keys[j]]] <- ddirch[[keys[j]]]  + dpois(currBlocks + 1, lambda = nBlocksPrior, log = T) + 
         log(nMembPrior/(nMembPrior*(currBlocks + 1)+ sum(blockMembNumbers[,membNum])))  
         
     }else{
       
       ddirch[[keys[j]]] <- ddirch[[keys[j]]] +dpois(currBlocks, lambda = nBlocksPrior, log = T) +
         log((nMembPrior + blockMembNumbers[list(keys[[j]]), membNum])/(nMembPrior*(currBlocks +1)+ sum(blockMembNumbers[,membNum])))
           
         }
       }
   #Adjustment for the all-topic block
     ddirch[[keys[length(keys)]]] = ddirch[[keys[length(keys)]]] + log(1/currBlocks)
   infremove <-  !(unlist(ddirch) == -Inf)
   if(length(which(infremove)) == 0){
     print("Removed All likelihoods, assigning to first")
     
     ddirch <- -1
     keys <- keys[1]
   }else{
     ddirch <- unlist(ddirch)[infremove]
     keys <- keys[infremove]   
   }
#     print(as.vector(unlist(ddirch)))
#     print(as.numeric(keys))
#     print(typeof(fullLinkData))
#     print(typeof(currentNetworkParams))
#     print(typeof(blockMembVector))
#     print(typeof(i))
   output  <-  BlockAssignment(as.matrix(fullLinkData), currentNetworkParams, as.vector(unlist(ddirch)),as.numeric(keys), topicInterestDyads, 
                   blockMembVector, domain = i)
   blockMembVector[i] <- output[[2]]
   blockMembNumbers[list(blockMembVector[i]), membNum:=blockMembNumbers[list(blockMembVector[i]), membNum] + 1] 
   fullLinkData <- output[[1]]
   }


  
  return(list(blockMembVector, fullLinkData))
}


interestHash <- function(topicInterestVector){
  
  return(paste(topicInterestVector, sep = "0",collapse = "0"))
  
}

###Not Used
block.concentration.sampler <- function(block.concentration.vect,domain.topic.interest.matrix, block.average.matrix, block.memb,prop.scale,prior.concentration, 
                                        prior.scale, n.blocks, n.domains){
  
  require(truncnorm)
  for(i in 1:n.blocks){
    prop.concentration <- rtruncnorm(1,a = 0, mean = block.concentration[i], sd = prop.scale)
    
    
    curr.likelihood <- log(ddirichlet(t(domain.topic.interest.matrix[,which(block.memb == i)]), block.average.matrix[,i]*block.concentration.vect[i]))
    prop.likelihood <- log(ddirichlet(t(domain.topic.interest.matrix[,which(block.memb == i)]), block.average.matrix[,i]*prop.concentration))
    toggle = 0
    if(any(curr.likelihood > 0)){
      curr.likelihood = 0 
    }
    if(any(prop.likelihood > 0)){
      prop.likelihood = 0
    }
    
    
    curr.prior <- dgamma(block.concentration.vect[i], shape = prior.concentration/prior.scale, scale = prior.scale, log = T)
    prop.prior <- dgamma(prop.concentration, shape = prior.concentration/prior.scale, scale = prior.scale, log = T)
    
    ratio <- log(dtruncnorm(block.concentration[i], a = 0, mean = prop.concentration, sd = prop.scale)) -
      log(dtruncnorm(prop.concentration, a = 0, mean = block.concentration[i], sd = prop.scale))
    
    a <- prop.likelihood + prop.prior - (curr.likelihood - curr.prior) + ratio
    
    if(log(runif(1,0,1)) > a){
      block.concentration[i] <- prop.concentration
      
    }
  }
  return(block.concentration)
}

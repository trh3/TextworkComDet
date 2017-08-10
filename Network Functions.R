###Linkage Information
###Updated documentation 8/10/2017

sourceCpp("LinkDataProcess.cpp")
sourceCpp("DomainInterestDyad.cpp")

##This function preprocesses linkage information into the correct format
##Parameters are:
##Linkcube: an N by N by T array with 0/1 entries. A 1 in the [i,j,t] represents that individual i linked to individual j on day t.
##lagMat: This is a N by T matrix with 0/1 entries. A 1 in the [i,t] represents that individual i has been linked to within the last 
##X number of days (User determines this)
##covarMat: This is a N by N matrix with numeric entries that represents the block/interest relationship between each individual.
##This is primarily used within the sampler, for starting values for your preprocesser
##one should randomly allocate individuals to K communities, then assign a 0/1 for same community in this matrix.
##nDays: Number of time points.
linkDataForRegression <- function(linkCube, lagMat, covarMat, nDays){
  require(data.table)
  require(statnet)
  data1 <- list()
  data2 <- list()
  data3 <-  list()
  for(i in 1:nDays){
    print(paste("Compiling Day", i))
    net <-as.network(linkCube[,,i])
    net%v%"lag"<- as.numeric(lagMat[,i] > 0)
    dataSet <-  ergmMPLE(net ~ edges() + nodeicov("lag") + edgecov(covarMat) + nodeocov("vertex.names")+ nodeicov("vertex.names"));
    data1[[i]] <- dataSet[[1]]
    data2[[i]] <- dataSet[[2]]
    data3[[i]] <- dataSet[[3]]
  }
  rm(dataset)
  rm(net)
  gc()
  reworkedLinks <<- as.data.table(as.matrix(cbind(do.call("c", data1),do.call("rbind", data2),do.call("c", data3))))
  reworkedLinks[, V7:= sum(V7), by = list(V1, edges, nodeicov.lag, edgecov.covarMat, nodeocov.vertex.names, nodeicov.vertex.names)]
  setkey(reworkedLinks)
  testTable <- unique(reworkedLinks)
  rm(reworkedLinks)
  gc()
  return(as.matrix(testTable))
  
}
linkRegressionLikelihood <- function(fullLinkData, currentParams){
  
  match.prop.likelihood <- log((exp((fullLinkData[,2]*currentParams[1]+fullLinkData[,3]*currentParams[2]+fullLinkData[,4]*currentParams[3]+fullLinkData[,8]*currentParams[4]+fullLinkData[,9]*currentParams[5]))))*(fullLinkData[,1])*fullLinkData[,7]
  match.prop.normalize <- log((1+exp((fullLinkData[,2]*currentParams[1]+fullLinkData[,3]*currentParams[2]+fullLinkData[,4]*currentParams[3]+fullLinkData[,8]*currentParams[4]+fullLinkData[,9]*currentParams[5]))))*fullLinkData[,7]
  prop.likelihood <- match.prop.likelihood-match.prop.normalize
  return(sum(prop.likelihood, na.rm = T))
}
linkRegression <- function(fullLinkData, currentParams, proposalSds, priorMeans, priorSd){
  
  propInt <- rnorm(1,currentParams[1], proposalSds[1])
  proplag <- rnorm(1,currentParams[2], proposalSds[2])
  propBlock <- rnorm(1,currentParams[3], proposalSds[3])
  print(currentParams)
  propOut <- rnorm(1,currentParams[4], proposalSds[4])
  propIn <- rnorm(1,currentParams[5], proposalSds[5])
  ##IntSampler
  
  propParams = currentParams
  propParams[1] = propInt
  
  currLik <- linkRegressionLikelihood(fullLinkData,currentParams)
  propLik <- linkRegressionLikelihood(fullLinkData,propParams)
  
  priorCurr <- dnorm(currentParams[1], priorMeans[1], priorSd[1])
  priorProp <- dnorm(propParams[1], priorMeans[1], priorSd[1])
  
  a <- (propLik + priorProp) - (currLik + priorCurr)
  
  if(log(runif(1,0,1)) < a){
    
    currentParams = propParams
  }
  
  print(currLik)
  
  
  propParams = currentParams
  propParams[2] = proplag
  
  currLik <- linkRegressionLikelihood(fullLinkData,currentParams)
  propLik <- linkRegressionLikelihood(fullLinkData,propParams)
  
  priorCurr <- dnorm(currentParams[2], priorMeans[2], priorSd[2])
  priorProp <- dnorm(propParams[2], priorMeans[2], priorSd[2])
  
  a <- (propLik + priorProp) - (currLik + priorCurr)
  
  if(log(runif(1,0,1)) < a){
    
    currentParams = propParams
  }
  ##Block
  propParams = currentParams
  propParams[3] = propBlock
  currLik <- linkRegressionLikelihood(fullLinkData,currentParams)
  propLik <- linkRegressionLikelihood(fullLinkData,propParams)
  
  priorCurr <- dnorm(currentParams[3], priorMeans[3], priorSd[3])
  priorProp <- dnorm(propParams[3], priorMeans[3], priorSd[3])
  
  a <- (propLik + priorProp) - (currLik + priorCurr)
  
  if(log(runif(1,0,1)) < a){
    
    currentParams = propParams
  }

  ##Outdegree
  propParams = currentParams
  propParams[4] = propOut
  currLik <- linkRegressionLikelihood(fullLinkData,currentParams)
  propLik <- linkRegressionLikelihood(fullLinkData,propParams)
  
  priorCurr <- dnorm(currentParams[4], priorMeans[4], priorSd[4])
  priorProp <- dnorm(propParams[4], priorMeans[4], priorSd[4])
  
  a <- (propLik + priorProp) - (currLik + priorCurr)
  
  if(log(runif(1,0,1)) < a){
    
    currentParams = propParams
  }
  
  ##Outdegree
  propParams = currentParams
  propParams[5] = propIn
  currLik <- linkRegressionLikelihood(fullLinkData,currentParams)
  propLik <- linkRegressionLikelihood(fullLinkData,propParams)
  
  priorCurr <- dnorm(currentParams[5], priorMeans[5], priorSd[5])
  priorProp <- dnorm(propParams[5], priorMeans[5], priorSd[5])
  
  a <- (propLik + priorProp) - (currLik + priorCurr)
  
  if(log(runif(1,0,1)) < a){
    
    currentParams = propParams
  }
  
  
  
  
  return(currentParams)
}




####Legacy
# linkage.list.to.data.frame <- function(linkage.list, block.memb, domain.interest.matrix, n.domains, n.days, lag = 1){
#   library(data.table)
#   total.size <- (n.domains*n.domains-n.domains)*n.days
#   link.data <- data.table(edges = rep(0, total.size), sender.domain = rep(0,total.size), receiver.domain = rep(0,total.size),
#                           block.match = rep(0, total.size), topic.interest = rep(0, total.size), lag = rep(0, total.size), day =  rep(0, total.size) )
#   
#   count = 0
#   for(i in 1:n.days){
#     for(j in 1:n.domains){
#       toggle = 0
#       if(i != 1){
#         if(i-lag <= 0){end.point = 1}else{end.point = i-lag}
#         
#         for(p in end.point:(i-1)){
#           if(any(linkage.list[[p]][-j,j] > 0)){
#             toggle = 1
#           }
#         }}
#       for(k in (1:n.domains)[-j]){
#         count = count + 1
#         if(linkage.list[[i]][j,k] > 0){set(link.data, as.integer(count), 1L, 1)}else{set(link.data, as.integer(count), 1L, 0)}
#         set(link.data, as.integer(count), 2L, j)
#         set(link.data, as.integer(count), 3L, k)
#         if(block.memb[j] == block.memb[k]){
#           
#           set(link.data, as.integer(count), 4L, 1)}else{
#             
#             set(link.data, as.integer(count), 4L, 0)}
#         set(link.data, as.integer(count), 5L, sum(domain.interest.matrix[,j] * domain.interest.matrix[,k]))
#         set(link.data, as.integer(count), 6L, toggle)
#         set(link.data, as.integer(count), 7L, i)
#         
#         
#       }
#     }
#   }
#   setkey(link.data, sender.domain, receiver.domain, lag)
#   link.data[,state:=edges]
#   sum.link.data <- link.data[,sum(edges),by = .(sender.domain, receiver.domain, lag, block.match, topic.interest)]
#   sum.link.data[,V2 := (n.days - V1)]
#   sum.link.data[, topic.interest:= sum.link.data[,topic.interest]*(1-sum.link.data[,block.match])]
#   setkey(sum.link.data, sender.domain, receiver.domain)
#   return(sum.link.data)
#   
# }
# 
# link.regression.likelihood <- function(link.data.table, intercept, block.match.p, topic.int, lag){
#   match.prop.likelihood <- exp((intercept +link.data.table[,block.match]*block.match.p + link.data.table[,topic.interest]*topic.int+
#                                   link.data.table[,lag]*lag))^link.data.table[,V1]
#   
#   match.prop.normalize <- (1+ exp(intercept+ link.data.table[,block.match]*block.match.p + link.data.table[,topic.interest]*topic.int+
#                                     link.data.table[,lag]*lag))^(link.data.table[,V1]+link.data.table[,V2])
#   prop.likelihood <- log(match.prop.likelihood/match.prop.normalize)
#   #   print(prop.likelihood)
#   #   print(sum(prop.likelihood, na.rm = T))
#   return(sum(prop.likelihood, na.rm = T))
#   
# }
# 
# link.regression.sampler <- function(link.data.table, int, block.match.param, topic.interest.param, lag.param,
#                                     reg.prior.mean.vect,reg.prior.sd.vect, prop.sd.vect){
#   
#   match.prop <- rnorm(1, mean = block.match.param, sd = prop.sd.vect[1])
#   interest.prop <- rnorm(1, mean = topic.interest.param, sd = prop.sd.vect[2])
#   lag.prop <- rnorm(1, mean = lag.param, sd = prop.sd.vect[3])
#   int.prop <- rnorm(1, int, prop.sd.vect[4])
#   
#   #MH step for block.match
#   prop.likelihood <- link.regression.likelihood(link.data.table,int , match.prop, topic.interest.param, lag.param)
#   curr.likelihood <- link.regression.likelihood(link.data.table,int,  block.match.param, topic.interest.param, lag.param)
#   
#   prop.prior <- dnorm(match.prop, reg.prior.mean.vect[1], reg.prior.sd.vect[1], log = T)
#   curr.prior <- dnorm(block.match.param, reg.prior.mean.vect[1], reg.prior.sd.vect[1], log = T)
#   
#   ratio <- dnorm(block.match.param, match.prop, prop.sd.vect[1], log = T)- dnorm(match.prop, block.match.param, prop.sd.vect[1], log = T)
#   
#   a <- prop.likelihood-curr.likelihood + prop.prior- curr.prior + ratio
#   
#   if(log(runif(1,0,1)) < a){
#     block.match.param <- match.prop
#   }
#   
#   
#   #MH stop for interest.param
#   prop.likelihood <- link.regression.likelihood(link.data.table, int, block.match.param, interest.prop, lag.param)
#   curr.likelihood <- link.regression.likelihood(link.data.table, int, block.match.param, topic.interest.param, lag.param)
#   
#   prop.prior <- dnorm(interest.prop, reg.prior.mean.vect[2], reg.prior.sd.vect[2], log = T)
#   curr.prior <- dnorm(topic.interest.param, reg.prior.mean.vect[2], reg.prior.sd.vect[2], log = T)
#   
#   ratio <- dnorm(topic.interest.param, interest.prop, prop.sd.vect[2], log = T)- dnorm(interest.prop, topic.interest.param, prop.sd.vect[2], log = T)
#   
#   a <- prop.likelihood-curr.likelihood + prop.prior- curr.prior + ratio
#   
#   if(log(runif(1,0,1)) < a){
#     topic.interest.param <- interest.prop
#   }
#   
#   #MH stop for lag
#   prop.likelihood <- link.regression.likelihood(link.data.table,int, block.match.param, topic.interest.param, lag.prop)
#   curr.likelihood <- link.regression.likelihood(link.data.table,int, block.match.param, topic.interest.param, lag.param)
#   
#   prop.prior <- dnorm(lag.prop, reg.prior.mean.vect[3], reg.prior.sd.vect[3], log = T)
#   curr.prior <- dnorm(topic.interest.param, reg.prior.mean.vect[3], reg.prior.sd.vect[3], log = T)
#   
#   ratio <- dnorm(lag.param, lag.prop, prop.sd.vect[3], log = T)- dnorm(lag.prop, lag.param, prop.sd.vect[3], log = T)
#   
#   a <- prop.likelihood-curr.likelihood + prop.prior- curr.prior + ratio
#   
#   if(log(runif(1,0,1)) < a){
#     lag.param <- lag.prop
#   }
#   
#   #MH stop for int
#   prop.likelihood <- link.regression.likelihood(link.data.table,int.prop, block.match.param, topic.interest.param, lag.param)
#   curr.likelihood <- link.regression.likelihood(link.data.table,int, block.match.param, topic.interest.param, lag.param)
#   
#   prop.prior <- dnorm(lag.prop, reg.prior.mean.vect[4], reg.prior.sd.vect[4], log = T)
#   curr.prior <- dnorm(topic.interest.param, reg.prior.mean.vect[3], reg.prior.sd.vect[3], log = T)
#   
#   ratio <- dnorm(lag.param, lag.prop, prop.sd.vect[4], log = T)- dnorm(lag.prop, lag.param, prop.sd.vect[4], log = T)
#   
#   a <- prop.likelihood-curr.likelihood + prop.prior- curr.prior + ratio
#   
#   if(log(runif(1,0,1)) < a){
#     int <- int.prop
#   }
#   
#   
#   return(c(int,block.match.param, topic.interest.param, lag.param))
# }





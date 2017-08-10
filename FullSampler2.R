  ######Full Sampler Script
library(data.table)
library(gtools)
library(statnet)
library(truncnorm)
library(Rcpp)
library(RcppArmadillo)



source("Text Functions.R")
source("Domain Specific Functions.R")
source("DailyActivationScripts.R")
source("Network Functions.R")
source("BlockAssignmentFunctions.R") 



topicNetworkSampler <- function(wordData, processedLinkData, nTopics, nBlocks, control, iterList, priorList, miscList, runName){
#######Begin By Initializing The Vocab
  dir.create(paste("./", runName,"/", sep = ""))
  print("Initializing Some Starting Values")
  domainVector <<- wordData[,2]
  dayVector <<- wordData[,1]
  docLength <<- rowSums(wordData[,-1:-2])
  initialDomainInterest <<- t(rdirichlet(control[["nDomains"]], rep(1, nTopics)));
  initialActivation <<- matrix(rbinom(control[["nDays"]]*nTopics,1, .2), control[["nDays"]], nTopics)
  
  initialExciteParams <<- rep(.25, nTopics)
  initialDomainRates <<- rep(4, control[["nDomains"]])
  blockMembOutput <<- BlockGenerator(nTopics)   
  blockNames <<- blockMembOutput[[2]][,keys]
  blockMembNumbers <<- blockMembOutput[[2]]
  blockMat <<- blockMembOutput[[1]]
  initialBlockAssign <<- sample(blockNames, control[["nDomains"]], T)
  blockMembNumbers <<- initialMembershipAlloc(blockMembNumbers, initialBlockAssign)
  networkParams <<- c(-7,1,1,1,1)
  print("Initializing Topic Assignment")
  vocab <<- post.assign.full.initialization(wordData, nTopics, 
                                                   priorList[["alpha"]], priorList[["beta"]], 
                                                   miscList[["lag.window"]], iterList[["dayone.iter"]], c(control[["nDays"]], control[["nDomains"]], nTopics)) 
   


  initialAssign <<- do.call("c", vocab[[1]])
  assign <- initialAssign
  
  initialProcessedData <<- PostRegressionProcess(initialAssign, domainVector,
                                                dayVector, initialDomainInterest,
                                                initialActivation, control[["nDays"]], 
                                                control[["nDomains"]], nTopics)[[1]];

    postRegParams<<- post.to.topic.regression(initialProcessedData, nTopics, initialExciteParams,
                             initialDomainRates, priorList[["PostRegressionMean"]], priorList[["PostRegressionSD"]],
                             miscList[["PostRegressionProp"]])
   
   exciteParams <<- postRegParams

   domainRates <<- domain.base.rate.regression(initialProcessedData, control[["nDomains"]], exciteParams,
                                initialDomainRates, priorList[["DomainRateMean"]], priorList[["DomainRateSD"]], miscList[["DomainRateProp"]])
#   
    domainInterest <<- domain.specific.topic.interest.sampler(initialProcessedData, initialDomainInterest, exciteParams, domainRates, 
                                                            priorList[["chi"]], priorList[["domainInterest"]], blockMat, initialBlockAssign, priorList[["blockIntConcentration"]], control[["nDomains"]])

  processedData <- PostRegressionProcess(initialAssign, domainVector,
                                         dayVector, domainInterest,
                                         initialActivation, control[["nDays"]], 
                                         control[["nDomains"]], nTopics)[[1]];
    activation <<- dailyActivationFullSampler(processedData, exciteParams,  domainRates, priorList[["activationPrior"]], control[["nDays"]], control[["nDomains"]], nTopics)

  processedData <<- PostRegressionProcess(initialAssign, domainVector,
                                         dayVector, domainInterest,
                                         initialActivation, control[["nDays"]], 
                                         control[["nDomains"]], nTopics)[[1]]
  
 print("Completed Initializing Topic Data")
 

 processedWordDataList <<- InitializeData(wordData, c(control[["nDays"]], control[["nDomains"]], nTopics))

 data.list <<- processedWordDataList[[1]]
 domainList <<- processedWordDataList[[2]]

 rm(processedWordDataList)
 rm(wordData)
gc()



 

 
 print("Initialization Done, Beginning Iterations")
 blockAssign <<- initialBlockAssign
 for(iter in 1:iterList[["totalIter"]]){
   iter <<- iter
   time <- proc.time()
   print(paste("Iteration:", iter))
   for(netIter in 1:iterList[["netRegIter"]]){
      networkParams <<- linkRegression(processedLinkData, networkParams, miscList[["netProposals"]], priorList[["netRegMeans"]], priorList[["netRegSD"]]) 
   }
   domainDyads <<- DomainInterestProcess(domainInterest, c(control[["nDays"]], control[["nDomains"]], nTopics))[[1]]
     

  
     domainTopicCounts <<-  processedDataToDomainTopicCount(processedData)
     linkList <<- blockLikelihoodAssign(domainTopicCounts, domainInterest, blockMembNumbers, blockAssign,
                           processedLinkData, networkParams, control, domainDyads, priorList)
     processedLinkData <<- linkList[[2]]
     blockAssign <<- linkList[[1]]
     
   
   processedData <<- PostRegressionProcess(assign, domainVector,
                                          dayVector, domainInterest,
                                          activation, control[["nDays"]], 
                                          control[["nDomains"]], nTopics)[[1]];
  
   postRegParams<<- post.to.topic.regression(processedData, nTopics, exciteParams, 
                                            domainRates, priorList[["PostRegressionMean"]], priorList[["PostRegressionSD"]],
                                            miscList[["PostRegressionProp"]])
   
   exciteParams <<- postRegParams

   
   domainRates <<- domain.base.rate.regression(processedData, control[["nDomains"]], exciteParams, 
                                              domainRates, priorList[["DomainRateMean"]], priorList[["DomainRateSD"]], miscList[["DomainRateProp"]])
   
   domainInterest <<- domain.specific.topic.interest.sampler(processedData, domainInterest, exciteParams, domainRates, 
                                                             priorList[["chi"]], priorList[["domainInterest"]], blockMat, blockAssign, priorList[["blockIntConcentration"]], control[["nDomains"]])
   
   processedData <<- PostRegressionProcess(assign, domainVector,
                                          dayVector, domainInterest,
                                          activation, control[["nDays"]], 
                                          control[["nDomains"]], nTopics)[[1]];
   
   activation <<- dailyActivationFullSampler(processedData, exciteParams,  domainRates, priorList[["activationPrior"]], control[["nDays"]], control[["nDomains"]], nTopics)
   
   domainTopicLik <- domainTopicLikelihoods(domainInterest, activation, exciteParams, domainRates,c(control[["nDays"]], control[["nDomains"]], nTopics))[[1]]
   
   vocab <<- post.assign.full(data.list, domainList, vocab[[2]], vocab[[1]], vocab[[3]], domainTopicLik, priorList[["alpha"]], priorList[["beta"]], miscList[["lag.window"]], c(control[["nDays"]], control[["nDomains"]], nTopics) )
  # KLDiverge <- vocabKLDecomp(vocab, docLength, control[["nDays"]], control[["nTopics"]])
   
   assign <- do.call("c", vocab[[1]])
   gc()
   print(time-proc.time())
   if((iter %% 10 == 0) | (iter == 1)){
     dir.create(paste("./",runName,"/",iter,"/", sep = ""))
     time = proc.time()
     print("SAVING")
     save(list= c("vocab", "networkParams", "exciteParams", "domainRates", "activation","domainInterest", "blockAssign"),file = paste("./",runName,"/",iter,"/readout",iter,".Rdata", sep = ""), compress = T)
     print(proc.time() - time )
   }
   gc()
 }
 
}

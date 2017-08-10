
###Install the relevant packages
setwd("/netscr/trhenry/Duke Project/RunningFolder")
#Load in Data
#Word data can be any number of words, but the final word data matrix has first column date (in integer form), second column domain (in integer form), then the word counts.
wordData <- read.csv("../wordDataTrimmed.csv")[,-1]

#Due to memory use issues, this matrix is specially constructed. It is provided for the whole year of linkage information.
linkData <- read.csv("processedLinkData.csv")[,-1]
processedLinkData <- linkData
#change them to matrices
linkData <- as.matrix(linkData)
wordData <- as.matrix(wordData)

#This loads in all functions and compiles all relevant c++ code
source("FullSampler2.R")
gc()


##These are your control vectors, pretty self explanatory. If you change the number of topics, you need to change the domainInterest Prior
control = list(nDomains = 467, nDays = 366)
nTopics <- 21
priorList = list(alpha =.1, beta=.1, chi = 1, PostRegressionMean =0, PostRegressionSD = 1000, DomainRateMean = 4, DomainRateSD = 1000, domainInterest = rep(1, 21), activationPrior = .2, 
  netRegMeans = c(0,0,0,0,0), netRegSD = c(1000,1000,1000,1000,1000), nBlocksPrior = 25, nMembPrior = 1, blockIntConcentration = 50)
iterList = list(dayone.iter = 10, totalIter = 1000, netRegIter = 10, blockAssignIter = 10)
miscList = list(lag.window = 62, PostRegressionProp = .25, DomainRateProp = .5, netLag = 7, netProposals = c(1,.25,.25,.25,.25))

#Parameter order is: Word data, link data, nTopics, nBlocks, controlVector, IterationList, prior, MiscList, and the directory name for the runs.
topicNetworkSampler(wordData, linkData, 21, 35, control, iterList, priorList, miscList, "Run21")

#This function will produce the information as compressed Rdata files in special folders inside the run directory. You can load these files with load(name of the rData file)

linkData <- read.csv("processedLinkData.csv")[,-1]

linkData[,"outdegree"] <- 0
linkData[,"indegree"] <- 0
for(i in 1:467){
  target <- which(linkData[,5] == i)
  
  outdeg <- sum(linkData[target, 1]*linkData[target,7])/366
  
  linkData[target, 8] <- outdeg
  
}

for(i in 1:467){
  target <- which(linkData[,6] == i)
  
  indeg <- sum(linkData[target, 1]*linkData[target,7])/366
  
  linkData[target, 9] <- indeg
  
}

write.csv(linkData, "processedLinkData.csv")
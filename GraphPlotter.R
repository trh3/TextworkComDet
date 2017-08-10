###Processing Script for Gephi

networkData <- read.csv("newDataOriginalFromDerek.csv")
domainList <- read.csv("domainList.csv", stringsAsFactors = F)

networkData  <- networkData[-1]
networkData  <- networkData[-84]


blogNames <- names(networkData)[-1:-2]
mapping2 <- 1:81
names(mapping2) <- blogNames
domainList[,2] <- gsub("-", ".",x = domainList[,2],fixed = T)

networkData[,1:3]

aggNetData <- aggregate(networkData, by = list(networkData[,2]), FUN = sum)

net <- NULL
for(i in 1:dim(aggNetData)[[1]]){
  for(j in 4:dim(aggNetData)[[2]]){
    print(i)
    net <- rbind(net, c(aggNetData[i,1], domainList[which(domainList[,2] == names(mapping2)[j-3]),1],aggNetData[i,j]))
    
  }
}

head(aggNetData)

net <- net[-which(net[,3] == 0),]

nodes <- domainList
nodes[,"blockAssign"] <- blockAssign
names(nodes) <- c("Id", "Label", "Block")
write.csv(nodes,"NodeList.csv", row.names = F)

net <- as.data.frame(net)
names(net) <- c("Source", "Target", "Weight")
write.csv(net,"EdgesList.csv", row.names = F)


library(igraph)

netGraph <- graph.edgelist(as.matrix(net[,1:2]))
E(netGraph)$weight <- net[,3]
V(netGraph)$label <- domainList[,2]
E(netGraph)$width <- log(E(netGraph)$weight)+1
V(netGraph)$color <- as.numeric(blockAssignF)
netGraph2 <- delete.vertices(netGraph, which(degree(netGraph) == 0))

V(netGraph2)$size = log(degree(netGraph2,mode = "in")+1)+1
netGraph2<- simplify(netGraph2)
E(netGraph2)$width <- log(E(netGraph2)$weight)+1
V(netGraph2)$color
png("NetworkVisualization.png", height = 10000, width = 10000)
plot(netGraph2,layout=layout.fruchterman.reingold(netGraph2, area = 380^5))
dev.off()

which(degree(netGraph, mode = "in")>0)
test <- NULL

blockAssignF <- as.factor(blockAssign)
as.numeric(blockAssignF)



for(i in 1:81){
  
  test <- c(test,which(domainList[,2] == names(mapping2)[i]))
  
}



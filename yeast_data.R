# YEAST-DATA
SwitchOn <- read.csv2("../Data/Average_Switched_On.csv")
SwitchOff <- read.csv2("../Data/Average_Switched_Off.csv")

##Remove Washing period and SE
SwitchOn <- SwitchOn[2:16,1:6]
SwitchOff <- SwitchOff[2:21,1:6]

##Combine data
yeast <- as.data.frame(rbind(SwitchOn,SwitchOff))
rownames(yeast) <- 1:dim(yeast)[1]
 
yeast.gene.scale <- data.frame(scale(yeast[,2:6]),"Time"=as.numeric(as.character(yeast[,1])),"Switch"=c(rep("On",15),rep("Off",20)))
yeast.gene.mean.scale <- data.frame(t(t(yeast[,2:6])-apply(yeast[,2:6],2,mean)),"Time"=as.numeric(as.character(yeast[,1])),"Switch"=c(rep("On",15),rep("Off",20)))
yeast.scale <- data.frame(matrix(scale(as.numeric(as.matrix(yeast[,2:6]))),ncol=5),"Time"=as.numeric(as.character(yeast[,1])),"Switch"=c(rep("On",15),rep("Off",20)))
colnames(yeast.scale)[1:5] <- colnames(yeast)[2:6] 

yeast <- data.frame(yeast[,2:6],"Time"=as.numeric(as.character(yeast[,1])),"Switch"=c(rep("On",15),rep("Off",20)))

##Time dependency
yeast.time <- yeast.gene.scale.time <- yeast.gene.mean.scale.time <- yeast.scale.time <- list()
for(i in 1:5){
	yeast.gene.scale.time[[i]] <- as.data.frame(cbind(yeast.gene.scale[c(2:15,17:35),i],yeast.gene.scale[c(1:14,16:34),c(-i,-6,-7)],yeast.gene.scale[c(2:15,17:35),6:7]))
	colnames(yeast.gene.scale.time[[i]])[1] <- as.character(paste(colnames(yeast)[i]))
	
	yeast.gene.mean.scale.time[[i]] <- as.data.frame(cbind(yeast.gene.mean.scale[c(2:15,17:35),i],yeast.gene.mean.scale[c(1:14,16:34),c(-i,-6,-7)],yeast.gene.mean.scale[c(2:15,17:35),6:7]))
	colnames(yeast.gene.mean.scale.time[[i]])[1] <- as.character(paste(colnames(yeast)[i]))
	
	yeast.scale.time[[i]] <- as.data.frame(cbind(yeast.scale[c(2:15,17:35),i],yeast.scale[c(1:14,16:34),c(-i,-6,-7)],yeast.scale[c(2:15,17:35),6:7]))
	colnames(yeast.scale.time[[i]])[1] <- as.character(paste(colnames(yeast)[i]))
	
	yeast.time[[i]] <- as.data.frame(cbind(yeast[c(2:15,17:35),i],yeast[c(1:14,16:34),c(-i,-6,-7)],yeast[c(2:15,17:35),6:7]))
	colnames(yeast.time[[i]])[1] <- as.character(paste(colnames(yeast)[i]))
}




PreProcessedOn <- t(readMat("../Data/data_on.mat")$data.on)
PreProcessedOff <- t(readMat("../Data/data_off.mat")$data.off)
colnames(PreProcessedOn) <- colnames(PreProcessedOff) <- c("CBF1","GAL4","SWI5","GAL80","ASH1")
PreProcessed.time <- PreProcessed <- list()
for(i in 1:5){
	PreProcessed.time[[i]] <- data.frame(c(PreProcessedOn[2:15,i],PreProcessedOff[2:20,i]),rbind(PreProcessedOn[1:14,-i],PreProcessedOff[1:19,-i]),"Switch"=c(rep("On",14),rep("Off",19)))
	PreProcessed[[i]] <- data.frame(c(PreProcessedOn[1:15,i],PreProcessedOff[1:20,i]),rbind(PreProcessedOn[1:15,-i],PreProcessedOff[1:20,-i]),"Switch"=c(rep("On",15),rep("Off",20)))
	colnames(PreProcessed.time[[i]])[1] <- colnames(PreProcessed[[i]])[1] <- colnames(PreProcessedOn)[i]
}

PreProcessedMult.y <- c()
PreProcessedMult.X <- matrix(0,ncol=25,nrow=165)
for(i in 2:15){
	PreProcessedMult.y <- c(PreProcessedMult.y,PreProcessedOn[i,])
	PreProcessedMult.X[1 + 5*(i-2),1:25] <- c(1,PreProcessedOn[i-1,-1],rep(0,20))
	PreProcessedMult.X[2 + 5*(i-2),1:25] <- c(rep(0,5),1,PreProcessedOn[i-1,-2],rep(0,15))
	PreProcessedMult.X[3 + 5*(i-2),] <- c(rep(0,10),1,PreProcessedOn[i-1,-3],rep(0,10))
	PreProcessedMult.X[4 + 5*(i-2),] <- c(rep(0,15),1,PreProcessedOn[i-1,-4],rep(0,5))
	PreProcessedMult.X[5 + 5*(i-2),] <- c(rep(0,20),1,PreProcessedOn[i-1,-5])
}
for(i in 2:20){
	PreProcessedMult.y <- c(PreProcessedMult.y,PreProcessedOff[i,])
	PreProcessedMult.X[71 + 5*(i-2),] <- c(1,PreProcessedOff[i-1,-1],rep(0,20))
	PreProcessedMult.X[72 + 5*(i-2),] <- c(rep(0,5),1,PreProcessedOff[i-1,-2],rep(0,15))
	PreProcessedMult.X[73 + 5*(i-2),] <- c(rep(0,10),1,PreProcessedOff[i-1,-3],rep(0,10))
	PreProcessedMult.X[74 + 5*(i-2),] <- c(rep(0,15),1,PreProcessedOff[i-1,-4],rep(0,5))
	PreProcessedMult.X[75 + 5*(i-2),] <- c(rep(0,20),1,PreProcessedOff[i-1,-5])
}


adjmTrue <- matrix(c(NA,1,0,0,0,0,NA,1,1,0,1,0,NA,1,1,0,1,0,NA,0,1,0,0,0,NA),ncol=5,byrow=T)
rownames(adjmTrue) <- colnames(adjmTrue) <- c("CBF1","GAL4","SWI5","GAL80","ASH1")

RegressionPosibilities <- matrix(c(2,0,0,0,0,3,0,0,0,0,4,0,0,0,0,5,2,3,0,0,2,0,4,0,2,0,0,5,0,3,4,0,0,3,0,5,0,0,4,5,2,3,4,0,2,3,0,5,2,0,4,5,0,3,4,5,2,3,4,5),ncol=4,byrow=T)


gTrue <- graph.adjacency(adjmTrue)
CL <- layout.circle(gTrue)

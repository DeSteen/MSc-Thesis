####Creating artificial data for model comparison (Chapter 5)
library(MASS)
library(xtable)
library(igraph)
library(PRROC)

source("auroc.R")
source("BaumWelchReg.R")
source("BwdAlgReg.R")
source("ExtractSNSModels.R")
source("FwdAlgReg.R")
source("gplotyeast.R")
source("MultQSNS.R")
source("SNS.R")
source("StepGraph.R")
source("ViterbiReg.R")

set.seed(42)
NodeAHOM <- NodeAHMM <- NodeACPM <- NodeACLA <- 1
NodeBHOM <- NodeBHMM <- NodeBCPM <- NodeBCLA <- 1
NodeCHOM <- NodeCHMM <- NodeCCPM <- NodeCCLA <- 1
NodeDHOM <- NodeDHMM <- NodeDCPM <- NodeDCLA <- 1
NodeEHOM <- NodeEHMM <- NodeECPM <- NodeECLA <- 1
sigma <- 1
beta <- matrix(rnorm(30,0,4*sigma),ncol=10,nrow=3,byrow=T)
for(i in 1:3){
	beta[i,2:3] <- beta[i,2:3]/sqrt(sum(beta[i,2:3]^2))
	beta[i,4:5] <- beta[i,4:5]/sqrt(sum(beta[i,4:5]^2))
	beta[i,6:8] <- beta[i,6:8]/sqrt(sum(beta[i,6:8]^2))
	beta[i,9:10] <- beta[i,9:10]/sqrt(sum(beta[i,9:10]^2))
}

RegressionPosibilities <- matrix(c(0,0,0,0,2,0,0,0,0,3,0,0,0,0,4,0,0,0,0,5,2,3,0,0,2,0,4,0,2,0,0,5,0,3,4,0,0,3,0,5,0,0,4,5,2,3,4,0,2,3,0,5,2,0,4,5,0,3,4,5,2,3,4,5),ncol=4,byrow=T)

gplotad(adjmTrue <- matrix(c(NA,0,0,0,0,1,NA,0,0,0,1,0,NA,0,0,1,0,1,NA,0,0,1,0,0,NA),ncol=5))


sHOM <- rep(1,4*12); sCPM <- rep(c(1,2,3),rep(16,3)); sHMM <- rep(c(1,2,3,1,2,3),c(8,8,8,8,8,8)); sCLA <-c()
A <- abs(beta[1:3,1:3]); A[1,] <- A[1,]/sum(A[1,]); A[2,] <- A[2,]/sum(A[2,]); A[3,] <- A[3,]/sum(A[3,])
p <- A[,2]/sum(A[,2])
sCLA[1] <- 1
for(k in 2:48){
	sCLA[k] <- sample(1:3,1,prob=A[sHMM[k-1],])
	NodeAHOM[k] <- beta[sHOM[k],1] + rnorm(1,0,sigma)
	NodeBHOM[k] <- c(1,NodeAHOM[k-1])%*%beta[sHOM[k],2:3] + rnorm(1,0,sigma)
	NodeCHOM[k] <- c(1,NodeAHOM[k-1])%*%beta[sHOM[k],4:5] + rnorm(1,0,sigma)
	NodeDHOM[k] <- c(1,NodeAHOM[k-1],NodeCHOM[k-1])%*%beta[sHOM[k],6:8] + rnorm(1,0,sigma)
	NodeEHOM[k] <- c(1,NodeBHOM[k-1])%*%beta[sHOM[k],9:10] + rnorm(1,0,sigma)
	
	NodeAHMM[k] <- beta[sHMM[k],1] + rnorm(1,0,sigma)
	NodeBHMM[k] <- c(1,NodeAHMM[k-1])%*%beta[sHMM[k],2:3] + rnorm(1,0,sigma)
	NodeCHMM[k] <- c(1,NodeAHMM[k-1])%*%beta[sHMM[k],4:5] + rnorm(1,0,sigma)
	NodeDHMM[k] <- c(1,NodeAHMM[k-1],NodeCHMM[k-1])%*%beta[sHMM[k],6:8] + rnorm(1,0,sigma)
	NodeEHMM[k] <- c(1,NodeBHMM[k-1])%*%beta[sHMM[k],9:10] + rnorm(1,0,sigma)
	
	NodeACPM[k] <- beta[sCPM[k],1] + rnorm(1,0,sigma)
	NodeBCPM[k] <- c(1,NodeACPM[k-1])%*%beta[sCPM[k],2:3] + rnorm(1,0,sigma)
	NodeCCPM[k] <- c(1,NodeACPM[k-1])%*%beta[sCPM[k],4:5] + rnorm(1,0,sigma)
	NodeDCPM[k] <- c(1,NodeACPM[k-1],NodeCCPM[k-1])%*%beta[sCPM[k],6:8] + rnorm(1,0,sigma)
	NodeECPM[k] <- c(1,NodeBCPM[k-1])%*%beta[sCPM[k],9:10] + rnorm(1,0,sigma)
} 

HOMDat <- data.frame(NodeAHOM,NodeBHOM,NodeCHOM,NodeDHOM,NodeEHOM)
HMMDat <- data.frame(NodeAHMM,NodeBHMM,NodeCHMM,NodeDHMM,NodeEHMM)
CPMDat <- data.frame(NodeACPM,NodeBCPM,NodeCCPM,NodeDCPM,NodeECPM)

hom.dat <- hmm.dat <- cpm.dat <- cla.dat <- list()
for(i in 1:5){
	hom.dat[[i]] <- cbind(HOMDat[2:48,i],HOMDat[1:47,-i])
	hmm.dat[[i]] <- cbind(HMMDat[2:48,i],HMMDat[1:47,-i])
	cpm.dat[[i]] <- cbind(CPMDat[2:48,i],CPMDat[1:47,-i])
}


### Linear Models
#Homogeneous sequence
hom.lin.mod <- list()
hom.lin.mod_score <- hom.lin.mod_BIC <- matrix(rep(0,80),ncol=16)
	for(i in 1:5){
		hom.lin.mod[[i]] <- list()
		for(j in 1:16){
			hom.lin.mod[[i]][[j]] <- lm(as.numeric(hom.dat[[i]][,1])~cbind(1,as.matrix(hom.dat[[i]][,c(RegressionPosibilities[j,])])) -1 )
			hom.lin.mod_BIC[i,j] <- BIC(hom.lin.mod[[i]][[j]])
		}
		hom.lin.mod_score[i,which.min(hom.lin.mod_BIC[i,])] <- 1 
	}

hom.lin.mod_adjm <- matrix(rep(NA,25),ncol=5)
for(i in 1:5){
	hom.lin.mod_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*hom.lin.mod_score[i,],2,sum))
}
par(mfrow=c(2,2))
gplotad(t(hom.lin.mod_adjm))
quartz()
auroc(hom.lin.mod_adjm,t(adjmTrue),"hom.lin.mod")


#Changepoint sequence
cpm.lin.mod <- list()
cpm.lin.mod_score <- cpm.lin.mod_BIC <- matrix(rep(0,80),ncol=16)
	for(i in 1:5){
		cpm.lin.mod[[i]] <- list()
		for(j in 1:16){
			cpm.lin.mod[[i]][[j]] <- lm(as.numeric(cpm.dat[[i]][,1])~cbind(1,as.matrix(cpm.dat[[i]][,c(RegressionPosibilities[j,])])) -1 )
			cpm.lin.mod_BIC[i,j] <- BIC(cpm.lin.mod[[i]][[j]])
		}
		cpm.lin.mod_score[i,which.min(cpm.lin.mod_BIC[i,])] <- 1 
	}

cpm.lin.mod_adjm <- matrix(rep(NA,25),ncol=5)
for(i in 1:5){
	cpm.lin.mod_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*cpm.lin.mod_score[i,],2,sum))
}
par(mfrow=c(2,2))
gplotad(t(cpm.lin.mod_adjm))
auroc(cpm.lin.mod_adjm,t(adjmTrue),"cpm.lin.mod")


#Hidden Markov Sequence
hmm.lin.mod <- list()
hmm.lin.mod_score <- hmm.lin.mod_BIC <- matrix(rep(0,80),ncol=16)
	for(i in 1:5){
		hmm.lin.mod[[i]] <- list()
		for(j in 1:16){
			hmm.lin.mod[[i]][[j]] <- lm(as.numeric(hmm.dat[[i]][,1])~cbind(1,as.matrix(hmm.dat[[i]][,c(RegressionPosibilities[j,])])) -1 )
			hmm.lin.mod_BIC[i,j] <- BIC(hmm.lin.mod[[i]][[j]])
		}
		hmm.lin.mod_score[i,which.min(hmm.lin.mod_BIC[i,])] <- 1 
	}

hmm.lin.mod_adjm <- matrix(rep(NA,25),ncol=5)
for(i in 1:5){
	hmm.lin.mod_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*hmm.lin.mod_score[i,],2,sum))
}
par(mfrow=c(2,2))
gplotad(t(hmm.lin.mod_adjm))
auroc(hmm.lin.mod_adjm,t(adjmTrue),"hmm.lin.mod")
matplot(1:16,t(hmm.lin.mod_BIC),col=1:5,xlab="Parent option",ylab="BIC",type="p",xaxt="n")
axis(1,at=1:16)



##Changepoint Models
#Linear sequence
hom.cpm.mod.BIC <- hom.cpm.mod.MDL <- list()
hom.cpm.mod.BIC_score <- hom.cpm.mod.BIC_BIC <- hom.cpm.mod.MDL_score <- hom.cpm.mod.MDL_BIC <- matrix(rep(0,80),ncol=16)
for(i in 1:5){
	hom.cpm.mod.BIC[[i]] <- hom.cpm.mod.MDL[[i]] <- list()
	for(j in 1:16){
		hom.cpm.mod.BIC[[i]][[j]] <- SNS(as.numeric(hom.dat[[i]][,1]),cbind(1,as.matrix(hom.dat[[i]][,c(RegressionPosibilities[j,])])),8)
		hom.cpm.mod.MDL[[i]][[j]] <- SNS(as.numeric(hom.dat[[i]][,1]),cbind(1,as.matrix(hom.dat[[i]][,c(RegressionPosibilities[j,])])),8,mincrit="MDL")
			
		hom.cpm.mod.BIC_BIC[i,j] <- hom.cpm.mod.BIC[[i]][[j]]$BIC
		hom.cpm.mod.MDL_BIC[i,j] <- hom.cpm.mod.MDL[[i]][[j]]$BIC
			
		cat("\r node ",i," posibility ",formatC(j,width=2))
	}
	hom.cpm.mod.BIC_score[i,which.min(hom.cpm.mod.BIC_BIC[i,])] <- 1 
	hom.cpm.mod.MDL_score[i,which.min(hom.cpm.mod.MDL_BIC[i,])] <- 1
}

hom.cpm.mod.BIC_adjm <- hom.cpm.mod.MDL_adjm <- matrix(rep(NA,25),ncol=5)
for(i in 1:5){
	hom.cpm.mod.BIC_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*hom.cpm.mod.BIC_score[i,],2,sum))
	hom.cpm.mod.MDL_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*hom.cpm.mod.MDL_score[i,],2,sum))
}

hom.cpm.mod.glob.BIC <- StepGraph(hom.dat,type="SNSFull",InitMat=t(adjmTrue))
hom.cpm.mod.glob.MDL <- StepGraph(hom.dat,type="SNSFull",InitMat=hom.cpm.mod.glob.BIC$Adjency,mincrit="MDL")

par(mfrow=c(2,2))
gplotad(t(hom.cpm.mod.BIC_adjm))
auroc(hom.cpm.mod.BIC_adjm,t(adjmTrue),"hom.cpm.mod.BIC")
gplotad(t(hom.cpm.mod.glob.BIC$Adjency))
auroc(hom.cpm.mod.glob.BIC$Adjency,t(adjmTrue),"hom.cpm.mod.glob.BIC")

par(mfrow=c(2,2))
gplotad(t(hom.cpm.mod.MDL_adjm))
auroc(hom.cpm.mod.MDL_adjm,t(adjmTrue),"hom.cpm.mod.BIC")
gplotad(t(hom.cpm.mod.glob.MDL$Adjency))
auroc(hom.cpm.mod.glob.MDL$Adjency,t(adjmTrue),"hom.cpm.mod.glob.MDL")

hom.cpm.mod.MDL_score
hom.cpm.mod.MDL[[5]][[1]]$Tau

confint(lm(as.numeric(hom.dat[[4]][,1])~cbind(1,as.matrix(hom.dat[[4]][,c(RegressionPosibilities[7,])]))-1))


# Changepoint sequence
cpm.cpm.mod.BIC <- cpm.cpm.mod.MDL <- list()
cpm.cpm.mod.BIC_score <- cpm.cpm.mod.BIC_BIC <- cpm.cpm.mod.MDL_score <- cpm.cpm.mod.MDL_BIC <- matrix(rep(0,80),ncol=16)
for(i in 1:5){
	cpm.cpm.mod.BIC[[i]] <- cpm.cpm.mod.MDL[[i]] <- list()
	for(j in 1:16){
		cpm.cpm.mod.BIC[[i]][[j]] <- SNS(as.numeric(cpm.dat[[i]][,1]),cbind(1,as.matrix(cpm.dat[[i]][,c(RegressionPosibilities[j,])])),8)
		cpm.cpm.mod.MDL[[i]][[j]] <- SNS(as.numeric(cpm.dat[[i]][,1]),cbind(1,as.matrix(cpm.dat[[i]][,c(RegressionPosibilities[j,])])),8,mincrit="MDL")
			
		cpm.cpm.mod.BIC_BIC[i,j] <- cpm.cpm.mod.BIC[[i]][[j]]$BIC
		cpm.cpm.mod.MDL_BIC[i,j] <- cpm.cpm.mod.MDL[[i]][[j]]$BIC
			
		cat("\r node ",i," posibility ",formatC(j,width=2))
	}
	cpm.cpm.mod.BIC_score[i,which.min(cpm.cpm.mod.BIC_BIC[i,])] <- 1 
	cpm.cpm.mod.MDL_score[i,which.min(cpm.cpm.mod.MDL_BIC[i,])] <- 1
}

cpm.cpm.mod.BIC_adjm <- cpm.cpm.mod.MDL_adjm <- matrix(rep(NA,25),ncol=5)
for(i in 1:5){
	cpm.cpm.mod.BIC_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*cpm.cpm.mod.BIC_score[i,],2,sum))
	cpm.cpm.mod.MDL_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*cpm.cpm.mod.MDL_score[i,],2,sum))
	print(cpm.cpm.mod.BIC[[i]][[which.min(cpm.cpm.mod.BIC_BIC[i,])]]$Tau)
	print(cpm.cpm.mod.MDL[[i]][[which.min(cpm.cpm.mod.MDL_BIC[i,])]]$Tau)
}

cpm.cpm.mod.glob.BIC <- StepGraph(cpm.dat,type="SNSFull",InitMat=t(adjmTrue))
cpm.cpm.mod.glob.MDL <- StepGraph(cpm.dat,type="SNSFull",InitMat=cpm.cpm.mod.glob.BIC$Adjency,mincrit="MDL")

par(mfrow=c(2,2))
gplotad(t(cpm.cpm.mod.BIC_adjm))
auroc(cpm.cpm.mod.BIC_adjm,t(adjmTrue),"cpm.cpm.mod.BIC")
gplotad(t(cpm.cpm.mod.glob.BIC$Adjency))
auroc(cpm.cpm.mod.glob.BIC$Adjency,t(adjmTrue),"cpm.cpm.mod.glob.BIC")

par(mfrow=c(2,2))
gplotad(t(cpm.cpm.mod.MDL_adjm))
auroc(cpm.cpm.mod.MDL_adjm,t(adjmTrue),"cpm.cpm.mod.BIC")
gplotad(t(cpm.cpm.mod.glob.MDL$Adjency))
auroc(cpm.cpm.mod.glob.MDL$Adjency,t(adjmTrue),"cpm.cpm.mod.glob.MDL")

cpm.cpm.mod.glob.MDL

opt <- c(1,1,2,7,3)
for(i in 1:5){
	a <- lm(as.numeric(cpm.dat[[i]][,1])~cbind(1,as.matrix(cpm.dat[[i]][,c(RegressionPosibilities[opt[i],])]))-1,subset=1:15)
	b <- lm(as.numeric(cpm.dat[[i]][,1])~cbind(1,as.matrix(cpm.dat[[i]][,c(RegressionPosibilities[opt[i],])]))-1,subset=16:31)
	d <- lm(as.numeric(cpm.dat[[i]][,1])~cbind(1,as.matrix(cpm.dat[[i]][,c(RegressionPosibilities[opt[i],])]))-1,subset=32:47)
	print(cbind(coef(a),confint(a)))
	print(cbind(coef(b),confint(b)))
	print(cbind(coef(d),confint(d)))
}


#Hidden Markov Sequence
hmm.cpm.mod.BIC <- hmm.cpm.mod.MDL <- list()
hmm.cpm.mod.BIC_score <- hmm.cpm.mod.BIC_BIC <- hmm.cpm.mod.MDL_score <- hmm.cpm.mod.MDL_BIC <- matrix(rep(0,80),ncol=16)
for(i in 1:5){
	hmm.cpm.mod.BIC[[i]] <- hmm.cpm.mod.MDL[[i]] <- list()
	for(j in 1:16){
		hmm.cpm.mod.BIC[[i]][[j]] <- SNS(as.numeric(hmm.dat[[i]][,1]),cbind(1,as.matrix(hmm.dat[[i]][,c(RegressionPosibilities[j,])])),8)
		hmm.cpm.mod.MDL[[i]][[j]] <- SNS(as.numeric(hmm.dat[[i]][,1]),cbind(1,as.matrix(hmm.dat[[i]][,c(RegressionPosibilities[j,])])),8,mincrit="MDL")
			
		hmm.cpm.mod.BIC_BIC[i,j] <- hmm.cpm.mod.BIC[[i]][[j]]$BIC
		hmm.cpm.mod.MDL_BIC[i,j] <- hmm.cpm.mod.MDL[[i]][[j]]$BIC
			
		cat("\r node ",i," posibility ",formatC(j,width=2))
	}
	hmm.cpm.mod.BIC_score[i,which.min(hmm.cpm.mod.BIC_BIC[i,])] <- 1 
	hmm.cpm.mod.MDL_score[i,which.min(hmm.cpm.mod.MDL_BIC[i,])] <- 1
}

hmm.cpm.mod.BIC_adjm <- hmm.cpm.mod.MDL_adjm <- matrix(rep(NA,25),ncol=5)
for(i in 1:5){
	hmm.cpm.mod.BIC_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*hmm.cpm.mod.BIC_score[i,],2,sum))
	hmm.cpm.mod.MDL_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*hmm.cpm.mod.MDL_score[i,],2,sum))
	print(hmm.cpm.mod.BIC[[i]][[which.min(hmm.cpm.mod.BIC_BIC[i,])]]$Tau)
	print(hmm.cpm.mod.MDL[[i]][[which.min(hmm.cpm.mod.MDL_BIC[i,])]]$Tau)
}

hmm.cpm.mod.glob.BIC <- StepGraph(hmm.dat,type="SNSFull",InitMat=t(adjmTrue))
hmm.cpm.mod.glob.MDL <- StepGraph(hmm.dat,type="SNSFull",InitMat=hmm.cpm.mod.glob.BIC$Adjency,mincrit="MDL")

par(mfrow=c(2,2))
gplotad(t(hmm.cpm.mod.BIC_adjm))
auroc(hmm.cpm.mod.BIC_adjm,t(adjmTrue),"hmm.cpm.mod.BIC")
gplotad(t(hmm.cpm.mod.glob.BIC$Adjency))
auroc(hmm.cpm.mod.glob.BIC$Adjency,t(adjmTrue),"hmm.cpm.mod.glob.BIC")

par(mfrow=c(2,2))
gplotad(t(hmm.cpm.mod.MDL_adjm))
auroc(hmm.cpm.mod.MDL_adjm,t(adjmTrue),"hmm.cpm.mod.BIC")
gplotad(t(hmm.cpm.mod.glob.MDL$Adjency))
auroc(hmm.cpm.mod.glob.MDL$Adjency,t(adjmTrue),"hmm.cpm.mod.glob.MDL")




##Hidden Markov Model
#Changepoint sequence
cpm.hmm.mod <- cpm.hmm.mod_EqVar <- cpm.hmm.mod_EqVar_LR <- cpm.hmm.mod_LR <- list()
for(n in 1:100){
	cpm.hmm.mod[[n]] <- cpm.hmm.mod_EqVar[[n]] <- cpm.hmm.mod_LR[[n]] <- cpm.hmm.mod_EqVar_LR[[n]] <- list()
	for(i in 1:5){
		cpm.hmm.mod[[n]][[i]] <- cpm.hmm.mod_EqVar[[n]][[i]] <- cpm.hmm.mod_LR[[n]][[i]] <- cpm.hmm.mod_EqVar_LR[[n]][[i]] <- list()
		for(j in 1:16){
			cat("\r iteration",formatC(n,width=3),"node",i,"Regression Posibility",formatC(j,width=2))
			cpm.hmm.mod[[n]][[i]][[j]] <- BaumWelchReg(as.numeric(cpm.dat[[i]][,1]),cbind(1,as.matrix(cpm.dat[[i]][,c(RegressionPosibilities[j,])])),3,-1e9,1)
			cpm.hmm.mod_EqVar[[n]][[i]][[j]] <- BaumWelchReg(as.numeric(cpm.dat[[i]][,1]),cbind(1,as.matrix(cpm.dat[[i]][,c(RegressionPosibilities[j,])])),3,-1e9,1,eq.var=T)
			cpm.hmm.mod_LR[[n]][[i]][[j]] <- BaumWelchReg(as.numeric(cpm.dat[[i]][,1]),cbind(1,as.matrix(cpm.dat[[i]][,c(RegressionPosibilities[j,])])),3,-1e9,1,LR=T)
			cpm.hmm.mod_EqVar_LR[[n]][[i]][[j]] <- BaumWelchReg(as.numeric(cpm.dat[[i]][,1]),cbind(1,as.matrix(cpm.dat[[i]][,c(RegressionPosibilities[j,])])),3,-1e9,1,LR=T,eq.var=T)
		}
	}
}

cpm.hmm.mod_BIC <- cpm.hmm.mod_EqVar_BIC <- cpm.hmm.mod_LR_BIC <- cpm.hmm.mod_EqVar_LR_BIC <- array(0,dim=c(5,16,100))
cpm.hmm.mod_score <- cpm.hmm.mod_EqVar_score <- cpm.hmm.mod_LR_score <- cpm.hmm.mod_EqVar_LR_score <- matrix(0,ncol=16,nrow=5)
cpm.hmm.mod_adjm <- cpm.hmm.mod_EqVar_adjm <- cpm.hmm.mod_LR_adjm <- cpm.hmm.mod_EqVar_LR_adjm <- matrix(rep(NA,25),ncol=5)
for(i in 1:5){
	for(n in 1:100){
		for(j in 1:16){
			cpm.hmm.mod_BIC[i,j,n] <- cpm.hmm.mod[[n]][[i]][[j]]$BIC
			cpm.hmm.mod_EqVar_BIC[i,j,n] <- cpm.hmm.mod_EqVar[[n]][[i]][[j]]$BIC
			cpm.hmm.mod_LR_BIC[i,j,n] <- cpm.hmm.mod_LR[[n]][[i]][[j]]$BIC
			cpm.hmm.mod_EqVar_LR_BIC[i,j,n] <- cpm.hmm.mod_EqVar_LR[[n]][[i]][[j]]$BIC
		}
	}
	cpm.hmm.mod_score[i,] <- apply(cpm.hmm.mod_BIC[i,,] == min(cpm.hmm.mod_BIC[i,,]),1,sum) 
	cpm.hmm.mod_EqVar_score[i,] <- apply(cpm.hmm.mod_EqVar_BIC[i,,] == min(cpm.hmm.mod_EqVar_BIC[i,,]),1,sum) 
	cpm.hmm.mod_LR_score[i,] <- apply(cpm.hmm.mod_LR_BIC[i,,] == min(cpm.hmm.mod_LR_BIC[i,,]),1,sum) 
	cpm.hmm.mod_EqVar_LR_score[i,] <- apply(cpm.hmm.mod_EqVar_LR_BIC[i,,] == min(cpm.hmm.mod_EqVar_LR_BIC[i,,]),1,sum) 
		
	cpm.hmm.mod_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*(cpm.hmm.mod_score[i,]!=0),2,sum))
	cpm.hmm.mod_EqVar_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*(cpm.hmm.mod_EqVar_score[i,]!=0),2,sum))
	cpm.hmm.mod_LR_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*(cpm.hmm.mod_LR_score[i,]!=0),2,sum))
	cpm.hmm.mod_EqVar_LR_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*(cpm.hmm.mod_EqVar_LR_score[i,]!=0),2,sum))
}

par(mfrow=c(2,2))
gplotad(t(cpm.hmm.mod_adjm))
auroc(cpm.hmm.mod_adjm,t(adjmTrue),"cpm.hmm.mod")
gplotad(t(cpm.hmm.mod_LR_adjm))
auroc(cpm.hmm.mod_LR_adjm,t(adjmTrue),"cpm.hmm.mod")

cpm.hmm.mod_seq1 <- matrix(0,nrow=12,ncol=47); colnames(cpm.hmm.mod_seq1) <- 2:48
for(i in 1:5){
opt <- which(cpm.hmm.mod_BIC[i,,] == min(cpm.hmm.mod_BIC[i,,]),arr.ind=T)
cpm.hmm.mod_seq1[i,] <- ViterbiReg(as.numeric(cpm.dat[[i]][,1]),cbind(1,as.matrix(cpm.dat[[i]][,c(RegressionPosibilities[opt[1,1],])])),cpm.hmm.mod[[opt[1,2]]][[i]][[opt[1,1]]]$A,cpm.hmm.mod[[opt[1,2]]][[i]][[opt[1,1]]]$B,cpm.hmm.mod[[opt[1,2]]][[i]][[opt[1,1]]]$pi)
}

for(i in 1:5){
opt <- which(cpm.hmm.mod_LR_BIC[i,,] == min(cpm.hmm.mod_LR_BIC[i,,]),arr.ind=T)
cpm.hmm.mod_seq1[i+7,] <- ViterbiReg(as.numeric(cpm.dat[[i]][,1]),cbind(1,as.matrix(cpm.dat[[i]][,c(RegressionPosibilities[opt[1,1],])])),cpm.hmm.mod_LR[[opt[1,2]]][[i]][[opt[1,1]]]$A,cpm.hmm.mod_LR[[opt[1,2]]][[i]][[opt[1,1]]]$B,cpm.hmm.mod_LR[[opt[1,2]]][[i]][[opt[1,1]]]$pi)
}

xtable(t(cpm.hmm.mod_seq1),digits=0)

for(i in 1:5){
opt <- which(cpm.hmm.mod_EqVar_LR_BIC[i,,] == min(cpm.hmm.mod_EqVar_LR_BIC[i,,]),arr.ind=T)
print(ViterbiReg(as.numeric(cpm.dat[[i]][,1]),cbind(1,as.matrix(cpm.dat[[i]][,c(RegressionPosibilities[opt[1,1],])])),cpm.hmm.mod_EqVar_LR[[opt[1,2]]][[i]][[opt[1,1]]]$A,cpm.hmm.mod_EqVar_LR[[opt[1,2]]][[i]][[opt[1,1]]]$B,cpm.hmm.mod_EqVar_LR[[opt[1,2]]][[i]][[opt[1,1]]]$pi))
}

boxplot(t(cpm.hmm.mod_BIC[1,,]))
quartz()
boxplot(t(cpm.hmm.mod_EqVar_BIC[1,,]))

par(mfrow=c(2,2))
gplotad(t(cpm.hmm.mod_EqVar_adjm))
auroc(cpm.hmm.mod_EqVar_adjm,t(adjmTrue),"cpm.hmm.mod")
gplotad(t(cpm.hmm.mod_EqVar_LR_adjm))
auroc(cpm.hmm.mod_EqVar_LR_adjm,t(adjmTrue),"cpm.hmm.mod")

boxplot(t(cpm.hmm.mod_LR_BIC[1,,]))
quartz()
boxplot(t(cpm.hmm.mod_EqVar_LR_BIC[1,,]))


#Hidden Markov Sequence
hmm.hmm.mod <- hmm.hmm.mod_EqVar <- list()
for(n in 1:100){
	hmm.hmm.mod[[n]] <- hmm.hmm.mod_EqVar[[n]] <- list()
	for(i in 1:5){
		hmm.hmm.mod[[n]][[i]] <- hmm.hmm.mod_EqVar[[n]][[i]] <- list()
		for(j in 1:16){
			cat("\r iteration",formatC(n,width=3),"node",i,"Regression Posibility",formatC(j,width=2))
			hmm.hmm.mod[[n]][[i]][[j]] <- BaumWelchReg(as.numeric(hmm.dat[[i]][,1]),cbind(1,as.matrix(hmm.dat[[i]][,c(RegressionPosibilities[j,])])),3,-1e9,1,coef.traj=T)
			hmm.hmm.mod_EqVar[[n]][[i]][[j]] <- BaumWelchReg(as.numeric(hmm.dat[[i]][,1]),cbind(1,as.matrix(hmm.dat[[i]][,c(RegressionPosibilities[j,])])),3,-1e9,1,eq.var=T,coef.traj=T)
		}
	}
}

hmm.hmm.mod_BIC <- hmm.hmm.mod_EqVar_BIC <- array(0,dim=c(5,16,100))
hmm.hmm.mod_score <- hmm.hmm.mod_EqVar_score <- hmm.hmm.mod_score_outlier <- hmm.hmm.mod_EqVar_score_outlier <- matrix(0,ncol=16,nrow=5)
hmm.hmm.mod_adjm <- hmm.hmm.mod_EqVar_adjm <- hmm.hmm.mod_adjm_outlier <- hmm.hmm.mod_EqVar_adjm_outlier <- matrix(rep(NA,25),ncol=5)
for(i in 1:5){
	for(n in 1:100){
		for(j in 1:16){
			hmm.hmm.mod_BIC[i,j,n] <- hmm.hmm.mod[[n]][[i]][[j]]$BIC
			hmm.hmm.mod_EqVar_BIC[i,j,n] <- hmm.hmm.mod_EqVar[[n]][[i]][[j]]$BIC
		}
	}
	hmm.hmm.mod_score[i,] <- apply(hmm.hmm.mod_BIC[i,,] == min(hmm.hmm.mod_BIC[i,,]),1,sum) 
	hmm.hmm.mod_EqVar_score[i,] <- apply(hmm.hmm.mod_EqVar_BIC[i,,] == min(hmm.hmm.mod_EqVar_BIC[i,,]),1,sum) 
	
	hmm.hmm.mod_score_outlier[i,] <- apply(apply(hmm.hmm.mod_BIC[i,,],1,RemoveOutlier) == min(apply(hmm.hmm.mod_BIC[i,,],1,RemoveOutlier),na.rm=T),2,sum,na.rm=T) 
	hmm.hmm.mod_EqVar_score_outlier[i,] <- apply(apply(hmm.hmm.mod_EqVar_BIC[i,,],1,RemoveOutlier) == min(apply(hmm.hmm.mod_EqVar_BIC[i,,],1,RemoveOutlier),na.rm=T),2,sum,na.rm=T)
	
	hmm.hmm.mod_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*(hmm.hmm.mod_score[i,]!=0),2,sum))
	hmm.hmm.mod_EqVar_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*(hmm.hmm.mod_EqVar_score[i,]!=0),2,sum))
	
	hmm.hmm.mod_adjm_outlier[i,-i] <- (apply((RegressionPosibilities!=0)*(hmm.hmm.mod_score_outlier[i,]!=0),2,sum))
	hmm.hmm.mod_EqVar_adjm_outlier[i,-i] <- (apply((RegressionPosibilities!=0)*(hmm.hmm.mod_EqVar_score_outlier[i,]!=0),2,sum))
}

par(mfrow=c(2,2))
gplotad(t(hmm.hmm.mod_adjm))
auroc(hmm.hmm.mod_adjm,t(adjmTrue),"hmm.hmm.mod")
gplotad(t(hmm.hmm.mod_EqVar_adjm))
auroc(hmm.hmm.mod_EqVar_adjm,t(adjmTrue),"hmm.hmm.mod")

hmm.hmm.mod_seq1 <- matrix(0,nrow=12,ncol=47); colnames(hmm.hmm.mod_seq1) <- 2:48
for(i in 1:5){
opt <- which(hmm.hmm.mod_BIC[i,,] == min(hmm.hmm.mod_BIC[i,,]),arr.ind=T)
hmm.hmm.mod_seq1[i,] <- ViterbiReg(as.numeric(hmm.dat[[i]][,1]),cbind(1,as.matrix(hmm.dat[[i]][,c(RegressionPosibilities[opt[1,1],])])),hmm.hmm.mod[[opt[1,2]]][[i]][[opt[1,1]]]$A,hmm.hmm.mod[[opt[1,2]]][[i]][[opt[1,1]]]$B,hmm.hmm.mod[[opt[1,2]]][[i]][[opt[1,1]]]$pi)
}

for(i in 1:5){
opt <- which(hmm.hmm.mod_EqVar_BIC[i,,] == min(hmm.hmm.mod_EqVar_BIC[i,,]),arr.ind=T)
hmm.hmm.mod_seq1[i+7,] <- ViterbiReg(as.numeric(hmm.dat[[i]][,1]),cbind(1,as.matrix(hmm.dat[[i]][,c(RegressionPosibilities[opt[1,1],])])),hmm.hmm.mod_EqVar[[opt[1,2]]][[i]][[opt[1,1]]]$A,hmm.hmm.mod_EqVar[[opt[1,2]]][[i]][[opt[1,1]]]$B,hmm.hmm.mod_EqVar[[opt[1,2]]][[i]][[opt[1,1]]]$pi)
}
xtable(t(hmm.hmm.mod_seq1),digits=0)


#Normal HMM
which(hmm.hmm.mod_BIC[3,2,] ==min(hmm.hmm.mod_BIC[3,2,]),arr.ind=T)
hmm.hmm.mod[[76]][[3]][[2]]$Coef.Traj$B[[length(hmm.hmm.mod[[76]][[3]][[2]]$Coef.Traj$B)]]
ViterbiReg(hmm.dat[[3]][,1],cbind(1,as.matrix(hmm.dat[[3]][,RegressionPosibilities[2,]])),hmm.hmm.mod[[76]][[3]][[2]]$A,hmm.hmm.mod[[76]][[3]][[2]]$B,hmm.hmm.mod[[76]][[3]][[2]]$pi)

#Equal Variance
which(hmm.hmm.mod_EqVar_BIC[3,,] == min(hmm.hmm.mod_EqVar_BIC[3,,]),arr.ind=T)
hmm.hmm.mod_EqVar[[75]][[3]][[2]]$Coef.Traj$B[[length(hmm.hmm.mod_EqVar[[75]][[3]][[2]]$Coef.Traj$B)]]
ViterbiReg(hmm.dat[[3]][,1],cbind(1,as.matrix(hmm.dat[[3]][,RegressionPosibilities[2,]])),hmm.hmm.mod_EqVar[[75]][[3]][[2]]$A,hmm.hmm.mod_EqVar[[75]][[3]][[2]]$B,hmm.hmm.mod_EqVar[[75]][[3]][[2]]$pi)

xN <- 1:length(hmm.hmm.mod[[76]][[3]][[2]]$Coef.Traj$B)
beta1hatN <- beta2hatN <- beta3hatN <- sigmahatN <- matrix(0,nrow=length(xN),ncol=3)
for(i in xN){
	beta1hatN[i,1:3] <-  hmm.hmm.mod[[76]][[3]][[2]]$Coef.Traj$B[[i]][,1]
	beta2hatN[i,1:3] <-  hmm.hmm.mod[[76]][[3]][[2]]$Coef.Traj$B[[i]][,2]
	sigmahatN[i,1:3] <-  log(hmm.hmm.mod[[76]][[3]][[2]]$Coef.Traj$B[[i]][,3])
}

xE <- 1:length(hmm.hmm.mod_EqVar[[75]][[3]][[2]]$Coef.Traj$B)
beta1hatE <- beta2hatE <- sigmahatE <- matrix(0,nrow=length(xE),ncol=3)
for(i in xE){
	beta1hatE[i,1:3] <-  hmm.hmm.mod_EqVar[[75]][[3]][[2]]$Coef.Traj$B[[i]][,1]
	beta2hatE[i,1:3] <-  hmm.hmm.mod_EqVar[[75]][[3]][[2]]$Coef.Traj$B[[i]][,2]
	sigmahatE[i,1:3] <-  hmm.hmm.mod_EqVar[[75]][[3]][[2]]$Coef.Traj$B[[i]][,3]
}

par(mfrow=c(2,2),mar=c(5.1,4.6,4.1,2.1))
plot(xN,beta1hatN[,1],ylim=c(-1,1.5),col=1,type="l",ylab=expression(hat(beta)[paste(C0,",",s)]),xlab="# Iterations",main="Unequal Variances Among Groups")
lines(xN,beta1hatN[,2],col=3)
lines(xN,beta1hatN[,3],col=2)
curve(beta[1,4] + 0*x,from=0,to=29,add=T,lty=2,col=3)
curve(beta[2,4] + 0*x,from=0,to=29,add=T,lty=2,col=1)
curve(beta[3,4] + 0*x,from=0,to=29,add=T,lty=2,col=2)

par(mfg=c(2,1))
plot(xN,beta2hatN[,1],ylim=c(-1,2),col=1,type="l",ylab=expression(hat(beta)[paste(CA,",",s)]),xlab="# Iterations")
lines(xN,beta2hatN[,2],col=3)
lines(xN,beta2hatN[,3],col=2)
curve(beta[1,5] + 0*x,from=0,to=29,add=T,lty=2,col=3)
curve(beta[2,5] + 0*x,from=0,to=29,add=T,lty=2,col=1)
curve(beta[3,5] + 0*x,from=0,to=29,add=T,lty=2,col=2)

par(mfg=c(1,2))
plot(xE,beta1hatE[,1],ylim=c(-1,1.5),col=1,type="l",ylab=expression(hat(beta)[paste(C0,",",s)]),xlab="# Iterations",main="Equal Variances Among Groups")
lines(xE,beta1hatE[,2],col=2)
lines(xE,beta1hatE[,3],col=3)
curve(beta[1,4] + 0*x,from=0,to=65,add=T,lty=2,col=3)
curve(beta[2,4] + 0*x,from=0,to=65,add=T,lty=2,col=1)
curve(beta[3,4] + 0*x,from=0,to=65,add=T,lty=2,col=2)

par(mfg=c(2,2))
plot(xE,beta2hatE[,1],ylim=c(-1,2),col=1,type="l",ylab=expression(hat(beta)[paste(CA,",",s)]),xlab="# Iterations")
lines(xE,beta2hatE[,2],col=2)
lines(xE,beta2hatE[,3],col=3)
curve(beta[1,5] + 0*x,from=0,to=65,add=T,lty=2,col=3)
curve(beta[2,5] + 0*x,from=0,to=65,add=T,lty=2,col=1)
curve(beta[3,5] + 0*x,from=0,to=65,add=T,lty=2,col=2)


par(mfrow=c(2,2))
plot(xN,sigmahatN[,1],ylim=c(-37,1),col=1,type="l",ylab=expression(paste("log(",hat(sigma)[paste(s)],")")),xlab="# Iterations")
lines(xN,sigmahatN[,2],col=2)
lines(xN,sigmahatN[,3],col=3)
curve(0*x,from=0,to=29,add=T,lty=2)

par(mfg=c(2,1))
plot(xN,hmm.hmm.mod[[76]][[3]][[2]]$Log[-1],type="l",ylab="Log-Likelihood",xlab="# Iterations",ylim=c(-95,-5))

par(mfg=c(1,2))
plot(xE,log(sigmahatE[,1]),ylim=c(-37,1),col=1,type="l",ylab=expression(paste("log(",hat(sigma)[paste(s)],")")),xlab="# Iterations")
curve(0*x,from=0,to=length(x),add=T,lty=2,col=2)

par(mfg=c(2,2))
plot(xE,hmm.hmm.mod_EqVar[[75]][[3]][[2]]$Log[-1],type="l",ylab="Log-Likelihood",xlab="# Iterations",ylim=c(-95,-5))

hmm.hmm.mod_EqVar_adjm
boxplot(t(hmm.hmm.mod_BIC[3,,]),ylab="BIC",xlab="Parent Option")
quartz()
boxplot(t(hmm.hmm.mod_EqVar_BIC[3,,]),ylab="BIC",xlab="Parent Option")

which(hmm.hmm.mod_EqVar_BIC[3,5,] == min(hmm.hmm.mod_EqVar_BIC[3,5,]),arr.ind=T)
hmm.hmm.mod_EqVar[[94]][[3]][[5]]
ViterbiReg(hmm.dat[[3]][,1],cbind(1,as.matrix(hmm.dat[[3]][,RegressionPosibilities[2,]])),hmm.hmm.mod_EqVar[[75]][[3]][[2]]$A,hmm.hmm.mod_EqVar[[75]][[3]][[2]]$B,hmm.hmm.mod_EqVar[[75]][[3]][[2]]$pi)


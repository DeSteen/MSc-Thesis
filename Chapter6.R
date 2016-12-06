install.packages("igraph")
install.packages("R.matlab")

setwd("/Users/dennissteenhuis/Dropbox/Msc Thesis/Rscripts")
library(MASS)
library(xtable)
library(igraph)
library(PRROC)
library(R.matlab)
source("auroc.R")
source("auc-pr.R")
source("BaumWelchReg.R")
source("BwdAlgReg.R")
source("ExtractSNSModels.R")
source("FwdAlgReg.R")
source("gplotyeast.R")
source("MultQSNS.R")
source("SNS.R")
source("StepGraph.R")
source("ViterbiReg.R")

##Use same data as in Grzegorczyk2016
source("yeast_data.R")

yeast.data <- PreProcessed.time

gplotyeast(adjmTrue)

##Homogeneous Model (Control, 1 group)
LinearModel <- list()
LinearModel_score <- LinearModel_BIC <- matrix(rep(0,75),ncol=15)
	for(i in 1:5){
		LinearModel[[i]] <- list()
		for(j in 1:15){
			LinearModel[[i]][[j]] <- lm(as.numeric(yeast.data[[i]][,1])~cbind(1,as.matrix(yeast.data[[i]][,c(RegressionPosibilities[j,])])) -1 )
			LinearModel_BIC[i,j] <- BIC(LinearModel[[i]][[j]])
		}
		LinearModel_score[i,which.min(LinearModel_BIC[i,])] <- 1 
	}

LinearModel_adjm <- matrix(rep(NA,25),ncol=5)
for(i in 1:5){
	LinearModel_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*LinearModel_score[i,],2,sum))
}
auroc(LinearModel_adjm,t(adjmTrue),"Homogeneous model")
auc.pr(LinearModel_adjm,t(adjmTrue),"Homogeneous model")
quartz()
gplotyeast(t(LinearModel_adjm))


LinearModelFull1_adjm <- StepGraph(yeast.data,direc="Rand",type="lm",InitMat = t(adjmTrue))
LinearModelFull2_adjm <- StepGraph(yeast.data,direc="Rand",type="lm",InitMat = LinearModel_adjm)


##Changepoint models
PreciseCPM_BIC <- matrix(rep(0,75),ncol=15)
PreciseCPM_score <- matrix(rep(0,75),ncol=15)
PreciseCPM_adjm <- matrix(rep(NA,25),ncol=5)
for(i in 1:5){
	for(j in 1:15){
		PreciseCPM_BIC[i,j] <- -2*logLik(lm(yeast.data[[i]][,1] ~ as.matrix(yeast.data[[i]][,c(RegressionPosibilities[j,])]),subset=(yeast.data[[i]][,6]=="On"))) -2*logLik(lm(yeast.data[[i]][,1] ~ as.matrix(yeast.data[[i]][,c(RegressionPosibilities[j,])]),subset=(yeast.data[[i]][,6]=="Off"))) + (2+2+2*sum(RegressionPosibilities[j,]!=0))*log(length(yeast.data[[i]][,1]))
	}
	PreciseCPM_score[i,] <- PreciseCPM_BIC[i,] == min(PreciseCPM_BIC[i,])
	PreciseCPM_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*PreciseCPM_score[i,],2,sum))
}
auc.pr(PreciseCPM_adjm,t(adjmTrue),"Precise CPM")
gplotyeast(t(PreciseCPM_adjm))

PreciseCPMFull1_adjm <- StepGraph(yeast.data,direc="Rand",type="lmsub",InitMat = t(adjmTrue))
PreciseCPMFull2_adjm <- StepGraph(yeast.data,direc="Rand",type="lmsub",InitMat = PreciseCPM_adjm)

auc.pr(PreciseCPMFull1_adjm[[1]],t(adjmTrue),"PreciseCPM Full")
gplotyeast(t(PreciseCPMFull1_adjm[[1]]))

##Fixed parents, infering 1 changepoints
ParentCPM1lFix <- SNS(yeast.data[[1]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3,5)])),1,type="fixed",Qout=T)
ParentCPM2lFix <- SNS(yeast.data[[2]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(2,4)])),1,type="fixed",Qout=T)
ParentCPM3lFix <- SNS(yeast.data[[3]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3)])),1,type="fixed",Qout=T)
ParentCPM4lFix <- SNS(yeast.data[[4]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3,4)])),1,type="fixed",Qout=T)
ParentCPM5lFix <- SNS(yeast.data[[5]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(4)])),1,type="fixed",Qout=T)
ParentCPMFulllFix <- MultQSNS(list(ParentCPM1lFix,ParentCPM2lFix,ParentCPM3lFix,ParentCPM4lFix,ParentCPM5lFix),"BIC")

ParentCPM1bFix <- SNS(yeast.data[[1]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3,5)])),1,mincrit="MDL",type="fixed",Qout=T)
ParentCPM2bFix <- SNS(yeast.data[[2]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(2,4)])),1,mincrit="MDL",type="fixed",Qout=T)
ParentCPM3bFix <- SNS(yeast.data[[3]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3)])),1,mincrit="MDL",type="fixed",Qout=T)
ParentCPM4bFix <- SNS(yeast.data[[4]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3,4)])),1,mincrit="MDL",type="fixed",Qout=T)
ParentCPM5bFix <- SNS(yeast.data[[5]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(4)])),1,mincrit="MDL",type="fixed",Qout=T)
ParentCPMFullbFix <- MultQSNS(list(ParentCPM1bFix,ParentCPM2bFix,ParentCPM3bFix,ParentCPM4bFix,ParentCPM5bFix),"MDL")

##Fixed parents, infering best number of changepoints
ParentCPM1b <- SNS(yeast.data[[1]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3,5)])),7,Qout=T)
ParentCPM2b <- SNS(yeast.data[[2]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(2,4)])),7,Qout=T)
ParentCPM3b <- SNS(yeast.data[[3]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3)])),7,Qout=T)
ParentCPM4b <- SNS(yeast.data[[4]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3,4)])),7,Qout=T)
ParentCPM5b <- SNS(yeast.data[[5]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(4)])),7,Qout=T)
ParentCPMFullb <- MultQSNS(list(ParentCPM1b,ParentCPM2b,ParentCPM3b,ParentCPM4b,ParentCPM5b),"BIC")

ParentCPM1a <- SNS(yeast.data[[1]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3,5)])),7,mincrit="AIC",Qout=T)
ParentCPM2a <- SNS(yeast.data[[2]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(2,4)])),7,mincrit="AIC",Qout=T)
ParentCPM3a <- SNS(yeast.data[[3]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3)])),7,mincrit="AIC",Qout=T)
ParentCPM4a <- SNS(yeast.data[[4]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3,4)])),7,mincrit="AIC",Qout=T)
ParentCPM5a <- SNS(yeast.data[[5]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(4)])),7,mincrit="AIC",Qout=T)
ParentCPMFulla <- MultQSNS(list(ParentCPM1a,ParentCPM2a,ParentCPM3a,ParentCPM4a,ParentCPM5a),"AIC")

ParentCPM1m <- SNS(yeast.data[[1]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3,5)])),7,mincrit="MDL",Qout=T)
ParentCPM2m <- SNS(yeast.data[[2]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(2,4)])),7,mincrit="MDL",Qout=T)
ParentCPM3m <- SNS(yeast.data[[3]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3)])),7,mincrit="MDL",Qout=T)
ParentCPM4m <- SNS(yeast.data[[4]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(3,4)])),7,mincrit="MDL",Qout=T)
ParentCPM5m <- SNS(yeast.data[[5]][,1],cbind(rep(1,dim(yeast.data[[1]])[1]),as.matrix(yeast.data[[1]][,c(4)])),7,mincrit="MDL",Qout=T)
ParentCPMFullm <- MultQSNS(list(ParentCPM1m,ParentCPM2m,ParentCPM3m,ParentCPM4m,ParentCPM5m),"MDL")


# ## CPM Modeling, use load("CPM.Rda")
# CPMbic <- CPMaic <- CPMmdl <- CPMbic1 <- CPMmdl1 <- list()
# for(i in 1:5){
	 # CPMaic[[i]] <- CPMbic[[i]] <- CPMmdl[[i]] <- CPMbic1[[i]] <- CPMmdl1[[i]] <- list()
	# for(j in 1:15){
		# CPMbic[[i]][[j]] <- SNS(as.numeric(yeast.data[[i]][,1]),cbind(1,as.matrix(yeast.data[[i]][,c(RegressionPosibilities[j,])])),7,type="min")
		# CPMaic[[i]][[j]] <- SNS(as.numeric(yeast.data[[i]][,1]),cbind(1,as.matrix(yeast.data[[i]][,c(RegressionPosibilities[j,])])),7,"AIC",type="min")
		# CPMmdl[[i]][[j]] <- SNS(as.numeric(yeast.data[[i]][,1]),cbind(1,as.matrix(yeast.data[[i]][,c(RegressionPosibilities[j,])])),7,"MDL",type="min")
		# CPMbic1[[i]][[j]] <- SNS(as.numeric(yeast.data[[i]][,1]),cbind(1,as.matrix(yeast.data[[i]][,c(RegressionPosibilities[j,])])),1,type="fixed")
		# CPMmdl1[[i]][[j]] <- SNS(as.numeric(yeast.data[[i]][,1]),cbind(1,as.matrix(yeast.data[[i]][,c(RegressionPosibilities[j,])])),1,"MDL",type="fixed")
		# cat("\r Node ",i," Posibility ",j)
	# }
# }
# save(CPMaic,CPMbic,CPMmdl,CPMbic1,CPMmdl1,file="CPM.Rda")
load("CPM.Rda")

CPMaic_BIC <- CPMbic_BIC <- CPMmdl_BIC <- CPMbic1_BIC <- CPMmdl1_BIC <- matrix(rep(0,75),ncol=15)
CPMaic_score <- CPMbic_score <- CPMmdl_score <- CPMbic1_score <- CPMmdl1_score <- matrix(rep(0,75),ncol=15)
CPMaic_adjm <- CPMbic_adjm <- CPMmdl_adjm <- CPMbic1_adjm <- CPMmdl1_adjm <- matrix(rep(NA,25),ncol=5)
for(i in 1:5){
	for(j in 1:15){
		CPMaic_BIC[i,j] <- CPMaic[[i]][[j]]$BIC
		CPMbic_BIC[i,j] <- CPMbic[[i]][[j]]$BIC
		CPMmdl_BIC[i,j] <- CPMmdl[[i]][[j]]$BIC
		CPMbic1_BIC[i,j] <- CPMbic1[[i]][[j]]$BIC
		CPMmdl1_BIC[i,j] <- CPMmdl1[[i]][[j]]$BIC
	}
	CPMaic_score[i,] <- CPMaic_BIC[i,] == min(CPMaic_BIC[i,])
	CPMbic_score[i,] <- CPMbic_BIC[i,] == min(CPMbic_BIC[i,])
	CPMmdl_score[i,] <- CPMmdl_BIC[i,] == min(CPMmdl_BIC[i,])
	CPMbic1_score[i,] <- CPMbic1_BIC[i,] == min(CPMbic1_BIC[i,])
	CPMmdl1_score[i,] <- CPMmdl1_BIC[i,] == min(CPMmdl1_BIC[i,])
	
	CPMaic_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*CPMaic_score[i,],2,sum))
	CPMbic_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*CPMbic_score[i,],2,sum))
	CPMmdl_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*CPMmdl_score[i,],2,sum))
	CPMbic1_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*CPMbic1_score[i,],2,sum))
	CPMmdl1_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*CPMmdl1_score[i,],2,sum))
}

##Both minimizing LL and BIC give twice the same graph.
par(mfrow=c(2,2))
gplotyeast(t(CPMbic1_adjm))
auc.pr(CPMbic1_adjm,t(adjmTrue),"CPM, 1 changepoint, LL")
gplotyeast(t(CPMmdl1_adjm))
auroc(CPMbic1_adjm,t(adjmTrue),"CPM, 1 changepoint, BIC")

CPMbic1.full <- StepGraph(yeast.data,type="SNS",InitMat = matrix(c(NA,F,T,F,T,T,NA,F,F,F,T,F,NA,T,T,F,F,T,NA,T,F,T,T,T,NA),byrow=T,ncol=5),SNStype="fixed",Ncp=1)
CPMmdl1.full <- StepGraph(yeast.data,type="SNS",InitMat = matrix(c(NA,F,T,F,T,T,NA,F,F,F,T,F,NA,T,T,F,F,T,NA,T,F,T,T,T,NA),byrow=T,ncol=5),mincrit="MDL",SNStype="fixed",Ncp=1)

par(mfrow=c(2,2))
gplotyeast(t(CPMbic1.full$AdjencyMatrix))
auc.pr(CPMbic1.full$AdjencyMatrix,t(adjmTrue),"CPM, 1 changepoint, LL, Full BIC")
gplotyeast(t(CPMmdl1.full$AdjencyMatrix))
auroc(CPMmdl1.full$AdjencyMatrix,t(adjmTrue),"CPM, 1 changepoint, BIC, Full BIC")


CPMbic1.full.glob <- StepGraph(yeast.data,type="SNSFull",InitMat = CPMbic1_adjm,SNStype="fixed",Ncp=1)
CPMmdl1.full.glob <- StepGraph(yeast.data,type="SNSFull",InitMat = CPMbic1_adjm,mincrit="MDL",SNStype="fixed",Ncp=1)

par(mfrow=c(2,2))
gplotyeast(t(CPMbic1.full.glob$AdjencyMatrix))
auc.pr(CPMbic1.full.glob$AdjencyMatrix,t(adjmTrue),"CPM, 1 changepoint, LL, Full BIC")
gplotyeast(t(CPMmdl1.full.glob$AdjencyMatrix))
auc.pr(CPMmdl1.full.glob$AdjencyMatrix,t(adjmTrue),"CPM, 1 changepoint, BIC, Full BIC")

CPMaic_score
CPMaic[[1]][[1]][[1]]; CPMaic[[2]][[9]][[1]]; CPMaic[[3]][[5]][[1]]; CPMaic[[4]][[2]][[1]]; CPMaic[[5]][[14]][[1]]

CPMbic_score
CPMbic[[1]][[4]][[1]]; CPMbic[[2]][[9]][[1]]; CPMbic[[3]][[5]][[1]]; CPMbic[[4]][[2]][[1]]; CPMbic[[5]][[14]][[1]]

CPMmdl_score
CPMmdl[[1]][[1]][[1]]; CPMmdl[[2]][[9]][[1]]; CPMmdl[[3]][[5]][[1]]; CPMmdl[[4]][[2]][[1]]; CPMmdl[[5]][[14]][[1]]

quartz()
par(mfrow=c(2,2))
auroc(CPMaic_adjm,t(adjmTrue),"SNS, #CPs chosen by AIC")
auc.pr(CPMbic_adjm,t(adjmTrue),"SNS, #CPs chosen by BIC")
auc.pr(CPMmdl_adjm,t(adjmTrue),"SNS, #CPs chosen by MDL")

CPMaic.full <- StepGraph(yeast.data,type="SNS",InitMat = CPMaic_adjm,mincrit="AIC",Ncp=8,plotgraph=T)
CPMbic.full <- StepGraph(yeast.data,type="SNS",InitMat = t(adjmTrue),mincrit="BIC",Ncp=8,plotgraph=T)
CPMmdl.full <- StepGraph(yeast.data,type="SNS",InitMat = CPMaic_adjm,mincrit="MDL",Ncp=8,plotgraph=T)

par(mfrow=c(2,2))
auroc(CPMaic.full[[1]],t(adjmTrue),"SNS, #CPs chosen by AIC")
auc.pr(CPMbic.full[[1]],t(adjmTrue),"SNS, #CPs chosen by BIC")
auc.pr(CPMmdl.full[[1]],t(adjmTrue),"SNS, #CPs chosen by MDL")
gplotyeast(CPMbic.full[[1]])

CPMaic.full
CPMbic.full
CPMmdl.full


quartz()
par(mfrow=c(2,2))
gplotyeast(t(CPMaic_adjm))
gplotyeast(t(CPMbic_adjm))
gplotyeast(t(CPMmdl_adjm))


SNSFullAIC <- StepGraph(yeast.data,type="SNSFull",mincrit="AIC")
SNSFullBIC <- StepGraph(yeast.data,type="SNSFull")
SNSFullMDL <- StepGraph(yeast.data,type="SNSFull",mincrit="MDL")
for(i in 1:50){
	SNSFullAICnew <- StepGraph(yeast.data,type="SNSFull",mincrit="AIC")
	SNSFullBICnew <- StepGraph(yeast.data,type="SNSFull")
	SNSFullMDLnew <- StepGraph(yeast.data,type="SNSFull",mincrit="MDL")
	if(SNSFullAICnew[[2]] < SNSFullAIC[[2]]){SNSFullAIC <- SNSFullAICnew}
	if(SNSFullBICnew[[2]] < SNSFullBIC[[2]]){SNSFullBIC <- SNSFullBICnew}
	if(SNSFullMDLnew[[2]] < SNSFullMDL[[2]]){SNSFullMDL <- SNSFullMDLnew}
	cat("\r iteration ",i)
}
save(SNSFullAIC,SNSFullBIC,SNSFullMDL,file="SNSFulls.Rda")
load("SNSFulls.Rda")

SNSFullAIC
SNSFullBIC
SNSFullMDL

par(mfrow=c(3,2))
gplotyeast(t(SNSFullAIC[[1]])); au(SNSFullAIC[[1]],t(adjmTrue),"SNSFullAIC")
gplotyeast(t(SNSFullBIC[[1]])); auc.pr(SNSFullBIC[[1]],t(adjmTrue),"SNSFullBIC")
gplotyeast(t(SNSFullMDL[[1]])); auc.pr(SNSFullMDL[[1]],t(adjmTrue),"SNSFullMDL")

##HMM's
HMM <- HMM_EqVar <- HMM_LR <- HMM_LR_EqVar <- list()
for(n in 1:200){
	HMM[[n]] <- HMM_EqVar[[n]] <- HMM_LR[[n]] <- HMM_LR_EqVar[[n]] <- list()
	for(i in 1:5){
		HMM[[n]][[i]] <- HMM_EqVar[[n]][[i]] <- HMM_LR[[n]][[i]] <- HMM_LR_EqVar[[n]][[i]] <- list()
		for(j in 1:15){
			cat("\r node",i,"posibility",formatC(j,width=2),"in iteration",formatC(n,width=3))
			HMM[[n]][[i]][[j]] <- BaumWelchReg(yeast.data[[i]][,1],cbind(1,as.matrix(yeast.data[[i]][,RegressionPosibilities[j,]])),2,-1500,1)
			HMM_EqVar[[n]][[i]][[j]] <- BaumWelchReg(yeast.data[[i]][,1],cbind(1,as.matrix(yeast.data[[i]][,RegressionPosibilities[j,]])),2,-1500,1,eq.var=T)
			HMM_LR[[n]][[i]][[j]] <- BaumWelchReg(yeast.data[[i]][,1],cbind(1,as.matrix(yeast.data[[i]][,RegressionPosibilities[j,]])),2,-1500,1,LR=T)
			HMM_LR_EqVar[[n]][[i]][[j]] <- BaumWelchReg(yeast.data[[i]][,1],cbind(1,as.matrix(yeast.data[[i]][,RegressionPosibilities[j,]])),2,-1500,1,LR=T,eq.var=T)

		} 
	}
}


HMM_BIC <- HMM_EqVar_BIC <- HMM_LR_BIC <- HMM_LR_EqVar_BIC <- array(rep(0,75*200),dim=c(5,15,200)) 
HMM_score <- HMM_EqVar_score <- HMM_LR_score <- HMM_LR_EqVar_score <- matrix(rep(0,75),ncol=15)
HMM_score_outlier <- HMM_EqVar_score_outlier <- HMM_LR_score_outlier <- HMM_LR_EqVar_score_outlier <- matrix(rep(0,75),ncol=15)
HMM_adjm <- HMM_LR_adjm <- HMM_EqVar_adjm <- HMM_LR_EqVar_adjm <- matrix(rep(NA,25),ncol=5)
HMM_adjm_outlier <- HMM_LR_adjm_outlier <- HMM_EqVar_adjm_outlier <- HMM_LR_EqVar_adjm_outlier <- matrix(rep(NA,25),ncol=5)
for(i in 1:5){
	for(n in 1:200){
		for(j in 1:15){
			HMM_BIC[i,j,n] <- HMM[[n]][[i]][[j]]$BIC
			HMM_EqVar_BIC[i,j,n] <- HMM_EqVar[[n]][[i]][[j]]$BIC
			HMM_LR_BIC[i,j,n] <- HMM_LR[[n]][[i]][[j]]$BIC
			HMM_LR_EqVar_BIC[i,j,n] <- HMM_LR_EqVar[[n]][[i]][[j]]$BIC
		}
	}
	HMM_score[i,] <- apply(HMM_BIC[i,,] == min(HMM_BIC[i,,]),1,sum)		
	HMM_EqVar_score[i,] <- apply(HMM_EqVar_BIC[i,,] == min(HMM_EqVar_BIC[i,,]),1,sum)
	HMM_LR_score[i,] <- apply(HMM_LR_BIC[i,,] == min(HMM_LR_BIC[i,,]),1,sum)
	HMM_LR_EqVar_score[i,] <- apply(HMM_LR_EqVar_BIC[i,,] == min(HMM_LR_EqVar_BIC[i,,]),1,sum)
	
	HMM_score_outlier[i,] <- apply(apply(HMM_BIC[i,,],1,RemoveOutlier,k=1) == min(apply(HMM_BIC[i,,],1,RemoveOutlier,k=1),na.rm=T),2,sum,na.rm=T)		
	HMM_EqVar_score_outlier[i,] <- apply(apply(HMM_EqVar_BIC[i,,],1,RemoveOutlier) == min(apply(HMM_EqVar_BIC[i,,],1,RemoveOutlier),na.rm=T),2,sum,na.rm=T)	
	HMM_LR_score_outlier[i,] <- apply(apply(HMM_LR_BIC[i,,],1,RemoveOutlier,k=1) == min(apply(HMM_LR_BIC[i,,],1,RemoveOutlier,k=1),na.rm=T),2,sum,na.rm=T)	
	HMM_LR_EqVar_score_outlier[i,] <- apply(apply(HMM_LR_EqVar_BIC[i,,],1,RemoveOutlier) == min(apply(HMM_LR_EqVar_BIC[i,,],1,RemoveOutlier),na.rm=T),2,sum,na.rm=T)	
	
	HMM_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*(HMM_score[i,]!=0),2,sum))
	HMM_LR_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*(HMM_LR_score[i,]!=0),2,sum))
	HMM_EqVar_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*(HMM_EqVar_score[i,]!=0),2,sum))
	HMM_LR_EqVar_adjm[i,-i] <- (apply((RegressionPosibilities!=0)*(HMM_LR_EqVar_score[i,]!=0),2,sum))
	
	HMM_adjm_outlier[i,-i] <- (apply((RegressionPosibilities!=0)*(HMM_score_outlier[i,]!=0),2,sum))
	HMM_LR_adjm_outlier[i,-i] <- (apply((RegressionPosibilities!=0)*(HMM_LR_score_outlier[i,]!=0),2,sum))
	HMM_EqVar_adjm_outlier[i,-i] <- (apply((RegressionPosibilities!=0)*(HMM_EqVar_score_outlier[i,]!=0),2,sum))
	HMM_LR_EqVar_adjm_outlier[i,-i] <- (apply((RegressionPosibilities!=0)*(HMM_LR_EqVar_score_outlier[i,]!=0),2,sum))
	
}

par(mfrow=c(2,2))
auc.pr(HMM_adjm,t(adjmTrue),"HMM")
gplotyeast(t(HMM_adjm))
auroc(HMM_adjm_outlier,t(adjmTrue),"HMM")
gplotyeast(t(HMM_adjm_outlier))

par(mfrow=c(2,2))
boxplot(t(HMM_BIC[1,,]))
quartz()
boxplot(apply(HMM_BIC[1,,],1,RemoveOutlier,k=-0.27))

par(mfrow=c(2,2))
auc.pr(HMM_LR_adjm,t(adjmTrue),"Left-Right HMM")
gplotyeast(t(HMM_LR_adjm))
auroc(HMM_LR_adjm_outlier,t(adjmTrue),"Left-Right HMM")
gplotyeast(t(HMM_LR_adjm_outlier))

par(mfrow=c(2,2))
auroc(HMM_EqVar_adjm,t(adjmTrue),"HMM with equal variances")
auc.pr(HMM_EqVar_adjm,t(adjmTrue),"HMM with equal variances")
gplotyeast(t(HMM_EqVar_adjm))
auroc(HMM_EqVar_adjm_outlier,t(adjmTrue),"HMM with equal variances")
gplotyeast(t(HMM_EqVar_adjm_outlier))

par(mfrow=c(2,2))
auc.pr(HMM_LR_EqVar_adjm,t(adjmTrue),"Left-Right HMM with equal variances")
gplotyeast(t(HMM_LR_EqVar_adjm))
auroc(HMM_LR_EqVar_adjm_outlier,t(adjmTrue),"Left-Right HMM with equal variances")
gplotyeast(t(HMM_LR_EqVar_adjm_outlier))

par(mfrow=c(2,2))
auroc(HMM_adjm_outlier,t(adjmTrue),"HMM")
auroc(HMM_LR_adjm_outlier,t(adjmTrue),"Left-Right HMM")
auroc(HMM_EqVar_adjm_outlier,t(adjmTrue),"HMM with equal variances")
auroc(HMM_LR_EqVar_adjm_outlier,t(adjmTrue),"Left-Right HMM with equal variances")

save(HMM,file="HMM.Rda")
save(HMM_LR,file="HMM_LR.Rda")
save(HMM_EqVar,file="HMM_EqVar.Rda")
save(HMM_LR_EqVar,file="HMM_LR_EqVar.Rda")

par(mfrow=c(2,2))
a <- boxplot(t(HMM_EqVar_BIC[5,,]))
plot(apply(HMM_BIC[1,,],1,mean))
plot(apply(HMM_BIC[1,,],1,median))

HMM_seq1 <- matrix(0,nrow=5,ncol=33); colnames(HMM_seq1) <- c(2:15,17:35); rownames(HMM_seq1) <- c("CBF1","GAL4","SWI5","GAL80","ASH1")
for(i in 1:5){
opt <- which(HMM_LR_EqVar_BIC[i,,] == min(HMM_LR_EqVar_BIC[i,,]),arr.ind=T)
HMM_seq1[i,] <- ViterbiReg(as.numeric(yeast.data[[i]][,1]),cbind(1,as.matrix(yeast.data[[i]][,c(RegressionPosibilities[opt[1,1],])])),HMM_LR_EqVar[[opt[1,2]]][[i]][[opt[1,1]]]$A,HMM_LR_EqVar[[opt[1,2]]][[i]][[opt[1,1]]]$B,HMM_LR_EqVar[[opt[1,2]]][[i]][[opt[1,1]]]$pi)
}
xtable(t(HMM_seq1),digits=0)

HMM_seq2 <- matrix(0,nrow=5,ncol=33); colnames(HMM_seq2) <- c(2:15,17:35); rownames(HMM_seq2) <- c("CBF1","GAL4","SWI5","GAL80","ASH1")
for(i in 1:5){
opt <- which(HMM_EqVar_BIC[i,,] == min(HMM_EqVar_BIC[i,,]),arr.ind=T)
HMM_seq2[i,] <- ViterbiReg(as.numeric(yeast.data[[i]][,1]),cbind(1,as.matrix(yeast.data[[i]][,c(RegressionPosibilities[opt[1,1],])])),HMM_EqVar[[opt[1,2]]][[i]][[opt[1,1]]]$A,HMM_EqVar[[opt[1,2]]][[i]][[opt[1,1]]]$B,HMM_EqVar[[opt[1,2]]][[i]][[opt[1,1]]]$pi)
}
xtable(t(HMM_seq2),digits=0)

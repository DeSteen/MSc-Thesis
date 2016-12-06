StepGraph <- function(dat,direc="Rand",type="lm",InitMat=NaN,mincrit="BIC",SNStype="min",Ncp = 7,plotgraph=F){
	Npoints <- length(dat[[1]][,1])
	Nnodes <- length(dat)
Models <-list()
InitBIC <- 0
Nseg <- NaN
if(SNStype == "fixed"){Nseg <- Ncp + 1}


if(type == "lm"){
	modtype <- function(y,X){lm(y~X-1)}
	LLtype <- function(Mod){logLik(Mod)}
	npar <- function(Mod){1 + length(coef(Mod))}
}
if(type == "lmsub"){
	modtype <- function(y,X){ModeList <- list()
		for(i in 1:length(levels(dat[[1]][,(Nnodes+1)]))){
			ModeList[[i]] <- lm(y~X-1,subset=(dat[[1]][,(Nnodes+1)] == levels(dat[[1]][,(Nnodes+1)])[i]))
		}
		return(ModeList)
		}
	LLtype <- function(Mod){sum(sapply(Mod,logLik))}
	npar <- function(Mod){length(Mod) * (1 + length(coef(Mod[[1]])))}
}
if(type == "SNS"){
	modtype <- function(y,X){SNS(y,X,Ncp,mincrit,SNStype)}
	LLtype <- function(Mod){Mod$LogLikelihood}
	npar <- function(Mod){(1+length(coef(Mod$Models[[1]])))*(sum(Mod$Tau != 0)+1) + sum(Mod$Tau != 0)}
}
if(type == "SNSFull"){
	modtype <- function(y,X){SNS(y,X,Ncp,mincrit,SNStype,Qout=T)}
	LLtype <- function(Mods){MultQSNS(Mods,mincrit,Nseg)}
	npar <- function(Mod){(1+length(coef(Mod$Models[[1]])))*(sum(Mod$Tau != 0)+1) + sum(Mod$Tau != 0)}
}

if(direc=="Fwd" || direc=="Bwd"){BestModelMat <- matrix(rep(direc=="Bwd",Nnodes^2),ncol=Nnodes,nrow=Nnodes)}
if(direc=="Rand" && is.nan(InitMat) == T){BestModelMat <- matrix(as.logical(rbinom(Nnodes^2,1,0.5)),ncol=Nnodes,nrow=Nnodes)}
if(direc=="Rand" && is.nan(InitMat) == F){BestModelMat <- matrix(F,Nnodes,Nnodes)
	for(i in 1:Nnodes){
		BestModelMat[i,-1] <- as.logical(InitMat[i,-i])
	}
}

BestModelMat[,1] <- F

if(type != "SNSFull"){
	if(direc == "Fwd"){
		for(i in 1:Nnodes){
			Models[[i]] <- modtype(as.numeric(dat[[i]][,1]),rep(1,Npoints))
			InitBIC <- InitBIC - 2*LLtype(Models[[i]]) + npar(Models[[i]])*log(Nnodes*Npoints)
		} 
	}

	if(direc == "Bwd"){
		for(i in 1:Nnodes){
			Models[[i]] <- modtype(as.numeric(dat[[i]][,1]),cbind(rep(1,Npoints),as.matrix(dat[[i]][,2:Nnodes])))
			InitBIC <- InitBIC - 2*LLtype(Models[[i]]) + npar(Models[[i]])*log(Nnodes*Npoints)
		}	
	}

	if(direc == "Rand"){
		for(i in 1:Nnodes){
			Models[[i]] <- modtype(as.numeric(dat[[i]][,1]),cbind(rep(1,Npoints),as.matrix(dat[[i]][,BestModelMat[i,]])))
			InitBIC <- InitBIC - 2*LLtype(Models[[i]]) + npar(Models[[i]])*log(Nnodes*Npoints)
		}
	}
}

if(type == "SNSFull"){
	if(direc == "Fwd"){
		for(i in 1:Nnodes){
			Models[[i]] <- modtype(as.numeric(dat[[i]][,1]),rep(1,Npoints))
		} 
		InitBIC <- MultQSNS(Models,mincrit)$bic
	}

	if(direc == "Bwd"){
		for(i in 1:Nnodes){
			Models[[i]] <- modtype(as.numeric(dat[[i]][,1]),cbind(rep(1,Npoints),as.matrix(dat[[i]][,2:Nnodes])))
		}	
		InitBIC <- MultQSNS(Models,mincrit)$bic
	}

	if(direc == "Rand"){
		for(i in 1:Nnodes){
			Models[[i]] <- modtype(as.numeric(dat[[i]][,1]),cbind(rep(1,Npoints),as.matrix(dat[[i]][,BestModelMat[i,]])))
		}
		InitBIC <- MultQSNS(Models,mincrit)$bic
	}
}

print(InitBIC)
BestBIC <- InitBIC + 1
NewBestBIC <- InitBIC
counter <- 0
BestNewModelMat <- BestModelMat
while(NewBestBIC < BestBIC){
	BestBIC <- NewBestBIC
	BestOpt <- "Nothing"
	#print(BestModelMat)
	for(i in 1:Nnodes){
		for(j in 2:Nnodes){
			#print(c(i,j))
				NewModelMat <- BestModelMat
				NewModels <- Models
				NewModelMat[i,j] <- as.logical(1-NewModelMat[i,j])
				if(NewModelMat[i,j] == T){deladd <- "adding"}
				if(NewModelMat[i,j] == F){deladd <- "deleting"}
				if(i < j){to <- j}
				if(i >= j){to <- j - 1}
				NewMod <- modtype(as.numeric(dat[[i]][,1]),cbind(rep(1,Npoints),as.matrix(dat[[i]][,NewModelMat[i,]])))
				#print(c(i,j))
				if(type != "SNSFull"){
					NewBIC <- BestBIC + 2*LLtype(Models[[i]]) - 2*LLtype(NewMod) + log(Nnodes*Npoints)*(npar(NewMod) - npar(Models[[i]]))
				}
				if(type == "SNSFull"){
					NewModels[[Nnodes + 1]] <- NewMod
					#print(class(NewModels))
					MQ <- MultQSNS(NewModels[-i],mincrit,Nseg)
					NewBIC <- MQ$bic
				}
				#print(NewBIC)
				if(NewBIC < NewBestBIC){NewBestBIC <- NewBIC; BestNewMod <- NewMod; BestNewModelMat <- NewModelMat; BestOpt <- paste(deladd,"the edge from",to,"to",i) }
				else{NewBestBIC <- NewBestBIC; }
		}
	}
	counter <- counter + 1
	print(paste("Iteration",counter,"with BIC =",NewBestBIC,BestOpt))
	NodeChange <- which(BestNewModelMat != BestModelMat, arr.ind=T)[1]
	if(is.na(NodeChange) == T){
		adjm <- matrix(NA,ncol=Nnodes,nrow=Nnodes)
		for(i in 1:Nnodes){
			adjm[i,-i] <- BestModelMat[i,2:Nnodes]
		}
		if(type == "lm"){if(plotgraph==T){gplotyeast(t(adjm))}; return(list("AdjencyMatrix"=adjm,"BIC-Value"=BestBIC[1],"Models"=Models))}
		if(type == "lmsub"){if(plotgraph==T){gplotyeast(t(adjm))}; return(list("AdjencyMatrix"=adjm,"BIC-Value"=BestBIC[1],"Models"=Models))}
		if(type == "SNS"){if(plotgraph==T){gplotyeast(t(adjm))}; return(list("AdjencyMatrix"=adjm,"BIC-Value"=BestBIC,"Models"=ExtractSNSModels(Models)))}
		if(type == "SNSFull"){
			MQ <- MultQSNS(Models,mincrit,Nseg)
			if(plotgraph==T){gplotyeast(t(adjm))}
			return(list("AdjencyMatrix"=adjm,"BIC-Value"=BestBIC,"Tau"=MQ$Tau))
		}
	}
	else{
		Models[[NodeChange]] <- BestNewMod
		BestModelMat <- BestNewModelMat
		if(plotgraph == T){
			adjm <- matrix(NA,ncol=Nnodes,nrow=Nnodes)
			for(i in 1:Nnodes){
				adjm[i,-i] <- BestModelMat[i,2:Nnodes]}
			gplotyeast(t(adjm))
		}
	}
	
	
}

}
InitVal <- function(Y,X,M,minLL,sigma,LR=F){
	X <- as.matrix(X)
	N <- length(Y)
	d <- dim(X)[2]
	## Setting initial values
	phat0 <- runif(M,0,1)
	phat0 <- phat0/sum(phat0)
	Ahat0 <- matrix(runif(M*M,0,1),nrow=M)
	if(LR==T){Ahat0[lower.tri(Ahat0)] <- 0; phat0 <- rep(0,M); phat0[1] <- 1}
	for(m in 1:M){
		Ahat0[m,] <- Ahat0[m,]/sum(Ahat0[m,])
	}
	Bhat0 <- matrix(rep(NA,M*(d+1)),nrow=M)
	
	LL0full <- -Inf
	while(LL0full < minLL){
		if(LR==T){
			GroupsY <- split(Y,ceiling(seq_along(Y)/(ceiling(N/M))))
			GroupsX <- list()
			s <- 1
			for(x in 1:length(GroupsY)){
				GroupsX[[x]] <- as.matrix(X[s:(s+length(GroupsY[[x]])-1),])
				s <- s + length(GroupsY[[x]])
			}
			for(m in 1:M){
				GroupM <- sample(1:length(GroupsY[[m]]),max(d+1,floor(length(GroupsY[[m]])/2)))
				Bhat0[m,1:d] <- coef(lm(GroupsY[[m]][GroupM] ~ GroupsX[[m]][GroupM,] - 1))
			}
		}
		else{ Bhat0[,1:d] <- rnorm(d*M,0,1) }
		
	Bhat0[,d+1] <- rep(sigma,M)
	
	
		
	alpha <- array(data=NA,dim=c(N,M+1))
	beta <- array(data=NA,dim=c(N,M))
	
	alpha <- FwdAlgReg(Y,X,Ahat0,Bhat0,phat0)
	#print(alpha)
	
	LL0full <- -sum(log(alpha[,M+1]))
	#print(LL0full)
	if(is.nan(LL0full)==T || is.na(LL0full)==T){LL0full <- -Inf}
	#print(LLoldfull)
	}## End initialisation loop
	return(list(p=phat0,A=Ahat0,B=Bhat0,LL=LL0full))
}


BaumWelchReg <- function(Y,X,M,minLL,sigma,LR=F,eq.var=F,coef.traj=F){
	Nitt <- ConvItt <- 0
	X <- as.matrix(X)
	N <- length(Y)
	d <- dim(X)[2]
	Models <- list()
	LLoldfull <- c()
	CoefTraj <- list()
	CoefTraj$A <- CoefTraj$B <- CoefTraj$pi <- list()
	
	ABp <- InitVal(Y,X,M,minLL,sigma,LR)
	Ahat0 <- CoefTraj$A[[1]] <- ABp$A
	Bhat0 <- CoefTraj$B[[1]] <- ABp$B
	phat0 <- CoefTraj$pi[[1]] <- ABp$p
	LL0full <- LLoldfull[Nitt+1] <- ABp$LL
	
	alpha <- array(data=NA,dim=c(N,M+1))
	beta <- array(data=NA,dim=c(N,M))
	
	LLerr <- 1
	##Starting EM-algorithm
	while(LLerr > 1e-8){
	
	alpha <- FwdAlgReg(Y,X,Ahat0,Bhat0,phat0)
	beta <- BwdAlgReg(Y,X,Ahat0,Bhat0,phat0,alpha[,M+1])
	
	
	##Creating gamma array
	gamma <- array(data=0,dim=c(N,M))
	for(n in 1:N){
		for(i in 1:M){
				gamma[n,i] <- alpha[n,i]*beta[n,i]
			}
		  gamma[n,] <- gamma[n,]/sum(gamma[n,])
	}
	
	##Creating zeta array
	zeta <- array(data=NA,dim=c(N,M,M))
	for(n in 1:(N-1)){
		for(i in 1:M){
			for(j in 1:M){
				zeta[n,i,j] <- alpha[n,i]*Ahat0[i,j]*dnorm(Y[n+1],Bhat0[j,1:d]%*%X[n+1,],Bhat0[j,d+1])*beta[n+1,j]
			}
		}
	zeta[n,,] <- zeta[n,,]/sum(zeta[n,,])
	}
	
	
	##Estimating new parameters
	#Estimating pi
	phat1 <- gamma[1,]
	
	#Estimating A
	Ahat1 <- matrix(rep(0,M*M),nrow=M)
	for(i in 1:M){
		for(j in 1:M){
			nominator <- 0
			denominator <- 0
				for(n in 1:(N-1)){
					nominator <- nominator + zeta[n,i,j]
					denominator <- denominator  + gamma[n,i]
				}
			  Ahat1[i,j] <- nominator/denominator 
		}
	}
	
	#Estimating beta
	Bhat1 <- array(data=0,c(M,d+1))
	for(m in 1:M){
				GAMMA <- diag(gamma[,m])
				Bhat1[m,1:d] <- ginv(t(X)%*%GAMMA%*%X)%*%t(X)%*%GAMMA%*%Y
				# Models[[m]] <- lm(Y ~ X - 1,weights=gamma[,m])
				# Bhat1[m,1:d] <- coef(Models[[m]])
        
	}
	
	#Estimating sigma
	if(eq.var == F){
		for(m in 1:M){
			nominator <- denominator <- 0
			for(n in 1:N){
				nominator <- nominator + gamma[n,m]*(Y[n]-Bhat1[m,1:d]%*%X[n,])^2
			}
			denominator <- sum(gamma[,m])
			  Bhat1[m,d+1] <- sqrt(nominator/denominator)
		}
	}
	if(eq.var == T){
		nominator <- 0
		for(m in 1:M){
			for(n in 1:N){
				nominator <- nominator + gamma[n,m]*(Y[n]-Bhat1[m,1:d]%*%X[n,])^2
			}
		}
		Bhat1[,d+1] <- sqrt(rep(nominator/N,M))
	}
	
	##Calculating full Log-Likelihood
	
	LLnewfull <- 0
		LLnewfull <- -sum(log(FwdAlgReg(Y,X,Ahat1,Bhat1,phat1)[,M+1]))
	# print(paste("Log Likelihood is:",LLnewfull,"in interation",Nitt))
	if(is.nan(LLnewfull) == T){
		ABp <- InitVal(Y,X,M,minLL,sigma,LR)
		Ahat1 <- ABp$A
		Bhat1 <- ABp$B
		phat1 <- ABp$p
		Nitt <- 0
		LLoldfull <- c()
		LL0full <- LLoldfull[Nitt+1] <- ABp$LL 
		LLerr <- 1
		ConvItt <- ConvItt + 1
	 }
	if(Nitt > 100 && LLnewfull < LLoldfull[Nitt+1]){
		ABp <- InitVal(Y,X,M,minLL,sigma,LR)
		Ahat1 <- ABp$A
		Bhat1 <- ABp$B
		phat1 <- ABp$p
		Nitt <- 0
		LLoldfull <- c()
		LL0full <- LLoldfull[Nitt+1] <- ABp$LL 
		LLerr <- 1
		ConvItt <- ConvItt + 1
	 }
	
	
	phat0 <- CoefTraj$pi[[Nitt+1]] <- phat1
	Ahat0 <- CoefTraj$A[[Nitt+1]] <- Ahat1
	Bhat0 <- CoefTraj$B[[Nitt+1]] <- Bhat1
	
	LLerr <- abs((LLnewfull - LLoldfull[Nitt+1])/(LLnewfull-LL0full))
	# print(paste("LL error is:",LLerr))
	if(is.nan(LLerr)==T || is.na(LLerr)==T){ABp <- InitVal(Y,X,M,minLL,sigma,LR)
		Ahat0 <- CoefTraj$A[[1]] <- ABp$A
		Bhat0 <- CoefTraj$B[[1]] <- ABp$B
		phat0 <- CoefTraj$pi[[1]] <- ABp$p
		Nitt <- 0
		LLoldfull <- c()
		LL0full <- LLoldfull[Nitt+1] <- ABp$LL 
		LLerr <- 1
		ConvItt <- ConvItt + 1
		}
		
	if(ConvItt > 50){
		bic <- Inf
		break
	}
	

	LLoldfull[Nitt+2] <- LLnewfull
	
	Nitt <- Nitt + 1
	bic <- -2*LLoldfull[Nitt+1] + (M*M - M - (sum(lower.tri(Ahat0)))*(LR == T) + (M - 1)*(LR == F) + d*M + M*(eq.var == F) + as.numeric(eq.var == T))*log(N)
	}## End while loop
	out <- list("pi"=phat0,"A"=Ahat0,"B"=Bhat0,"LogLikelihood"=LLoldfull,"BIC"=bic)
	if(coef.traj==T){out <- list("pi"=phat0,"A"=Ahat0,"B"=Bhat0,"LogLikelihood"=LLoldfull,"BIC"=bic,"Coef.Traj"=CoefTraj)}
	return(out)
}


SNS <- function(Y,X,M,mincrit="BIC",type="min",Qout=F){
	if(mincrit != "BIC" && mincrit != "AIC" && mincrit != "MDL"){
		stop("Wrong minimization criterion, must be 'BIC' (default), 'AIC' or 'MDL'")
	}
		
	if(type != "min" && type != "fixed"){
		stop("Wrong changepoint type, must be 'min' (default) or 'fixed'")
	}

M <- M+1
X <- as.matrix(X)
df <- dim(X)[2]
minlength <- df+1
N <- length(Y)
M <- min(M,floor(N/minlength))
Q <- array(Inf, dim=c(M,N,N))
Tau <- array(0, dim=c(M,M))


## Checking Cost function (overall minimization is AIC/BIC --> max logLik or MDL --> min BIC)
if(mincrit=="BIC" || mincrit=="AIC"){
	costfun <- "minell"
	minell <- function(x){-logLik(x)}
	}
if(mincrit=="MDL"){
	costfun <- "BIC"
}

## Except models that overfit the data
## Initialization
for(i in 1:(N-minlength)){
	for(j in (i+minlength):N){
		Mod <- lm(Y[i:j] ~ X[i:j,] - 1)
		Q[1,i,j] <- eval(call(costfun,Mod))
	}
}


## Running algorithm
for(m in 2:M){
	##Step 1
	for(j in 1:N){
		MinVec <- rep(Inf,j-1)
			for(k in 1:(j-1)){
				MinVec[k] <- Q[(m-1),1,k] + Q[1,(k+1),j]
			}
		Q[m,1,j] <- min(MinVec)
	}
	##Step 2
		MinTau <- rep(Inf,N)
			for(k in 1:(N-1)){
				MinTau[k] <- Q[(m-1),1,k] + Q[1,(k+1),N]
			}
		Tau[m,1] <- which.min(MinTau)
	##Step 3
	if(m>2){
	for(i in 2:(m-1)){
		MinTau2 <- rep(Inf,N)
			for(k in 1:(Tau[m,(i-1)]-1)){
				MinTau2[k] <- Q[(m-i),1,k] + Q[1,(k+1),Tau[m,(i-1)]]
			}
		Tau[m,i] <- which.min(MinTau2)
	}
	}
}

## Checking if fixed or best number of changepoints
	if(type=="min"){
		if(mincrit=="MDL"){
			Minis <- Q[,1,N]/2 + log(1:M) + (1:M)*log(N)
			Ncp <- which.min(Minis)
		}
		if(mincrit=="BIC"){
			Minis <- 2*Q[,1,N]+(1:M*minlength+0:(M-1))*log(N)
			Ncp <- which.min(Minis)
		}
		if(mincrit=="AIC"){
			Minis <- 2*Q[,1,N] + 2*(1:M*minlength+0:(M-1))
			Ncp <- which.min(Minis)
		}
	}

	if(type=="fixed"){
		if(mincrit=="MDL"){
			Minis <- Q[,1,N]/2 + log(1:M) + (1:M)*log(N)
			Ncp <- M
		}
		if(mincrit=="BIC"){
			Minis <- 2*Q[,1,N]+(1:M*minlength+0:(M-1))*log(N)
			Ncp <- M
		}
		if(mincrit=="AIC"){
			Minis <- 2*Q[,1,N] + 2*(1:M*minlength+0:(M-1))
			Ncp <- M
		}
	}
	
## Creating models for choosen number of changepoints	
	#print(Minis)
	#print(Tau[Ncp,])
	Models <- list()
	for(l in Ncp:1){
		if(l > 1){
		Models[[l]] <- lm(Y[(Tau[Ncp,l]+1):Tau[Ncp,l-1]] ~ X[(Tau[Ncp,l]+1):Tau[Ncp,l-1],]-1)
		}
		else{Models[[l]] <- lm(Y[(Tau[Ncp,l]+1):N] ~ X[(Tau[Ncp,l]+1):N,]-1)}
	}
	
## Calculating LogLikelihood, AIC/BIC value	
	if(mincrit=="BIC" || mincrit=="AIC"){
		LL <- -Q[Ncp,1,N]
		aic <- -2*LL + (Ncp*minlength+Ncp-1)*2
		bic <- -2*LL + (Ncp*minlength+Ncp-1)*log(N)
	}
	if(mincrit=="MDL"){
		dl <- 0
		for(i in Ncp:1){
			if(i>1){
				dl <- dl + log(Tau[Ncp,i-1] - Tau[Ncp,i])
			}
			else{dl <- dl + log(N - Tau[Ncp,i])}
		}
		LL <- (Q[Ncp,1,N] - minlength*dl)/-2
		aic <- -2*LL + (Ncp*minlength+Ncp-1)*2
		bic <- -2*LL + (Ncp*minlength+Ncp-1)*log(N)
	}	
	
## Calculate predicted values	
	ypred <- rep(0,N)
	for(i in 1:N){ypred[i] <- sum(coef(Models[[which(Tau[Ncp,] < i)[1]]])%*%X[i,])}
	
## Calculate Mean Squared Error
	mse <- mean((Y-ypred)^2)	
	
	if(Qout == F){
		Lijst <- list("Tau"=Tau[Ncp,1:Ncp], val=Minis, "LogLikelihood"=LL, "AIC"=aic, "BIC"=bic, "Models"=Models, "Predicted"=ypred, "MSE"=mse)
		names(Lijst)[2] <- paste(mincrit," values")
		return(Lijst)
	}
	
	if(Qout == T){
		Lijst <- list("Tau"=Tau[Ncp,1:Ncp], val=Minis, "LogLikelihood"=LL, "AIC"=aic, "BIC"=bic, "Models"=Models, "Predicted"=ypred, "MSE"=mse,"Q"=Q)
		names(Lijst)[2] <- paste(mincrit," values")
		return(Lijst)
	}

}

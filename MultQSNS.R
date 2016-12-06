MultQSNS <- function(SNSModsList,mincrit, Nseg = NaN){
	Nnodes <- length(SNSModsList)
	N <- dim(SNSModsList[[1]]$Q)[3]
	M <- Inf
	npar <- 0
	for(i in 1:Nnodes){
		M <- min(M,dim(SNSModsList[[i]]$Q)[1])
		npar <- npar + (1+length(coef(SNSModsList[[i]]$Models[[1]])))
	}
	Q <- array(0,dim=c(M,N,N))
	for(i in 1:Nnodes){
		Q <- Q + SNSModsList[[i]]$Q[1:M,,]
	}	
	if(is.nan(Nseg) != T){Nseg <- min(M,Nseg)}
	
	Tau <- array(0,dim=c(M,M))
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
	
	if(mincrit=="MDL"){
		Minis <- Q[,1,N]/2 + log(1:M) + (1:M)*log(N*Nnodes)
		Ncp <- which.min(Minis)
	}
	if(mincrit=="BIC"){
		Minis <- 2*Q[,1,N]+(1:M*npar+0:(M-1))*log(N*Nnodes)
		Ncp <- which.min(Minis)
	}
	if(mincrit=="AIC"){
		Minis <- 2*Q[,1,N] + 2*(1:M*npar+0:(M-1))
		Ncp <- which.min(Minis)
	}
	if(is.nan(Nseg) != T){
		Ncp <- Nseg
	}
		
	
	if(mincrit=="BIC" || mincrit=="AIC"){
		LL <- -Q[Ncp,1,N]
		aic <- -2*LL + (Ncp*npar+Ncp-1)*2
		bic <- -2*LL + (Ncp*npar+Ncp-1)*log(N*Nnodes)
	}
	if(mincrit=="MDL"){
		dl <- 0
		for(i in Ncp:1){
			if(i>1){
				dl <- dl + log(Tau[Ncp,i-1] - Tau[Ncp,i])
			}
			else{dl <- dl + log(N - Tau[Ncp,i])}
		}
		LL <- (Q[Ncp,1,N] - npar*dl)/-2
		aic <- -2*LL + (Ncp*npar+Ncp-1)*2
		bic <- -2*LL + (Ncp*npar+Ncp-1)*log(N*Nnodes)
	}
	
	return(list("bic"=bic,"Tau"=Tau[Ncp,1:Ncp]))
	
}
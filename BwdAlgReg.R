BwdAlgReg <- function(y,X,A,B,p,Cscaling){
	M <- length(p)
	d <- dim(B)[2]
	N <- sum(!is.na(y))
	X <- as.matrix(X)
	beta <- array(NA,c(N,M))
	beta[N,] <- rep(1,M) * Cscaling[N]
	for(n in (N-1):1){
		for(i in 1:M){
			beta[n,i] <- 0
			for(j in 1:M){
				beta[n,i] <- beta[n,i] + A[i,j]*dnorm(y[n+1],B[j,1:(d-1)]%*%X[n+1,],B[j,d])*beta[n+1,j]
			}
		}
		beta[n,] <- beta[n,]*Cscaling[n]
	}
	return(beta)
}

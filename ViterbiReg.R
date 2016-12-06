ViterbiReg <- function(y,X,A,B,p){
	M <- length(p)
	d <- dim(B)[2]
	N <- sum(!is.na(y))
	X <- as.matrix(X)
	delta <- vector("list",N)
	psi <- vector("list",N)
	shat <- rep(NA,N)
	##Step 1
	for(i in 1:M){
		delta[[1]][i] <- log(p[i]*dnorm(y[1],B[i,1:(d-1)]%*%X[1,],B[i,d]))
	}
	##Step 2
	for(n in 1:(N-1)){
		for(j in 1:M){
			ToMax <- delta[[n]] + log(A[,j])
			delta[[n+1]][j] <- log(dnorm(y[n+1],B[j,1:(d-1)]%*%X[n+1,],B[j,d])) + max(ToMax)
			psi[[n+1]][j] <- which.max(ToMax)
		}
	}
	##Step 3
	LL <- max(delta[[N]])
	shat[N] <- which.max(delta[[N]])
	
	##Step 4
	for(n in (N-1):1){
		shat[n] <- psi[[n+1]][shat[n+1]]
	}
	return(shat)
}

FwdAlgReg <- function(y,X,A,B,p){
	M <- length(p)
	d <- dim(B)[2]
	N <- sum(!is.na(y))
	X <- as.matrix(X)
	alpha <- array(NA,c(N,(M+1)))
	Cscaling <- rep(0,N)
	for(i in 1:M){
		alpha[1,i] <- p[i]*dnorm(y[1],B[i,1:(d-1)]%*%X[1,],B[i,d])
	}
  	Cscaling[1] <- 1/sum(alpha[1,1:M])
  
	alpha[1,1:M] <- Cscaling[1] * alpha[1,1:M]
	for(n in 1:(N-1)){
		for(j in 1:M){
			alpha[n+1,j] <- 0
			for(i in 1:M){
			  	alpha[n+1,j] <- alpha[n+1,j] + alpha[n,i]*A[i,j]
  
			}
			alpha[n+1,j] <- alpha[n+1,j]*dnorm(y[n+1],B[j,1:(d-1)]%*%X[n+1,],B[j,d])
		}
	    Cscaling[n+1] <- 1/sum(alpha[n+1,1:M])
      alpha[n+1,] <- Cscaling[n+1]*alpha[n+1,]

	}
	alpha[,M+1] <- Cscaling
	return(alpha)
}

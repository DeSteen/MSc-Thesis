EMAlg <- function(Y,M){
	mu <- sample(Y,M)
	sigma <- runif(M,0,0.5)
	theta <- cbind(mu,sigma)
	pihat <- runif(1,0,1)
	
	gamma <- c()
	LLnew <- c()
	LLnew[1] <- LL0 <- LLold <- sum(log((1-pihat)*dnorm(Y,mu[1],sigma[1])+pihat*dnorm(Y,mu[2],sigma[2])))
	err <- 1
	k <- 2
	while(err > 1e-8){
		for(i in 1:length(Y)){
			gamma[i] <- (pihat*dnorm(Y[i],mu[2],sigma[2]))/((1-pihat)*dnorm(Y[i],mu[1],sigma[1])+pihat*dnorm(Y[i],mu[2],sigma[2]))
		}
	
	
	mu[1] <- sum(Y - gamma*Y)/sum(1-gamma)
	mu[2] <- sum(gamma*Y)/sum(gamma)
	sigma[1] <- sqrt(sum((1-gamma)*(Y-mu[1])^2)/sum(1-gamma))
	sigma[2] <- sqrt(sum(gamma*(Y-mu[2])^2)/sum(gamma))
	pihat <- sum(gamma)/length(Y)
	LLnew[k] <- sum(log((1-pihat)*dnorm(Y,mu[1],sigma[1])+pihat*dnorm(Y,mu[2],sigma[2])))
	
	err <- abs((LLnew[k] - LLold)/(LLnew[k] - LL0))
	LLold <- LLnew[k]
	k <- k+1
		}
	
	return(list(mu,sigma,pihat,LLnew))	
	
}
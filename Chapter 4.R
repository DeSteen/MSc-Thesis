library(MASS)
library(xtable)
source("BaumWelchReg.R")
source("ViterbiReg.R")
source

##Stock Price example
set.seed(1)
Time <- seq(0,40,0.125)
y <- c()
y[1] <- rnorm(1,20,5)
beta <- matrix(c(20,0,-4,3,124,-5,-44,2,20,0),ncol=2,byrow=T)
for(i in 2:321){y[i] <- sum(beta[ceiling((i-1)/64),]%*%c(1,Time[i])) + rnorm(1,0,3)}

## Figure 1
symbols(c(2,4,6,8,2,4,6,8),c(1,1,1,1,2,2,2,2),circles=rep(0.5,8),inches=F,ylim=c(0,3),axes=F,xlab=NA,ylab=NA)
text(c(2,4,6,8,2,4,6,8),c(1,1,1,1,2,2,2,2),c(expression("Y"[1]),expression("Y"[2]),expression("Y"[3]),expression("Y"[4]),expression("S"[1]),expression("S"[2]),expression("S"[3]),expression("S"[4])))
arrows(x0=c(2.5,4.5,6.5,8.5),y0=rep(2,4),x1=c(3.5,5.5,7.5,9.5),y1=rep(2,4),length=0.15)
arrows(x0=c(2,4,6,8),y0=rep(1.77,4),x1=c(2,4,6,8),y1=rep(1.23,4),length=0.15)


### HMM Stock Price example 
# M1 <- BaumWelchReg(y,cbind(1,Time),4,-100000,1)
# for(i in 1:100){
	# newM <- BaumWelchReg(y,cbind(1,Time),4,-100000,1)
	# if(newM$BIC < M1$BIC){M1 <- newM}
	# else{M1 <- M1}
	# if(i %% 10 == 0){print(i)}
# }
# save(M1,file="Example_HMM4.Rda")
load("Example_HMM4.Rda")

M1seq <- ViterbiReg(y,cbind(1,Time),M1$A,M1$B,M1$pi)
ypredhmm1 <- c()
for(i in 1:length(M1seq)){
	ypredhmm1[i] <- sum(M1$B[M1seq[i],1:2]%*%c(1,Time[i])) 
}

M1
print(xtable(t(as.matrix(M1$pi)),align=rep("c",5)),include.rownames=F,include.colnames=F)
print(xtable(M1$A,align=rep("c",5),digits=3),include.rownames=F,include.colnames=F)
print(xtable(M1$B,align=rep("c",4),digits=3),include.rownames=F,include.colnames=F)
print(xtable(t(as.matrix(M1seq)),align=rep("c",282)),include.rownames=F,include.colnames=F)

#Figure 2
plot(Time,y,type="l",ylab="Stock Price",ylim=c(1,46))
title("Example Stock Prices")
lines(Time,ypredhmm1,type="l",col="red",lwd=2)
legend("topright",c("Actual Stock Price","HMM Regression line"),col=c("black","red"),lty=rep(1,2))


### HMM with equal variances Stock Price example
# M2 <- BaumWelchReg(y,cbind(1,Time),4,-100000,1,eq.var=T)
# for(i in 1:50){
	# newM <- BaumWelchReg(y,cbind(1,Time),4,-100000,1,eq.var=T)
	# if(newM$BIC < M2$BIC){M2 <- newM}
	# else{M2 <- M2}
	# if(i %% 10 == 0){print(i)}
# }
# save(M2,file="Example_HMM4_EqVar.Rda")
load("Example_HMM4_EqVar.Rda")

M2seq <- ViterbiReg(y,cbind(1,Time),M2$A,M2$B,M2$pi)
ypredhmm2 <- c()
for(i in 1:length(M2seq)){
	ypredhmm2[i] <- sum(M2$B[M2seq[i],1:2]%*%c(1,Time[i])) 
}

M2
print(xtable(t(as.matrix(M2$pi)),align=rep("c",5)),include.rownames=F,include.colnames=F)
print(xtable(M2$A,align=rep("c",5),digits=3),include.rownames=F,include.colnames=F)
print(xtable(M2$B,align=rep("c",4),digits=3),include.rownames=F,include.colnames=F)
print(xtable(t(as.matrix(M2seq)),align=rep("c",282)),include.rownames=F,include.colnames=F)

plot(Time,y,type="l",ylab="Stock Price",ylim=c(1,46))
title("Example Stock Prices")
lines(Time,ypredhmm2,type="l",col="red",lwd=2)
legend("topright",c("Actual Stock Price","HMM-EqVar Regression line"),col=c("black","red"),lty=rep(1,2))


#Figure 4 and 5
M1stuck <- BaumWelchReg(y,cbind(1,Time),4,-100000,1)
M2stuck <- BaumWelchReg(y,cbind(1,Time),4,-100000,1,eq.var=T)
M1stuckseq <- ViterbiReg(y,cbind(1,Time),M1stuck$A,M1stuck$B,M1stuck$pi)
M2stuckseq <- ViterbiReg(y,cbind(1,Time),M2stuck$A,M2stuck$B,M2stuck$pi)

ypredhmm1stuck <- ypredhmm2stuck <- c()
for(i in 1:length(M1stuckseq)){
	ypredhmm1stuck[i] <- sum(M1stuck$B[M1stuckseq[i],1:2]%*%c(1,Time[i])) 
	ypredhmm2stuck[i] <- sum(M2stuck$B[M2stuckseq[i],1:2]%*%c(1,Time[i])) 
}

plot(Time,y,type="l",ylab="Stock Price",ylim=c(1,46))
title("Example Stock Prices")
lines(Time,ypredhmm1stuck,type="l",col="red",lwd=2)
legend("topright",c("Actual Stock Price",paste("HMM Regression line")),col=c("black","red"),lty=rep(1,2))

plot(Time,y,type="l",ylab="Stock Price",ylim=c(1,46))
title("Example Stock Prices")
lines(Time,ypredhmm2stuck,type="l",col="red",lwd=2)
legend("topright",c("Actual Stock Price",paste("HMM-EqVar Regression line")),col=c("black","red"),lty=rep(1,2))



### Left-Right HMM Stock Price example
# M3 <- BaumWelchReg(y,cbind(1,Time),5,-100000,1,LR=T)
# for(i in 1:100){
	# newM <- BaumWelchReg(y,cbind(1,Time),5,-100000,1,LR=T)
	# if(newM$BIC < M3$BIC){ M3 <- newM}
	# else{M3 <- M3}
	# if(i %% 10 == 0){print(i)}
# }
# save(M3,file="Example_HMM5_LR.Rda")

load("Example_HMM5_LR.Rda")

M3seq <- ViterbiReg(y,cbind(1,Time),M3$A,M3$B,M3$pi)
ypredhmm3 <- c()
for(i in 1:length(M3seq)){
	ypredhmm3[i] <- sum(M3$B[M3seq[i],1:2]%*%c(1,Time[i])) 
}

#Figure 6
plot(Time,y,type="l",ylab="Stock Price",ylim=c(1,46))
title("Example Stock Prices")
lines(Time,StockPriceSNS4$Predicted,type="l",col="red",lwd=2)
lines(Time,ypredhmm3,type="l",col="green",lwd=2,lty=2)
legend("topright",c("Actual Stock Price","CPM Regression line","LR-HMM Regression line"),col=c("black","red","green"),lty=c(1,1,1))



### Left-Right HMM with equal variance Stock Price example
# M4 <- BaumWelchReg(y,cbind(1,Time),5,-100000,1,LR=T,eq.var=T)
# for(i in 1:100){
	# newM <- BaumWelchReg(y,cbind(1,Time),5,-100000,1,LR=T,eq.var=T)
	# if(newM$BIC < M4$BIC){ M4 <- newM}
	# else{M4 <- M4}
	# if(i %% 10 == 0){print(i)}
# }
# save(M4,file="Example_HMM5_LR_EqVar.Rda")

load("Example_HMM5_LR_EqVar.Rda")

M4seq <- ViterbiReg(y,cbind(1,Time),M4$A,M4$B,M4$pi)
ypredhmm4 <- c()
for(i in 1:length(M4seq)){
	ypredhmm4[i] <- sum(M4$B[M4seq[i],1:2]%*%c(1,Time[i])) 
}

#Figure 7
plot(Time,y,type="l",ylab="Stock Price",ylim=c(1,46))
title("Example Stock Prices")
lines(Time,StockPriceSNS4$Predicted,type="l",col="red",lwd=2)
lines(Time,ypredhmm4,type="l",col="green",lwd=2,lty=2)
legend("topright",c("Actual Stock Price","CPM Regression line","LR-HMM-EqVar Regression line"),col=c("black","red","green"),lty=c(1,1,1))

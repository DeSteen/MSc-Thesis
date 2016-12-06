library(MASS)
library(xtable)
source("BaumWelchReg.R")
source("ViterbiReg.R")

##Example HMM
A <- matrix(c(0.67,0.59,0.33,0.41),ncol=2)
Pi <- c(0.6,0.4)
B <- matrix(c(0.17,0.55,0.83,0.45),ncol=2)
y <- c(2,6,6,4,1); y <- (y!=6)+1

(alpha1Green <- Pi[1]*B[1,y[1]]); (alpha1Red <- Pi[2]*B[2,y[1]])
(alpha2Green <- B[1,y[2]]*sum(c(alpha1Green,alpha1Red)*A[,1])); (alpha2Red <- B[2,y[2]]*sum(c(alpha1Green,alpha1Red)*A[,2]))
(alpha3Green <- B[1,y[3]]*sum(c(alpha2Green,alpha2Red)*A[,1])); (alpha3Red <- B[2,y[3]]*sum(c(alpha2Green,alpha2Red)*A[,2]))
(alpha4Green <- B[1,y[4]]*sum(c(alpha3Green,alpha3Red)*A[,1])); (alpha4Red <- B[2,y[4]]*sum(c(alpha3Green,alpha3Red)*A[,2]))
(alpha5Green <- B[1,y[5]]*sum(c(alpha4Green,alpha4Red)*A[,1])); (alpha5Red <- B[2,y[5]]*sum(c(alpha4Green,alpha4Red)*A[,2]))

(Lalpha <- alpha5Green + alpha5Red)


beta5Green <- beta5Red <- 1
(beta4Green <- sum(c(beta5Green,beta5Red)*A[1,]*B[,y[5]])); (beta4Red <- sum(c(beta5Green,beta5Red)*A[2,]*B[,y[5]]))
(beta3Green <- sum(c(beta4Green,beta4Red)*A[1,]*B[,y[4]])); (beta3Red <- sum(c(beta4Green,beta4Red)*A[2,]*B[,y[4]]))
(beta2Green <- sum(c(beta3Green,beta3Red)*A[1,]*B[,y[3]])); (beta2Red <- sum(c(beta3Green,beta3Red)*A[2,]*B[,y[3]]))
(beta1Green <- sum(c(beta2Green,beta2Red)*A[1,]*B[,y[2]])); (beta1Red <- sum(c(beta2Green,beta2Red)*A[2,]*B[,y[2]]))

(Lbeta <- sum(c(beta1Green,beta1Red)*Pi*B[,y[1]]))


delta1Green <- alpha1Green; delta1Red <- alpha1Red
(delta2Green <- B[1,y[2]]*max(c(delta1Green,delta1Red)*A[,1])); (delta2Red <- B[2,y[2]]*max(c(delta1Green,delta1Red)*A[,2]))
psi2Green <- which.max(c(delta1Green,delta1Red)*A[,1]); psi2Red <- which.max(c(delta1Green,delta1Red)*A[,2])
(delta3Green <- B[1,y[3]]*max(c(delta2Green,delta2Red)*A[,1])); (delta3Red <- B[2,y[3]]*max(c(delta2Green,delta2Red)*A[,2]))
psi3Green <- which.max(c(delta2Green,delta2Red)*A[,1]); psi3Red <- which.max(c(delta2Green,delta2Red)*A[,2])
(delta4Green <- B[1,y[4]]*max(c(delta3Green,delta3Red)*A[,1])); (delta4Red <- B[2,y[4]]*max(c(delta3Green,delta3Red)*A[,2]))
psi4Green <- which.max(c(delta3Green,delta3Red)*A[,1]); psi4Red <- which.max(c(delta3Green,delta3Red)*A[,2])
(delta5Green <- B[1,y[5]]*max(c(delta4Green,delta4Red)*A[,1])); (delta5Red <- B[2,y[5]]*max(c(delta4Green,delta4Red)*A[,2]))
psi5Green <- which.max(c(delta4Green,delta4Red)*A[,1]); psi5Red <- which.max(c(delta4Green,delta4Red)*A[,2])

(shat5 <- which.max(c(delta5Green,delta5Red)))


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

BetaTableLR <- matrix(rep(0,36),ncol=4)
for(i in 1:4){
  BetaTableLR[1:2,i] <- coef(StockPriceSNS$Models[[5-i]])[,1]
  BetaTableLR[3,i] <- StockPriceSNS$Models[[5-i]]$sigma
  BetaTableLR[4:6,i] <- M3$B[i,1:3]
  BetaTableLR[7:9,i] <- M4$B[i,1:3]
}
colnames(BetaTableLR) <- c("Group 1","Group 2","Group 3","Group 4")
xtable(BetaTableLR)

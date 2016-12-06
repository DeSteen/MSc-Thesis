library(MASS)
source("SNS.R")

##Stock Price example
set.seed(1)
Time <- seq(0,40,0.125)
y <- c()
y[1] <- rnorm(1,20,5)
beta <- matrix(c(20,0,-4,3,124,-5,-44,2,20,0),ncol=2,byrow=T)
for(i in 2:321){y[i] <- sum(beta[ceiling((i-1)/64),]%*%c(1,Time[i])) + rnorm(1,0,3)}


# StockPriceSNS4 <- SNS(y,cbind(1,Time),4,type="fixed")
# StockPriceSNS8 <- SNS(y,cbind(1,Time),8,type="fixed")
# StockPriceSNSminAIC <- SNS(y,cbind(1,Time),20,"AIC")
# StockPriceSNSminBIC <- SNS(y,cbind(1,Time),20,"BIC")
# StockPriceSNSminMDL <- SNS(y,cbind(1,Time),20,"MDL")
# save(StockPriceSNS4,StockPriceSNS8,StockPriceSNSminAIC,StockPriceSNSminBIC,StockPriceSNSminMDL,file="Example_SNS.Rda")

load("Example_SNS.Rda")

#Figure 1
plot(Time,y,type="l",ylab="Stock Price",ylim=c(1,46))
title("Example Stock Prices")

#Figure 2
plot(Time,y,type="l",ylab="Stock Price",ylim=c(1,46))
title("Example Stock Prices")
lines(Time,StockPriceSNS4$Predicted,type="l",col="red",lwd=2)
legend("topright",c("Actual Stockprice","CPM Regression line"),col=c("black","red"),lty=c(1,1))
points(StockPriceSNS4$Tau[1:4]*0.125,StockPriceSNS4$Predicted[StockPriceSNS4$Tau],col="red",pch=16)
StockPriceSNS4$Tau[1:4]*0.125

#Figure 3
plot(0:(length(StockPriceSNSminAIC$'AIC  values')-1),StockPriceSNSminAIC$'AIC  values',type="l",xlab="Number of Changepoints",ylab="Minimization Criterion Value",ylim=c(750,1800),col="red",xaxt="n")
segments(x0=c(4,6),y0=c(650,650),y1=c(min(StockPriceSNSminBIC$'BIC  values'),min(StockPriceSNSminAIC$'AIC  values')),lty=2)
axis(1,at=c(0,4,6,10,15,20))
points(which.min(StockPriceSNSminAIC$'AIC  values')-1,min(StockPriceSNSminAIC$'AIC  values'),pch=16,col="red")
lines(0:(length(StockPriceSNSminBIC$'BIC  values')-1),StockPriceSNSminBIC$'BIC  values',col="green")
points(which.min(StockPriceSNSminBIC$'BIC  values')-1,min(StockPriceSNSminBIC$'BIC  values'),pch=16,col="green")
lines(0:(length(StockPriceSNSminMDL$'MDL  values')-1),StockPriceSNSminMDL$'MDL  values',col="blue")
points(which.min(StockPriceSNSminMDL$'MDL  values')-1,min(StockPriceSNSminMDL$'MDL  values'),pch=16,col="blue")
legend(17,1750,c("AIC","BIC", "MDL"),col=c("red","green","blue"),lty=c(1,1,1))

#Figure 4
plot(Time,y,type="l",ylab="Stock Price",ylim=c(1,46))
title("Example Stock Prices")
lines(Time,StockPriceSNSminAIC$Predicted,lwd=2,col="red")
lines(Time,StockPriceSNSminBIC$Predicted,lwd=2,col="green")
lines(Time,StockPriceSNSminMDL$Predicted,lwd=2,col="blue",lty=2)
legend("topright",c("AIC Regression line","BIC Regression line", "MDL Regression line"),col=c("red","green","blue"),lty=c(1,1,1))


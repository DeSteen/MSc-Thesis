library(MASS)
library(igraph)
source("EMAlg.R")

########## Section Dynamic Programming ###########
#Example city / figure 1
edges <- matrix(rep(0,7*7),ncol=7)
edges[upper.tri(edges)] <- c(1,1,0,1,0,0,0,1,1,1,0,1,1,1,0,0,0,0,0,1,1)
g <- graph.adjacency(edges,mode="undirected")
sq.lay <- matrix(c(0,0,0,1,0.5,0.5,1,0,1,2,2,1,2,2),ncol=2,byrow=T)
plot(g,layout=sq.lay,vertex.label=LETTERS[1:7],edge.label=c(5,2,4,7,8,8,12,9,8,5,2),vertex.color="cyan")

## Step 1
T_E <- 5
T_F <- 2

## Step 2
T_B <- min(c(7+T_E,8+T_F))
T_C <- min(c(8+T_E,12+T_F))
T_D <- min(c(9+T_E,8+T_F))

## Step 3
T_A <- min(c(5+T_B,2+T_C,4+T_D))

#Figure 2
plot(g,layout=sq.lay,vertex.label=LETTERS[1:7],edge.label=c(5,2,4,7,8,8,12,9,8,5,2),edge.color=c("red","red","green","red","green","green","red","red","green","green","green"),vertex.color="cyan")

########## Section EM-algorithm ##########
set.seed(4)
x <- sample(c(4,7),50,rep=T); p <- sum(x==7)/50
x <- x + rnorm(50,0,0.7)

#Figure 3
hist(x,prob=T,breaks=10,xlim=c(1,11),main="Histogram of example data")

#Figure 4
t <- seq(0,12,0.05)
d <- (1-p)*dnorm(t,4,0.7) + p*dnorm(t,7,0.7)
plot(t,d,type="l",xlim=c(1,11),ylim=c(0,0.4),ylab=expression(paste(f['Y'],'(t;',theta,')')))

# (EM_EST <- EMAlg(x,2))
# (EM_EST2 <- EMAlg(x,2))
# (EM_EST3 <- EMAlg(x,2))
# save(EM_EST, file="EM_EST.Rda")
# save(EM_EST2,file="EM_EST2.Rda")
# save(EM_EST3,file="EM_EST3.Rda")

load("EM_EST.Rda")
load("EM_EST2.Rda")
load("EM_EST3.Rda")

d_est <- (1-EM_EST[[3]])*dnorm(t,EM_EST[[1]][1],EM_EST[[2]][1]) + EM_EST[[3]]*dnorm(t,EM_EST[[1]][2],EM_EST[[2]][2])
points(t,d_est,type="l",col="red")

#Figure 5
plot(1:length(EM_EST[[4]]),EM_EST[[4]], ylab="Log-Likelihood", xlab="Number of iterations", type="l", ylim=c(-100,-75),xlim=c(0,65),col="blue")
lines(1:length(EM_EST2[[4]]),EM_EST2[[4]], col="orange")
lines(1:length(EM_EST3[[4]]),EM_EST3[[4]], col="red")
legend("topright",c("Converged in 23 iterations to -77.99","Converged in 64 iterations to -83.78","Converged in 16 iterations to -86.64"),col=c("blue","orange","red"),lty=c(1,1,1))
points(1:length(EM_EST[[4]]),EM_EST[[4]],col="blue",pch=19)
points(1:length(EM_EST2[[4]]),EM_EST2[[4]], col="orange",pch=19)
points(1:length(EM_EST3[[4]]),EM_EST3[[4]], col="red", pch=19)

#Figure 6
plot(t,d,type="l",xlim=c(1,11),ylim=c(0,0.55),ylab=expression(paste(f['Y'],'(t;',theta,')')))
d_est <- (1-EM_EST[[3]])*dnorm(t,EM_EST[[1]][1],EM_EST[[2]][1]) + EM_EST[[3]]*dnorm(t,EM_EST[[1]][2],EM_EST[[2]][2])
points(t,d_est,type="l",col="blue")
d_est2 <- (1-EM_EST2[[3]])*dnorm(t,EM_EST2[[1]][1],EM_EST2[[2]][1]) + EM_EST2[[3]]*dnorm(t,EM_EST2[[1]][2],EM_EST2[[2]][2])
points(t,d_est2,type="l",col="orange")
d_est3 <- (1-EM_EST3[[3]])*dnorm(t,EM_EST3[[1]][1],EM_EST3[[2]][1]) + EM_EST3[[3]]*dnorm(t,EM_EST3[[1]][2],EM_EST3[[2]][2])
points(t,d_est3,type="l",col="red")
legend("topright",c("True density","LL: -77.99", "LL: -83.78", "LL: -86.64"),col=c("black","blue","orange","red"),lty=c(1,1,1,1))

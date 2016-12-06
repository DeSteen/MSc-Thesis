RemoveOutlier <- function(x,k=1.5){
	q <- quantile(x,prob=c(0.25,0.75))
	lb <- q[1] - k*(q[2]-q[1])
	ub <- q[2] + k*(q[2]-q[1])
	y <- x
	y[((x > lb) + (x < ub))!=2] <- NA
	return(y)
}
gplotyeast <- function(adjmMat){
	CL  <- layout.circle(gr <- graph.adjacency(adjmMat))
plot(gr,layout=CL,vertex.label=c("CBF1","GAL4","SWI5","GAL80","ASH1"),vertex.size=25,edge.color="black",vertex.color="cyan",edge.curved=0.2)
}

gplotad <- function(adjmMat){
	adl <- matrix(c(1,0,2,0.5,1.5,0,2,-0.5,2.5,0.5),ncol=2,byrow=T)
plot(graph.adjacency(adjmMat),layout=adl,vertex.label=LETTERS[1:5],vertex.size=25,edge.color="black",vertex.color="cyan",edge.curved=0)
}
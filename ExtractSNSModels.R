ExtractSNSModels <- function(SNSModelList){
	out <- list()
	for(i in 1:length(SNSModelList)){
		out[[i]] <- SNSModelList[[i]][c(1,6)]
	}
	return(out)
}
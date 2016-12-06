auc.pr <- function(score.matrix,relation.matrix,type){

EMscore_v = as.vector(score.matrix) 
relation_v = as.vector(relation.matrix)
relation_v = relation_v[!is.na(EMscore_v)]
EMscore_v = EMscore_v[!is.na(EMscore_v)]
order_score = order(EMscore_v)
relation_v = relation_v[order(EMscore_v, partial = relation_v, decreasing = TRUE)]
EMscore_v = sort(EMscore_v, decreasing = TRUE)

rcc <- pr.curve(EMscore_v[relation_v==T],EMscore_v[relation_v==F],curve=T)
print(rcc)
plot(rcc)

}
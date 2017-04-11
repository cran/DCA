normscore.row <-
function(a)
{
    for(i in 1:nrow(a)) a[i,]<-normal_trafo(a[i,])
    return(a)
}

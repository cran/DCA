normscore.row.small <-
function(a)
{
    b<-t(apply(a, 1, normal_trafo))
    return(b)
}

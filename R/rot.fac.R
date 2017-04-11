rot.fac <-
function(x, B)	# proj is projected data, B is factor scores
{
    l<-sqrt(apply(x^2, 1,sum))
    x<-x/l
    y<-GPFoblq(x, method="oblimin")
    return(y$Th)
}

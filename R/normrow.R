normrow <-
function(a, mean.rm=T)
{
    if(mean.rm)
    {
        m<-apply(a,1,mean,na.rm=T)
        s<-apply(a,1,sd,na.rm=T)
        a<-(a-m)/s #/sqrt(ncol(a)-1)
    }else{
        a<-a-mean(a)
        ss<-apply(a^2,1,sum)
        a<-a/sqrt(ss)
    }
    return(a)
}

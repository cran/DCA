la.simu.gen <-
function(n, p, n.grp,n.noise.gene, rho, pwr)
{
    la.simu.one.group<-function(n, p, rho=0.5, pwr=0.5, sub.grps=10)
    {

        p<-p/sub.grps
        
        #z<-rnorm(n)
        #z<-sort(z)
        #z<-(z-min(z))/((max(z)-min(z)))
        #z<-(z-0.5)*1.98
        
        z<-runif(n,min=-1, max=1)
        
        for(m in 1:sub.grps)
        {
            
            dat<-matrix(0, ncol=n, nrow=p)
            for(j in 1:n)
            {
                gamma<-sign(z[j])*abs(z[j])^pwr
                xy<-rmvnorm(1, mean=c(0,0), sigma=matrix(c(1, gamma, gamma,1),ncol=2))
                dat[,j]<- c(rep(xy[1],p/2), rep(xy[2],p/2)) + rnorm(p, mean=0, sd=rho)
            }
            
            if(m == 1)
            {
                all.dat<-dat
                
            }else{
                all.dat<-rbind(all.dat,dat)
            }
            
        }
        
        dat<-all.dat
        o<-sample(n, n, replace=F)
        
        z<-z[o]
        dat<-dat[,o]
        
        z<-qnorm((z+1)/2)
        return(list(dat=dat, z=z))
    }
    
    r<-la.simu.one.group(n,p,rho,pwr)
    dat<-r$dat
    z<-r$z
    
    for(i in 2:n.grp)
    {
        r<-la.simu.one.group(n,p,rho,pwr)
        dat<-rbind(dat, r$dat)
        z<-cbind(z, r$z)
    }
    
    dat<-rbind(dat, matrix(rnorm(n.noise.gene*n), ncol=n))
    return(list(dat=dat, z=z))
    
}

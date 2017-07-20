
find.xy<-function(array, z, fdr.cut=0.05, normalization="standardize", center.z=FALSE, lac.percentile=0.8)
# z is the output from find.z(), each column is a factor
{
    
    if(normalization == "normal score") array<-normscore.row(array)
    if(normalization == "standardize") array<-normrow(array)
    
    la.rec<-new("list")
    
    ccc.abs<-cor(abs(t(array)))
    abs.ccc<-abs(cor(t(array)))
    ccc.diff<-ccc.abs-abs.ccc
    rm(ccc.abs)
    rm(abs.ccc)
    gc()
    
    
    allow.mat<-which(ccc.diff>=quantile(as.dist(ccc.diff), lac.percentile), arr.ind=TRUE)
    rm(ccc.diff)
    gc()
    
    allow.mat<-allow.mat[allow.mat[,1]>allow.mat[,2],]
    allow.mat.2<-allow.mat
    if(nrow(allow.mat)>50000) allow.mat.2<-allow.mat[sample(nrow(allow.mat), 50000, replace=FALSE),]
    
    #    ccc<-cor(t(array))
    #    ccc.range<-quantile(as.dist(ccc),c(0.5-(1-lac.percentile)/2, 0.5+(1-lac.percentile)/2))
    #    background.mat<-which(ccc > ccc.range[1] & ccc < ccc.range[2], arr.ind=TRUE)
    #    background.mat<-background.mat[background.mat[,1]>background.mat[,2],]
    #    if(nrow(background.mat)>50000) background.mat<-background.mat[sample(nrow(background.mat), 50000, replace=FALSE),]
    
    
    for(n in 1:ncol(z))
    {
        message("############### working on factor ", n, " ##################")
        if(center.z) z[,n]<-z[,n]-mean(z[,n])
        b<-one.gene.LA(array,z[,n])
        diag(b)<-0
        
        s<-b[allow.mat.2]
        #b.background<-b[background.mat]
        
        #s<-sample(which(allow.mat==1), 50000, replace=FALSE)
        #s<-b[s]
        l<-locfdr(s, nulltype=3)
        
        s2<-s[l$fdr<=fdr.cut]
        cut.low<-max(s2[s2<0])
        cut.high<-min(s2[s2>0])
        
        d<-cbind(allow.mat, b[allow.mat])
        
        d2<-approx(x=s, y=l$fdr, xout=d[,3], rule=2)
        d[,3]<-d2$y
        d<-d[which(d[,3]<=fdr.cut),]
        
        la.rec[[n]]<-d
        message(" ")
        
        gc()
    }
    la.rec
}


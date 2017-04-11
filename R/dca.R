dca <-
function(array, top.pairs.prop=0.95, max.pairs=1e6, n.fac=10, sumabsv=sqrt(max.pairs)/10, normalization="standardize", method="PCA")
{

    if(normalization == "normal score") array<-normscore.row(array)
    if(normalization == "standardize") array<-normrow(array)
    
    n.pairs<-min(choose(nrow(array),2)*(1-top.pairs.prop), max.pairs)
    
    # make a new matrix, each row is the X*Y from two random genes
    
    ccc.abs<-cor(abs(t(array)))
    abs.ccc<-abs(cor(t(array)))
    ccc.diff<-ccc.abs-abs.ccc
    rm(ccc.abs)
    rm(abs.ccc)
    gc()
    pos<-which(ccc.diff>quantile(ccc.diff, 1-2*n.pairs/nrow(array)/nrow(array)), arr.ind=T)
    pos<-pos[which(pos[,1] > pos[,2]),]
    #if(nrow(pos)>n.pairs) pos<-pos[sample(nrow(pos), n.pairs, replace=F),]
    rm(ccc.diff)
    gc()
    
    b<-array[pos[,1],] * array[pos[,2],]
    for(i in 1:nrow(b)) b[i,]<-normal_trafo(b[i,])
    
    
    # find the dominant factors in the new matrix
    # b<-normscore.row(b)
    
    if(method=="PCA")
    {
        ccc<-t(b) %*% b
        e<-eigen(ccc)
        fac<-e$vec[,1:n.fac]
    }else if(method=="SPCA"){
        #cv.out <- SPC.cv(t(b), sumabsvs=seq(1, sqrt(nrow(b))/10, len=12)[c(-1,-12)])
        #aa<-SPC(t(b),sumabsv=cv.out$bestsumabsv,K=n.fac)
        aa<-SPC(t(b),sumabsv=sumabsv,K=n.fac)
        fac<-aa$u
    }else if(method=="kmeans"){
        fac<-t(kmeans(b, centers=n.fac, iter.max=50, nstart=1)$centers)
    }
    
    if(method != "kmeans")
    {
        if(nrow(b)>1e5) b<-b[sample(nrow(b),1e5, replace=FALSE),]
        
        proj<-b %*% fac
        
        rot<-diag(ncol(fac))
        try(rot<-varimax(proj)$rot)
        fac2<-fac %*% rot
    }else{
        fac2<-fac
    }
    
    proj<-b %*% fac2
    proj.len<-apply(proj^2,2,sum)
    o<-order(-proj.len)
    fac2<-fac2[,o]
    ss.proj<-proj.len[o]/sum(proj.len)
    
    return(list(fac=fac, rotated=fac2, ss.proj=ss.proj))
}

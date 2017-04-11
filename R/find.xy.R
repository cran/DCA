find.xy <-
function(array, z, fdr.cut=0.05, normalization="standardize", center.z=FALSE, low.cor.percentile=c(0.2, 0.8))
# z is the output from find.z(), each column is a factor
{

    
    if(normalization == "normal score") array<-normscore.row(array)
    if(normalization == "standardize") array<-normrow(array)
    
    la.rec<-new("list")
    
    cor.mat<-cor(t(array))
    diag(cor.mat)<-NA
    cor.cuts<-quantile(cor.mat, low.cor.percentile,na.rm=T)
    allow.mat<-1*(cor.mat>=cor.cuts[1] & cor.mat<=cor.cuts[2])
    diag(allow.mat)<-0
    
    for(n in 1:ncol(z))
    {
        message("############### working on factor ", n, " ##################")
        if(center.z) z[,n]<-z[,n]-mean(z[,n])
        b<-one.gene.LA(array,z[,n])
        diag(b)<-0
        
		s<-which(allow.mat==1)
        if(length(s)>50000) s<-sample(s, 50000, replace=FALSE)
        s<-b[s]
        l<-locfdr(s, nulltype=3)
        
        s2<-s[l$fdr<=fdr.cut]
        cut.low<-max(s2[s2<0])
        cut.high<-min(s2[s2>0])
        
        sel<-which(b>=cut.high | b<=cut.low, arr.ind =T)
        sel<-sel[sel[,1]>sel[,2],]
        d<-cbind(sel, b[sel])
        d2<-approx(x=s, y=l$fdr, xout=d[,3], rule=2)
        d[,3]<-d2$y
        
        la.rec[[n]]<-d
        message(" ")
        
        gc()
    }
    la.rec
}

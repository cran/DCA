plot_la <-
function(x,y,z, use.locfdr=FALSE, cols=c("red","green","blue"),cex=0.5)
{
    
	addTrans<-function(color,trans)
	{
    
    if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
    if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
    if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
    
    num2hex <- function(x)
    {
        hex <- unlist(strsplit("0123456789ABCDEF",split=""))
        return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
    }
    rgb <- rbind(col2rgb(color),trans)
    res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
    return(res)
	}

    if(use.locfdr)
    {
        l<-locfdr(z, plot=0)
        grp<-rep(2, length(z))
        grp[z<=max(z[z<median(z) & l$fdr<=0.5])]<-1
        grp[z>=min(z[z>median(z) & l$fdr<=0.5])]<-3
    }else{
        cuts<-quantile(z, c(0.3333, 0.6666))
        grp<-rep(2, length(z))
        grp[z<cuts[1]]<-1
        grp[z>cuts[2]]<-3
    }
    
    plot(x,y, type="n")
    
    all.grp<-unique(grp)
    all.grp<-all.grp[order(all.grp)]
    
    r<-NULL

#   message(c("LA score: ", sum(x*y*z)))
    r<-c(r,paste("LA score: ", signif(sum(x*y*z),3)))
    
    for(i in 1:length(all.grp))
    {
        s<-which(grp == all.grp[i])
        points(x[s],y[s],col=addTrans(cols[i],120), pch=19, cex=cex)
        
        r<-c(r, paste("group:", i,"(", signif(median(z[s]),2), ")", cols[i], ", corr=", round(cor(x[s],y[s]), 3)))
    }
    r
}

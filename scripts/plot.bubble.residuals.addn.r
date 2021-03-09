plot.bubble.fun<-function(datao,datap,nlen=22,sampsize,maxres=5){
  ny<-length(as.numeric(unlist(datao[,1])))
#this is standardized pearson residual (o-p)/(sqrt(p(1-p)/sample size))
  sampsize[sampsize<1.0]=1.0
  res<-(datao[,2:(nlen+1)]-datap[,2:(nlen+1)])/sqrt((datap[,2:(nlen+1)]*(1-datap[,2:(nlen+1)]))/sampsize)
  res[is.na(res)]<-0

resp<-res
rangeres=range(resp)
resp[resp<0]=0
respstar=resp
respstar[respstar<=maxres]=0
resp[resp>maxres]=0
resn<-res
resn[resn>=0]=0
resnstar=resn
resnstar[resnstar<=-maxres]=0
resn[resn<(-maxres)]= 0  

  plot(rep(as.numeric(unlist(datao[,1])),rep(nlen,ny)),rep(seq(27.5,by=5,length=nlen),rep(ny,nlen)),
               type="n",ylab="Length bin",xlab="Year");

  symbols(rep(as.numeric(unlist(datao[,1])),rep(nlen,ny)),rep(seq(27.5,by=5,length=nlen),ny),
                circles=sqrt(as.vector(t(resp[,1:nlen]))),add=T,inches=0.1)
if(max(resp[,1:nlen])>.0000001){
                  symbols(rep(as.numeric(unlist(datao[,1])),rep(nlen,ny)),rep(seq(27.5,by=5,length=nlen),ny),
                circles=sqrt(-as.vector(t(resn[,1:nlen]))),add=T,
                            inches=0.1*(sqrt(max(-resn[,1:nlen]))/sqrt(max(resp[,1:nlen]))),bg=1)
}
x=rep(seq(27.5,by=5,length=nlen),ny)
x[as.vector(unlist(t(respstar)))<=maxres]=0
  points((rep(as.numeric(unlist(datao[,1])),rep(nlen,ny))),x,pch="*",cex=2.0)

 # title(sub=paste(" Standardized Pearson Residual Range ",round(min(rangeres),3),round(max(rangeres),3)),cex.sub=1.0,adj=0)
}

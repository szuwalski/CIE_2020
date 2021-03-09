library(dplyr)
library(ggplot2)
inYlim<-600
#==plot theme
.THEME    = theme_bw(base_size = 12, base_family = "") +
  theme(strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(color="white",fill="white")) 
.COL =  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

#==func for plot
ribbon_prep<-function(tmp,adj=1,Quantity_name="MMB",F_mort="NA")
{
  adj<-adj
  tmp<-apply(tmp,2,sort)
  med_ch<-as.numeric(tmp[nrow(tmp)/2,])/adj
  up_ch<-as.numeric(tmp[round(nrow(tmp)*0.025),])/adj
  dn_ch<-as.numeric(tmp[round(nrow(tmp)*0.975),])/adj
  years<-seq(from=min(colnames(tmp)),to=max(colnames(tmp)))
  out<-data.frame(MMB=med_ch,upper=up_ch,lower=dn_ch,Year=years,
                  Quantity=Quantity_name,F_mort=as.character(F_mort))
  return(out)
}

# organize data_2019
repfile<-scan("19_gmacs_proj_rec/mcoutSSB.REP")
hist_ssb<-matrix(repfile,ncol=37,byrow=T)
colnames(hist_ssb)<-seq(1982,2018)

innames<-c("Draw","Replicate","Treatment",paste("F",seq(1,length(M[[1]]$log_fbar))),
           "B35",paste("proj_ssb_",seq(2019,2025)))
projfile<-as.data.frame(matrix(scan("19_gmacs_proj_rec/mcoutPROJ.REP"),ncol=length(innames),byrow=TRUE))
colnames(projfile)<-innames

ind<-grep('ssb',colnames(projfile))
adj_proj<-filter(projfile,Replicate==1)
all_ssb<-cbind(hist_ssb,adj_proj[,ind])
colnames(all_ssb)<-seq(1982,2025)


# organize data_2020
repfile20<-scan("20_gmacs_proj_rec/mcoutSSB.REP")
hist_ssb20<-matrix(repfile20,ncol=length(seq(1982,2019)),byrow=T)
colnames(hist_ssb20)<-seq(1982,2019)

innames<-c("Draw","Replicate","Treatment",paste("F",seq(1,length(M[[1]]$log_fbar))),
           "B35",paste("proj_ssb_",seq(2020,2025)))
projfile20<-as.data.frame(matrix(scan("20_gmacs_proj_rec/mcoutPROJ.REP"),ncol=length(innames),byrow=TRUE))
colnames(projfile20)<-innames

ind<-grep('ssb',colnames(projfile20))
adj_proj20<-filter(projfile20,Replicate==1)
all_ssb20<-cbind(hist_ssb20,adj_proj20[,ind])
colnames(all_ssb20)<-seq(1982,2025)



#==manipulate 2019
ind<-grep('ssb',colnames(projfile))
adj_proj<-filter(projfile,Replicate==1 & Treatment==1)
all_ssb<-cbind(hist_ssb,adj_proj[,ind])
colnames(all_ssb)<-seq(1982,2025)
out<-ribbon_prep(tmp=all_ssb,F_mort=0)

adj_proj<-filter(projfile,Replicate==1 & Treatment==2)
all_ssb<-cbind(hist_ssb,adj_proj[,ind])
colnames(all_ssb)<-seq(1982,2025)
out2<-ribbon_prep(all_ssb,F_mort=1)

adj_proj<-filter(projfile,Replicate==1 & Treatment==3)
all_ssb<-cbind(hist_ssb,adj_proj[,ind])
colnames(all_ssb)<-seq(1982,2025)
out3<-ribbon_prep(all_ssb,F_mort=2)

in_dat_in<-as.data.frame(rbind(out,out2,out3))
in_dat<-filter(in_dat_in,F_mort==2)


#==manipulate 2020
ind<-grep('ssb',colnames(projfile20))
adj_proj<-filter(projfile20,Replicate==1 & Treatment==1)
all_ssb<-cbind(hist_ssb20,adj_proj[,ind])
colnames(all_ssb)<-seq(1982,2025)
out<-ribbon_prep(tmp=all_ssb,F_mort=0)

adj_proj<-filter(projfile20,Replicate==1 & Treatment==2)
all_ssb<-cbind(hist_ssb20,adj_proj[,ind])
colnames(all_ssb)<-seq(1982,2025)
out2<-ribbon_prep(all_ssb,F_mort=1)

adj_proj<-filter(projfile20,Replicate==1 & Treatment==3)
all_ssb<-cbind(hist_ssb20,adj_proj[,ind])
colnames(all_ssb)<-seq(1982,2025)
out3<-ribbon_prep(all_ssb,F_mort=2)

in_dat_in<-as.data.frame(rbind(out,out2,out3))
in_dat20<-filter(in_dat_in,F_mort==2)

all_dat<-rbind(in_dat,in_dat20)
all_dat$model<-rep(c("2019 catch","2020 catch"),each=nrow(in_dat))

both_proj<-ggplot() +
  geom_ribbon(data=all_dat,aes(x=Year,ymin = lower, ymax = upper, fill = model),alpha=0.1) +
  geom_line(data=all_dat,aes(x=Year,y=MMB,col=model)) +
  .THEME +
  theme(legend.position = c(.1,.85)) + 
  geom_hline(yintercept=M[[1]]$spr_bmsy,linetype=2) +
  expand_limits(y=0)+
  ylim(0,inYlim)
print(both_proj)



png("plots/projection_state.png",height=5,width=8,res=400,units='in')
print(both_proj)
dev.off()


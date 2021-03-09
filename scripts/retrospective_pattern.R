#==read in all of the retro runs via readADMB
#==plot the peels of 



#===PULL gmacs DATA AND outputs
mod_names <- seq(2018,2010)
.MODELDIR = c("./retro/2018/","./retro/2017/",
              "./retro/2016/","./retro/2015/",
              "./retro/2014/","./retro/2013/",
              "./retro/2012/","./retro/2011/","./retro/2010/")
              
# .MODELDIR = c("./retro_2014/2018/","./retro_2014/2017/",
#               "./retro_2014/2016/","./retro_2014/2015/",
#               "./retro_2014/2014/","./retro_2014/2013/",
#               "./retro_2014/2012/","./retro_2014/2011/","./retro_2014/2010/")
# 
# .MODELDIR = c("./retro_q/2018/","./retro_q/2017/",
#               "./retro_q/2016/","./retro_q/2015/",
#               "./retro_q/2014/","./retro_q/2013/",
#               "./retro_q/2012/","./retro_q/2011/")

fn       <- paste0(.MODELDIR, "gmacs")
M        <- lapply(fn, read_admb) #need .prj file to run gmacs and need .rep file here
names(M) <- mod_names

png("plots/retro_mmb.png",height=4.5,width=8,res=400,units='in')
par(mar=c(4,4,1,1))
plot(M[[1]]$ssb~M[[1]]$mod_yrs,type='l',ylim=c(0,max(M[[5]]$ssb)),las=1,lwd=2,ylab="MMB (1,000t)",xlab="Year")
for(x in 2:length(mod_names))
  lines(M[[x]]$ssb~M[[x]]$mod_yrs,lty=1,col=x,lwd=2)
dev.off()

in_name<-NULL
in_yr<-NULL
in_val<-NULL
in_peel<-NULL

for(x in 1:length(mod_names))
{
 in_name<-c(in_name,rep('ssb',length(M[[x]]$mod_yrs)),rep('rec',length(M[[x]]$mod_yrs)))
 in_peel<-c(in_peel,rep(mod_names[x],2*length(M[[x]]$mod_yrs)))
 in_val<-c(in_val,M[[x]]$ssb,unlist(M[[x]]$recruits[1,]))
 in_yr<-c(in_yr,M[[x]]$mod_yrs,M[[x]]$mod_yrs)
} 

out_data<-data.frame(quantity=in_name,peel=in_peel,value=in_val,year=in_yr)
write.csv(out_data,"opilio_tseries.csv")

B35<-NULL
for(x in 1:length(mod_names))
 B35<-c(B35,M[[x]]$spr_bmsy)

mod_names_2 <- seq(2018,2010)
.MODELDIR2 = c("./retro/2018_s/","./retro/2017_s/",
              "./retro/2016_s/","./retro/2015_s/",
              "./retro/2014_s/","./retro/2013_s/",
              "./retro/2012_s/","./retro/2011_s/","./retro/2010_s/")

# .MODELDIR2 = c("./retro_q/2018_s/","./retro_q/2017_s/",
#                "./retro_q/2016_s/","./retro_q/2015_s/",
#                "./retro_q/2014_s/","./retro_q/2013_s/",
#                "./retro_q/2012_s/","./retro_q/2011_s/")

fn       <- paste0(.MODELDIR2, "gmacs")
M2        <- lapply(fn, read_admb) #need .prj file to run gmacs and need .rep file here
names(M2) <- mod_names_2

B35_d<-NULL
for(x in 1:length(mod_names))
  B35_d<-c(B35_d,M2[[x]]$spr_bmsy)

BMSY<-data.frame(BMSY=c(B35,B35_d),peel=c(mod_names,mod_names_2),
                 survey=c(rep("whole",length(mod_names)),rep("drop",length(mod_names_2))))
write.csv(BMSY,"opilio_bmsy.csv")

in_name<-NULL
in_yr<-NULL
in_val<-NULL
in_peel<-NULL

for(x in 1:length(mod_names_2))
{
  in_name<-c(in_name,rep('ssb',length(M2[[x]]$mod_yrs)),rep('rec',length(M2[[x]]$mod_yrs)))
  in_peel<-c(in_peel,rep(mod_names_2[x],2*length(M2[[x]]$mod_yrs)))
  in_val<-c(in_val,M2[[x]]$ssb,unlist(M2[[x]]$recruits[1,]))
  in_yr<-c(in_yr,M2[[x]]$mod_yrs,M2[[x]]$mod_yrs)
} 

out_data<-data.frame(quantity=in_name,peel=in_peel,value=in_val,year=in_yr)
write.csv(out_data,"opilio_tseries_drop.csv")



OFL<-rep(NA,length(mod_names))
OFL2<-rep(NA,length(mod_names))


for(x in 1:length(OFL))
{
  OFL[x]<-M[[x]]$spr_cofl
  OFL2[x]<-M2[[x]]$spr_cofl
}

png("plots/retro_OFL.png",height=4.5,width=8,res=400,units='in')
plot(OFL~as.numeric(mod_names+1),type='l',lwd=2,las=1,xlab="Year",ylim=c(0,250),ylab="OFL (1,000t)")
lines(OFL2~as.numeric(mod_names+1),lty=2,col=2,lwd=2)
legend('bottomleft',bty='n',col=c(1,2),lty=c(1,2),lwd=2,
       legend=c("With last year of survey","Without last year of survey"))
dev.off()

ofl_diff<-(OFL2-OFL)/OFL
mean(ofl_diff)
median(ofl_diff)

boxplot(ofl_diff*100,las=1,frame=F)
abline(h=0,lty=2)
mtext(side=2,line=2.2,"Percent change in OFL")
 df<-data.frame(Change=ofl_diff*100,stock="Opilio")

p <- ggplot(data = df, aes(y = Change, x = stock)) + 
  geom_boxplot(aes(middle = mean(Change)))

print(p)

#==martin's figure
retro_ssb<-rep(NA,length(mod_names))
retro_ssb_no<-rep(NA,length(mod_names))

for(x in 1:length(mod_names))
{
  retro_ssb[x]<-M[[x]]$ssb[length(M[[x]]$ssb)]
  retro_ssb_no[x]<-M2[[x]]$ssb[length(M[[x]]$ssb)]
  }

png("plots/retro_term_yr_mmb.png",height=4.5,width=8,res=400,units='in')
plot(M[[1]]$ssb~M[[1]]$mod_yrs,type='b',ylim=c(0,320),xlim=c(2011,2018),
     las=1,lwd=2,ylab="MMB (1,000t)",xlab="Year")
lines(retro_ssb~mod_names,lty=1,col=2,lwd=2,type='b')
lines(retro_ssb_no~mod_names,lty=1,col=3,lwd=2,type='b')

legend('topright',bty='n',col=c(1,2,3),lty=1,lwd=2,legend=c("Most recent assessment","Standard retrospective",
                                                "Drop terminal survey retrospective"))
dev.off()

ref_ssb<-rev(M[[1]]$ssb[(length(M[[1]]$ssb)-length(mod_names)+1):length(M[[1]]$ssb)])
mean((retro_ssb-ref_ssb)/ref_ssb)
mean((retro_ssb_no-ref_ssb)/ref_ssb,na.rm=T)

#==ralston's sigma

sqrt((1/(length(ref_ssb)-1))*sum((log(retro_ssb)-log(ref_ssb))^2))
sqrt((1/(length(ref_ssb)-1))*sum((log(retro_ssb_no)-log(ref_ssb))^2))

sqrt(sum(log(retro_ssb_no)-log(ref_ssb))/(length(ref_ssb)-1))

#==martin's figure
retro_q<-rep(NA,length(mod_names))
retro_q_no<-rep(NA,length(mod_names))

retro_m<-rep(NA,length(mod_names))
retro_m_no<-rep(NA,length(mod_names))

for(x in 1:length(mod_names))
{
  retro_q[x]<-M[[x]]$survey_q[4]
  retro_q_no[x]<-M2[[x]]$survey_q[4]
  
  retro_m[x]<-unique(c(M[[x]]$M))[1]
  retro_m_no[x]<-unique(c(M2[[x]]$M))[1]
}

png("plots/retro_q_m.png",height=4.5,width=8,res=400,units='in')

plot(retro_q~mod_names,ylim=c(0,1),type='l',las=1,pch=15,col=3,lty=2,ylab='',xlab='')
lines(retro_q_no~mod_names,col=3,pch=15,lty=1)

lines(retro_m~mod_names,col=4,pch=16,lty=2,type='l')
lines(retro_m_no~mod_names,col=4,pch=16,lty=1)

mtext(side=1,"Year",line=2.5)
legend('bottom',bty='n',col=c(3,4),lty=c(2,1),legend=c("Terminal survey include","Terminal survey excluded"))

dev.off()


lines(retro_m~mod_names,type='b',col=1)
lines(retro_m_no~mod_names,type='b',col=2)
projyr<-c(mod_names+1,mod_names_2+1)
outs<-data.frame(Spp="Opilio",projYear=projyr,retro=c(rep('retr',length(mod_names)),rep('survRed_retr',length(mod_names))),
           peel=projyr-max(projyr),SSB=c(retro_ssb,retro_ssb_no),OFL=c(OFL,OFL2))
#write.csv(outs,"opilio_retro.csv")


png("plots/retro_mmb_2.png",height=8,width=8,res=400,units='in')
par(mfrow=c(2,1),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
plot(M[[1]]$ssb~M[[1]]$mod_yrs,type='l',ylim=c(0,max(M[[1]]$ssb,M2[[1]]$ssb)),las=1,lwd=2,ylab="MMB (1,000t)",xlab="Year",xaxt='n')
for(x in 2:length(mod_names))
  lines(M[[x]]$ssb~M[[x]]$mod_yrs,lty=1,col=x,lwd=2)
legend('topright',bty='n',legend=c("Terminal survey included",
                                   paste("Mohn's rho = ",
                                         round(mean((retro_ssb-ref_ssb)/ref_ssb),2))))

mean((retro_ssb-ref_ssb)/ref_ssb)
mean((retro_ssb_no-ref_ssb)/ref_ssb)


plot(M2[[1]]$ssb~M2[[1]]$mod_yrs,type='l',ylim=c(0,max(M[[1]]$ssb,M2[[1]]$ssb)),las=1,lwd=2,ylab="MMB (1,000t)",xlab="Year")
for(x in 2:length(mod_names))
  lines(M2[[x]]$ssb~M2[[x]]$mod_yrs,lty=1,col=x,lwd=2)
legend('topright',bty='n',legend=c("Terminal survey excluded",
                                   paste("Mohn's rho = ",
                                         round(mean((retro_ssb_no-ref_ssb)/ref_ssb,na.rm=T),2))))
mtext(side=2,outer=T,"MMB (1,000t)",line=2.5)
dev.off()

#==pull 
      
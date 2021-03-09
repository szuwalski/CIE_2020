library(dplyr)
#==MARE for growth, catch, cpue
#==MAE for size composition data

#===========================
# GMACS
#===========================
#==catch
catch_mare<-matrix(nrow=nrow(M[[1]]$obs_catch),ncol=length(M))
colnames(catch_mare)<-mod_names
rownames(catch_mare)<-c("Retained","Discard (male)","Discard (female)","Bycatch")
for(x in 1:ncol(catch_mare))
  catch_mare[,x]<-round(apply(abs((M[[x]]$pre_catch - M[[x]]$obs_catch) / M[[x]]$obs_catch),1,mean,na.rm=T),3)


#==indices of abundance
cpue_mare<-matrix(ncol=length(M),nrow=12)
colnames(cpue_mare)<-mod_names
rownames(cpue_mare)<-c("Survey MMB era 1","Survey MMB era 2",
                       "2009 BSFRF MMB","2009 NMFS MMB",
                       "2010 BSFRF MMB","2010 NMFS MMB",
                       "Survey FMB era 1","Survey FMB era 2",
                       "2009 BSFRF FMB","2009 NMFS FMB",
                       "2010 BSFRF FMB","2010 NMFS FMB")

for(x in 1:length(M))
{
  data<-data.frame(ARE=abs((M[[x]]$pre_cpue - M[[x]]$obs_cpue)/M[[x]]$obs_cpue),group=M[[x]]$dSurveyData[,4],sex=M[[x]]$dSurveyData[,5])
  tmp<-data %>% 
  group_by(sex,group) %>% 
  summarize(MARE = mean(ARE))
  cpue_mare[,x]<-round(tmp$MARE,3)
}





#===========================
# Status quo
#===========================
catch_mare_sq<-matrix(nrow=nrow(M[[1]]$obs_catch),ncol=length(snowad.rep))
colnames(catch_mare_sq)<-c("19.1","20.1")
rownames(catch_mare_sq)<-c("Retained","Discard (male)","Discard (female)","Bycatch")
for(x in 1:ncol(catch_mare_sq))
{
  obs<-snowad.rep[[x]]$"Observed retained catch biomass"
  pred<-snowad.rep[[x]]$"Predicted retained catch biomas"
  catch_mare_sq[1,x]<-round(mean(abs((pred - obs) / obs)),3)
  
  obs<-snowad.rep[[x]]$"observed male discard mortality biomass"
  pred<-snowad.rep[[x]]$"predicted discard male catch biomass"
  catch_mare_sq[2,x]<-round(mean(abs((pred - obs) / obs)),3)
  
  obs<-snowad.rep[[x]]$"observed female discard mortality biomass"
  pred<-snowad.rep[[x]]$"predicted female discard mortality biomass"
  catch_mare_sq[3,x]<-round(mean(abs((pred - obs) / obs)),3)
  
  obs<-snowad.rep[[x]]$"observed trawl catch biomass"
  obs<-obs[obs!=0]
  pred<-snowad.rep[[x]]$"predicted trawl catch biomass"
  pred<-pred[pred!=0]
  pred<-pred[-length(pred)]
  catch_mare_sq[4,x]<-round(mean(abs((pred - obs) / obs)),3)
  
}



#==indices of abundance
cpue_mare_sq<-matrix(ncol=length(snowad.rep),nrow=12)
colnames(cpue_mare_sq)<-"SQ"
rownames(cpue_mare_sq)<-c("Survey MMB era 1","Survey MMB era 2",
                       "2009 BSFRF MMB","2009 NMFS MMB",
                       "2010 BSFRF MMB","2010 NMFS MMB",
                       "Survey FMB era 1","Survey FMB era 2",
                       "2009 BSFRF FMB","2009 NMFS FMB",
                       "2010 BSFRF FMB","2010 NMFS FMB")

for(x in 1:ncol(catch_mare_sq))
{
  obs<-snowad.rep[[x]]$"Observed survey male spawning biomass"[1:7]
  pred<-snowad.rep[[x]]$"Predicted Male survey mature Biomass"[1:7]
  cpue_mare_sq[1,x]<-round(mean(abs((pred - obs) / obs)),3)
  
  obs<-snowad.rep[[x]]$"Observed survey male spawning biomass"[8:38]
  pred<-snowad.rep[[x]]$"Predicted Male survey mature Biomass"[8:38]
  cpue_mare_sq[2,x]<-round(mean(abs((pred - obs) / obs)),3)

  obs<-snowad.rep[[x]]$"Observed survey female spawning biomass"[1:7]
  pred<-snowad.rep[[x]]$"Predicted Female survey mature Biomass"[1:7]
  cpue_mare_sq[7,x]<-round(mean(abs((pred - obs) / obs)),3)
  
  obs<-snowad.rep[[x]]$"Observed survey female spawning biomass"[8:38]
  pred<-snowad.rep[[x]]$"Predicted Female survey mature Biomass"[8:38]
  cpue_mare_sq[8,x]<-round(mean(abs((pred - obs) / obs)),3)
  
  
  #BSFRF
  obs<-snowad.rep[[x]]$"Observed industry survey mature biomass"[1]
  pred<-snowad.rep[[x]]$"Predicted industry survey mature biomass"[1]
  cpue_mare_sq[3,x]<-round(mean(abs((pred - obs) / obs)),3)
  
  obs<-snowad.rep[[x]]$"Observed industry survey mature biomass"[2]
  pred<-snowad.rep[[x]]$"Predicted industry survey mature biomass"[2]
  cpue_mare_sq[4,x]<-round(mean(abs((pred - obs) / obs)),3)
  
  obs<-snowad.rep[[x]]$"Observed industry survey mature biomass"[3]
  pred<-snowad.rep[[x]]$"Predicted industry survey mature biomass"[3]
  cpue_mare_sq[5,x]<-round(mean(abs((pred - obs) / obs)),3)
  
  obs<-snowad.rep[[x]]$"Observed industry survey mature biomass"[4]
  pred<-snowad.rep[[x]]$"Predicted industry survey mature biomass"[4]
  cpue_mare_sq[6,x]<-round(mean(abs((pred - obs) / obs)),3)
  
  #BSFRF 2010
  obs<-snowad.rep[[x]]$"Observed 2010 industry survey mature biomass"[1]
  pred<-snowad.rep[[x]]$"Predicted 2010 industry survey mature biomass"[1]
  cpue_mare_sq[9,x]<-round(mean(abs((pred - obs) / obs)),3)
  
  obs<-snowad.rep[[x]]$"Observed 2010 industry survey mature biomass"[2]
  pred<-snowad.rep[[x]]$"Predicted 2010 industry survey mature biomass"[2]
  cpue_mare_sq[10,x]<-round(mean(abs((pred - obs) / obs)),3)
  
  obs<-snowad.rep[[x]]$"Observed 2010 industry survey mature biomass"[3]
  pred<-snowad.rep[[x]]$"Predicted 2010 industry survey mature biomass"[3]
  cpue_mare_sq[11,x]<-round(mean(abs((pred - obs) / obs)),3)
  
  obs<-snowad.rep[[x]]$"Observed 2010 industry survey mature biomass"[4]
  pred<-snowad.rep[[x]]$"Predicted 2010 industry survey mature biomass"[4]
  cpue_mare_sq[12,x]<-round(mean(abs((pred - obs) / obs)),3)

}  

  
  
  
#==========================================
# Plots for catch and surveys
#+==========================================

in_mat_1<-rbind(catch_mare,cpue_mare)
in_mat_2<-rbind(catch_mare_sq,cpue_mare_sq)

in_mat<-cbind(in_mat_2,in_mat_1)

in_mat<-in_mat[c(1,2,3,4,5,6,11,12,7,8,9,10,13,14,15,16),]
in_mat<-in_mat[,-1]
library(RColorBrewer)
useCol<-matrix(as.numeric(in_mat),ncol=ncol(in_mat))
change_cex<-0.8
for(x in 1:ncol(useCol))
{
    colrange <-seq(0,.35,0.001)
    col 	<-colorRampPalette(c('white', brewer.pal(3,"Set1")[2]))(length(colrange))
    
    for(i in 1:nrow(useCol))
      useCol[i,x] <- col[which(abs(colrange-in_mat[i,x]) == min(abs(colrange-in_mat[i,x])))][1] 
}
png("plots/cpue_MARE.png",height=8,width=8,res=400,units='in')
library(plotrix)
par(mar=c(4,10,3,1))
color2D.matplot(in_mat, 
                show.values = 2,
                axes = FALSE,
                xlab = "",
                ylab='',
                vcex = 0.8,
                vcol = "black",
                cellcolors = useCol ,las=1)
axis(side=2,at=seq(0.5,nrow(in_mat)),labels=rev(rownames(in_mat)),las=1)
axis(side=3,at=seq(0.5,ncol(in_mat)),labels=(colnames(in_mat)),las=1)
dev.off()

###################################################
# Size composition data GMACS

size_comp_MAE<-matrix(ncol=length(M),nrow=16)
colnames(size_comp_MAE)<-mod_names
for(x in 1:length(M))
{
  
  tmp<-abs(M[[x]]$d3_pre_size_comps_out-M[[x]]$d3_obs_size_comps_out)
  tmp2<-cbind(M[[x]]$d3_SizeComps_in[,3],M[[x]]$d3_SizeComps_in[,4],as.numeric(tmp))
  data<-data.frame(tmp2)
  colnames(data)<-c("fleet","sex",'AE')
  tmp<-data %>% 
    group_by(sex,fleet) %>% 
    summarize(MAE = mean(AE))
  size_comp_MAE[,x]<-round(tmp$MAE,3)
}
rownames(size_comp_MAE)<-c("Directed male","Trawl male","NMFS (1982-88) male","NMFS (1989-present) male",
                           "BSFRF 2009 male","NMFS 2009 male","BSFRF 2010 male","NMFS 2010 male",
                           "Directed female","Trawl female","NMFS (1982-88) female","NMFS (1989-present) female",
                           "BSFRF 2009 female","NMFS 2009 female","BSFRF 2010 female","NMFS 2010 female")

#==indices of abundance
size_comp_MAE_sq<-matrix(ncol=length(snowad.rep),nrow=16)
colnames(size_comp_MAE_sq)<-c("19.1","20.1")
rownames(size_comp_MAE_sq)<-c("Directed male","Trawl male","NMFS (1982-88) male","NMFS (1989-present) male",
                           "BSFRF 2009 male","NMFS 2009 male","BSFRF 2010 male","NMFS 2010 male",
                           "Directed female","Trawl female","NMFS (1982-88) female","NMFS (1989-present) female",
                           "BSFRF 2009 female","NMFS 2009 female","BSFRF 2010 female","NMFS 2010 female")

make_sum_1<-function(input)
{
  tmp<-input
  tots<-sum(tmp)
  obs<-tmp/tots
  obs[obs="NaN"]<-0
  return(obs)
}


for(x in 1:ncol(size_comp_MAE_sq))
{
   #males
  obs1<-(snowad.rep[[x]]$"Observed proportion fishery retained new male" +snowad.rep[[x]]$"Observed proportion fishery retained old male")[,-1]
  pred1<-(snowad.rep[[x]]$"Predicted proportion fishery retained new male" + snowad.rep[[x]]$"Predicted proportion fishery retained old male")[,-1]
  pred2<-(snowad.rep[[x]]$"Predicted proportion fishery total new male" + snowad.rep[[x]]$"Predicted proportion fishery total old male")[,-1]
  obs2<-(snowad.rep[[x]]$"Observed proportion fishery total new male" +snowad.rep[[x]]$"Observed proportion fishery total old male")[,-1]
  obs<-c(t(apply(obs1,1,make_sum_1)),t(apply(obs2,1,make_sum_1)))
  pred<-c(t(apply(pred1,1,make_sum_1)),t(apply(pred2,1,make_sum_1)))
  
  
  size_comp_MAE_sq[1,x]<-round(mean((abs(pred - obs) )),3)
  
   obs<-snowad.rep[[x]]$"Observed proportion trawl male"[,-1]
   pred<-snowad.rep[[x]]$"Predicted proportion trawl male" [,-1]
   obs<-t(apply(obs,1,make_sum_1))
   pred<-t(apply(pred,1,make_sum_1))
   size_comp_MAE_sq[2,x]<-round(mean((abs(pred - obs) )),3)  

   obs1<-snowad.rep[[x]]$"Observed proportion survey immature new male" [1:7,-1] + 
     snowad.rep[[x]]$"Observed proportion survey immature old male" [1:7,-1]
   pred1<-snowad.rep[[x]]$"Predicted proportion survey immature new male" [1:7,-1] + 
     snowad.rep[[x]]$"Predicted proportion survey immature old male" [1:7,-1]
   obs2<-snowad.rep[[x]]$"Observed proportion survey mature new male" [1:7,-1] + 
     snowad.rep[[x]]$"Observed proportion survey mature old male" [1:7,-1]
   pred2<-snowad.rep[[x]]$"Predicted proportion survey mature new male" [1:7,-1] + 
     snowad.rep[[x]]$"Predicted proportion survey mature old male" [1:7,-1]
   obs<-c(t(apply(obs1,1,make_sum_1)),t(apply(obs2,1,make_sum_1)))
   pred<-c(t(apply(pred1,1,make_sum_1)),t(apply(pred2,1,make_sum_1)))
   
   size_comp_MAE_sq[3,x]<-round(mean((abs(pred - obs) )),3)  

   obs1<-snowad.rep[[x]]$"Observed proportion survey immature new male" [8:38,-1] + 
     snowad.rep[[x]]$"Observed proportion survey immature old male" [8:38,-1]
   pred1<-snowad.rep[[x]]$"Predicted proportion survey immature new male" [8:38,-1] + 
     snowad.rep[[x]]$"Predicted proportion survey immature old male" [8:38,-1]
   obs2<-snowad.rep[[x]]$"Observed proportion survey mature new male" [8:38,-1] + 
     snowad.rep[[x]]$"Observed proportion survey mature old male" [8:38,-1]
   pred2<-snowad.rep[[x]]$"Predicted proportion survey mature new male" [8:38,-1] + 
     snowad.rep[[x]]$"Predicted proportion survey mature old male" [8:38,-1]
   obs<-c(t(apply(obs1,1,make_sum_1)),t(apply(obs2,1,make_sum_1)))
   pred<-c(t(apply(pred1,1,make_sum_1)),t(apply(pred2,1,make_sum_1)))
   
   size_comp_MAE_sq[4,x]<-round(mean((abs(pred - obs) )),3)  
   
   
   #females
   obs<-snowad.rep[[x]]$"Observed proportion fishery discard all female"[,-1]
   pred<-snowad.rep[[x]]$"Predicted proportion fishery discard all female" [,-1]
   obs<-t(apply(obs,1,make_sum_1))
   pred<-t(apply(pred,1,make_sum_1))
   size_comp_MAE_sq[9,x]<-round(mean((abs(pred - obs) )),3)
   
   obs<-snowad.rep[[x]]$"Observed proportion trawl female"[,-1]
   pred<-snowad.rep[[x]]$"Predicted proportion trawl female" [,-1]
   obs<-t(apply(obs,1,make_sum_1))
   pred<-t(apply(pred,1,make_sum_1))
   size_comp_MAE_sq[10,x]<-round(mean((abs(pred - obs) )),3)  
   
   obs1<-snowad.rep[[x]]$"Observed proportion survey immature new female" [1:7,-1] + 
     snowad.rep[[x]]$"Observed proportion survey immature old female" [1:7,-1]
   pred1<-snowad.rep[[x]]$"Predicted proportion survey immature new female" [1:7,-1] + 
     snowad.rep[[x]]$"Predicted proportion survey immature old female" [1:7,-1]
   obs2<-snowad.rep[[x]]$"Observed proportion survey mature new female" [1:7,-1] + 
     snowad.rep[[x]]$"Observed proportion survey mature old female" [1:7,-1]
   pred2<-snowad.rep[[x]]$"Predicted proportion survey mature new female" [1:7,-1] + 
     snowad.rep[[x]]$"Predicted proportion survey mature old female" [1:7,-1]
   obs<-c(t(apply(obs1,1,make_sum_1)),t(apply(obs2,1,make_sum_1)))
   pred<-c(t(apply(pred1,1,make_sum_1)),t(apply(pred2,1,make_sum_1)))
   
   size_comp_MAE_sq[11,x]<-round(mean((abs(pred - obs) )),3)  
   
   obs1<-snowad.rep[[x]]$"Observed proportion survey immature new female" [8:38,-1] + 
     snowad.rep[[x]]$"Observed proportion survey immature old female" [8:38,-1]
   pred1<-snowad.rep[[x]]$"Predicted proportion survey immature new female" [8:38,-1] + 
     snowad.rep[[x]]$"Predicted proportion survey immature old female" [8:38,-1]
   obs2<-snowad.rep[[x]]$"Observed proportion survey mature new female" [8:38,-1] + 
     snowad.rep[[x]]$"Observed proportion survey mature old female" [8:38,-1]
   pred2<-snowad.rep[[x]]$"Predicted proportion survey mature new female" [8:38,-1] + 
     snowad.rep[[x]]$"Predicted proportion survey mature old female" [8:38,-1]
   obs<-c(t(apply(obs1,1,make_sum_1)),t(apply(obs2,1,make_sum_1)))
   pred<-c(t(apply(pred1,1,make_sum_1)),t(apply(pred2,1,make_sum_1)))
   
   size_comp_MAE_sq[12,x]<-round(mean((abs(pred - obs) )),3)  
   
   #BSFRF data
   obs<-snowad.rep[[x]]$"Observed Prop industry survey male"
   pred<-snowad.rep[[x]]$"Predicted Prop industry survey male"
   obs<-t(make_sum_1(obs))
   pred<-t(make_sum_1(pred))
   size_comp_MAE_sq[5,x]<-round(mean((abs(pred - obs) )),3)  
   
   obs<-snowad.rep[[x]]$"Observed Prop industry nmfs survey male"
   pred<-snowad.rep[[x]]$"redicted Prop industry nmfs survey male"
   obs<-t(make_sum_1(obs))
   pred<-t(make_sum_1(pred))
   size_comp_MAE_sq[6,x]<-round(mean((abs(pred - obs) )),3)  
   
   obs<-snowad.rep[[x]]$"Observed Prop 2010 industry survey male"
   pred<-snowad.rep[[x]]$"Predicted Prop 2010 industry survey male"
   obs<-t(make_sum_1(obs))
   pred<-t(make_sum_1(pred))
   size_comp_MAE_sq[7,x]<-round(mean((abs(pred - obs) )),3)  
   
   obs<-snowad.rep[[x]]$"Observed Prop 2010 industry nmfs survey male"
   pred<-snowad.rep[[x]]$"Predicted Prop 2010 industry nmfs survey male"
   obs<-t(make_sum_1(obs))
   pred<-t(make_sum_1(pred))
   size_comp_MAE_sq[8,x]<-round(mean((abs(pred - obs) )),3) 
   
   #females
   obs<-snowad.rep[[x]]$"bserved Prop industry survey female"
   pred<-snowad.rep[[x]]$"Predicted Prop industry survey female"
   obs<-t(make_sum_1(obs))
   pred<-t(make_sum_1(pred))
   size_comp_MAE_sq[13,x]<-round(mean((abs(pred - obs) )),3)  
   
   obs<-snowad.rep[[x]]$"Observed Prop industry nmfs survey female"
   pred<-snowad.rep[[x]]$"Predicted Prop industry nmfs survey female"
   obs<-t(make_sum_1(obs))
   pred<-t(make_sum_1(pred))
   size_comp_MAE_sq[14,x]<-round(mean((abs(pred - obs) )),3)  
   
   obs<-snowad.rep[[x]]$"Observed Prop 2010 industry survey female"
   pred<-snowad.rep[[x]]$"Predicted Prop 2010 industry survey female"
   obs<-t(make_sum_1(obs))
   pred<-t(make_sum_1(pred))
   size_comp_MAE_sq[15,x]<-round(mean((abs(pred - obs) )),3)  
   
   obs<-snowad.rep[[x]]$"Observed Prop 2010 industry nmfs survey female"
   pred<-snowad.rep[[x]]$"Predicted Prop 2010 industry nmfs survey female"
   obs<-t(make_sum_1(obs))
   pred<-t(make_sum_1(pred))
   size_comp_MAE_sq[16,x]<-round(mean((abs(pred - obs) )),3) 
}
   
   
   
in_mat<-cbind(size_comp_MAE_sq,size_comp_MAE)

in_mat<-in_mat[c(1,2,9,10,3,4,5,6,7,8,11,12,13,14,15,16),]
in_mat<-in_mat[,-1]
library(RColorBrewer)
useCol<-matrix(as.numeric(in_mat),ncol=ncol(in_mat))
change_cex<-0.8
for(x in 1:ncol(useCol))
{
  colrange <-seq(0,.025,0.001)
  col 	<-colorRampPalette(c('white', brewer.pal(3,"Set1")[2]))(length(colrange))
  
  for(i in 1:nrow(useCol))
    useCol[i,x] <- col[which(abs(colrange-in_mat[i,x]) == min(abs(colrange-in_mat[i,x])))][1] 
}

png("plots/size_comp_MAE.png",height=8,width=8,res=400,units='in')
library(plotrix)
par(mar=c(4,12,3,1))
color2D.matplot(in_mat, 
                show.values = 3,
                axes = FALSE,
                xlab = "",
                ylab='',
                vcex = 0.8,
                vcol = "black",
                cellcolors = useCol ,las=1)
axis(side=2,at=seq(0.5,nrow(in_mat)),labels=rev(rownames(in_mat)),las=1)
axis(side=3,at=seq(0.5,ncol(in_mat)),labels=(colnames(in_mat)),las=1)
dev.off()






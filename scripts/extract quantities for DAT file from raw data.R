library(dplyr)

PlotLenComp<-function(Input,Input2="NA")
{
par(mfrow=c(7,7),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
for(x in 1:nrow(Input))
 {
 plot(Input[x,],xaxt='n',yaxt='n',type='l')
 if(Input2!="NA")
  lines(Input2[x,],col=2,lty=2)
}
}

#=============================================
# Total survey numbers male plus female
#=============================================
EBSdata	<-read.csv("new data 2019/EBSCrab_AB_Summary.csv")
names(EBSdata)[1]<-"SURVEY_YEAR"

SurveyTotals<- EBSdata%>%
 		group_by(SURVEY_YEAR) %>%
  		summarise(Totals = sum(ABUNDANCE))

plot(SurveyTotals$Totals~SurveyTotals$SURVEY_YEAR,las=1,type="l")
write.csv(SurveyTotals$Totals,"survey_totalN.csv")

#========================
#==Female Survey Length compositions==
#========================
EBSdata_in<-read.csv("new data 2019/EBSCrab_Abundance_Biomass.csv")
EBSdata<-EBSdata_in[EBSdata_in$SIZE_CLASS_MM>=25 & EBSdata_in$SEX=="FEMALE",]

Years				<-sort(unique(EBSdata$SURVEY_YEAR))
LengthBins			<-seq(25,135,5)
FemaleNewMature		<-matrix(nrow=length(Years),ncol=length(LengthBins))
FemaleOldMature		<-matrix(nrow=length(Years),ncol=length(LengthBins))
FemaleNewImmature		<-matrix(nrow=length(Years),ncol=length(LengthBins))
FemaleOldImmature		<-matrix(0,nrow=length(Years),ncol=length(LengthBins))

NewShellIndex<-3
head(EBSdata)
for(x in 1:length(Years))
{
 temp<-EBSdata[EBSdata$SURVEY_YEAR==Years[x],]
 for(y in 1:(length(LengthBins)-1))
 {
 FemaleNewMature[x,y]	<-sum(temp$ABUNDANCE[temp$SHELL_CONDITION<NewShellIndex & temp$MATURITY=="MATURE" & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
 FemaleOldMature[x,y]	<-sum(temp$ABUNDANCE[temp$SHELL_CONDITION>=NewShellIndex & temp$MATURITY=="MATURE" & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
 FemaleNewImmature[x,y]	<-sum(temp$ABUNDANCE[temp$SHELL_CONDITION<NewShellIndex & temp$MATURITY=="IMMATURE" & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
 #FemaleOldImmature[x,y]	<-sum(temp$ABUNDANCE[temp$SHELL_CONDITION>=NewShellIndex & temp$MATURITY=="IMMATURE" & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
}
}

FemaleNewMature<-FemaleNewMature/1000000
FemaleOldMature<-FemaleOldMature/1000000
FemaleNewImmature<-FemaleNewImmature/1000000


write.table(round(FemaleNewMature,4)[,1:(length(LengthBins)-1)], "FemNewMatSurv.txt",row.names=FALSE,col.names=F)
write.table(round(FemaleOldMature,4)[,1:(length(LengthBins)-1)], "FemOlMatSurv.txt",row.names=FALSE,col.names=F)
write.table(round(FemaleNewImmature,4)[,1:(length(LengthBins)-1)], "FemNewImmSurv.txt",row.names=FALSE,col.names=F)
write.table(round(FemaleOldImmature,4)[,1:(length(LengthBins)-1)], "FemOldImmSurv.txt",row.names=FALSE,col.names=F)

PlotLenComp(Input=FemaleNewMature)
PlotLenComp(Input=FemaleOldMature)
PlotLenComp(Input=FemaleNewImmature)

#========================
#==Male survey Length compositions==
#========================
EBSdata_in<-read.csv("new data 2019/EBSCrab_Abundance_Biomass_m.csv")
EBSdata<-EBSdata_in[EBSdata_in$SIZE_CLASS_MM>=25 & EBSdata_in$SEX=="MALE",]
MaturityData<-as.numeric(unlist(read.csv("new data 2019/ProbabilityMature.csv",header=F)))

Years				<-sort(unique(EBSdata$SURVEY_YEAR))
LengthBins			<-seq(25,155,5)
MaleNew			<-matrix(nrow=length(Years),ncol=length(LengthBins))
MaleOld			<-matrix(nrow=length(Years),ncol=length(LengthBins))
MaleNewMature		<-matrix(nrow=length(Years),ncol=length(LengthBins))
MaleOldMature		<-matrix(nrow=length(Years),ncol=length(LengthBins))
MaleNewImmature		<-matrix(nrow=length(Years),ncol=length(LengthBins))
MaleOldImmature		<-matrix(0,nrow=length(Years),ncol=length(LengthBins))
NewShellIndex		<-3

for(x in 1:length(Years))
{
 temp<-EBSdata[EBSdata$SURVEY_YEAR==Years[x],]
 for(y in 1:(length(LengthBins)-1))
 {
 MaleNew[x,y]	<-sum(temp$ABUNDANCE[temp$SHELL_CONDITION<NewShellIndex & temp$SEX=="MALE"  & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
 MaleOld[x,y]	<-sum(temp$ABUNDANCE[temp$SHELL_CONDITION>=NewShellIndex & temp$SEX=="MALE"  & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
}
}

#==cap the sizes
TakeInd<-which(LengthBins==135)-1
MaleNew[,TakeInd]<-apply(MaleNew[,TakeInd:length(LengthBins)],1,sum,na.rm=T)
MaleOld[,TakeInd]<-apply(MaleOld[,TakeInd:length(LengthBins)],1,sum,na.rm=T)
MaleNew<-MaleNew[,1:TakeInd]/1000000	
MaleOld<-MaleOld[,1:TakeInd]/1000000	

#==divide out to mature and immature
MaleNewMature 	<-sweep(MaleNew,2,MaturityData,FUN="*")
MaleNewImmature   <-sweep(MaleNew,2,1-MaturityData,FUN="*")

PlotLenComp(MaleNewMature)
PlotLenComp(MaleNewImmature)	
PlotLenComp(MaleOld)

write.table(round(MaleNewMature ,4), "MalNewMatSurv.txt",row.names=FALSE,col.names=F)
write.table(round(MaleOld,4), "MalOldMatSurv.txt",row.names=FALSE,col.names=F)
write.table(round(MaleNewImmature,4), "MalNewImmSurv.txt",row.names=FALSE,col.names=F)

#==make numbers if different survey bins are used
cutoff<-3
fem_n_m_tot<-apply(FemaleNewMature[,cutoff:LengthBinN],1,sum,na.rm=T)
fem_o_m_tot<-apply(FemaleOldMature[,cutoff:LengthBinN],1,sum,na.rm=T)
fem_n_i_tot<-apply(FemaleNewImmature[,cutoff:LengthBinN],1,sum,na.rm=T)
mal_n_m_tot<-apply(MaleNewMature[,cutoff:LengthBinN],1,sum,na.rm=T)
mal_o_m_tot<-apply(MaleOld[,cutoff:LengthBinN],1,sum,na.rm=T)
mal_n_i_tot<-apply(MaleNewImmature[,cutoff:LengthBinN],1,sum,na.rm=T)

big_tot<-cbind(mal_n_m_tot,mal_o_m_tot,mal_n_i_tot,fem_n_m_tot,fem_o_m_tot,fem_n_i_tot)

plot(apply(big_tot,1,sum),type='l')
lines(SurveyTotals$Totals/1000000,col=2)
#write.table(apply(big_tot,1,sum)[3:38],"total_survey_numbers_trimmed.txt",row.names=F,col.names=F)

(apply(big_tot,1,sum)-SurveyTotals$Totals/1000000)/(SurveyTotals$Totals/1000000)

#============================================================
#  Retained catch length frequencies
# modified from "Item 1_20XX-20XX snow crab dockside
#============================================================
Data<-read.csv("data/item1_dockside_size_comp.csv")
binwidth<-2.5
LengthBins<-seq(2.5,162.5,binwidth*2)
LengthBins<-seq(27.5,165,5)
TotalNnew<-rep(0,length(LengthBins))
TotalNold<-rep(0,length(LengthBins))

upperBnd<-132.5
lowerBnd<-27.5

for(x in 1:length(LengthBins))
{
 tmp<-Data$new[Data$size>(LengthBins[x]-binwidth) & Data$size<=(LengthBins[x]+binwidth)]
 if(length(tmp)>0)
  TotalNnew[x]<-sum(tmp)
 tmp<-Data$old[Data$size>(LengthBins[x]-binwidth) & Data$size<=(LengthBins[x]+binwidth)]
 if(length(tmp)>0)
  TotalNold[x]<-sum(tmp)
}

TotalNnew[which(LengthBins==lowerBnd)]<-sum(TotalNnew[1:(which(LengthBins==lowerBnd))])
TotalNold[which(LengthBins==lowerBnd)]<-sum(TotalNold[1:(which(LengthBins==lowerBnd))])

TotalNnew[which(LengthBins==upperBnd)]<-sum(TotalNnew[(which(LengthBins==upperBnd)):length(LengthBins)])
TotalNold[which(LengthBins==upperBnd)]<-sum(TotalNold[(which(LengthBins==upperBnd)):length(LengthBins)])
TotalNnew[is.na(TotalNnew)]<-0
TotalNold[is.na(TotalNold)]<-0

#==take the observed retained numbers, multiply to get the numbers at len by shell cond
#==observed retained numbers come from "Item 3 and 5...fish ticket snow crab.csv"
write.table(rbind(TotalNnew[which(LengthBins==lowerBnd):which(LengthBins==upperBnd)],TotalNold[which(LengthBins==lowerBnd):which(LengthBins==upperBnd)]),
 "retained_catch_length_comps_new_then_old.txt",
 row.names=FALSE,col.names=F)
par(mfrow=c(2,1))
barplot(TotalNnew)
barplot(TotalNold)

#================================================
# Total length frequencies
# modified from "item 2_20XX-20XX discard length frequency.xlx"
#================================================
all_data<-read.csv("data/item2_snow_crab_observer_size_comp.csv")

#=this number comes from "Item 4..."  number of observed Males caught in EBS snow fishery
total_crab_num<-77897274/1000
ret_crab_num<-28626114/1000
disc_num<-total_crab_num-ret_crab_num

disc_num/ret_crab_num


#==these are not the same, but they need to be for calculating discard length comps
sum(TotalNold,TotalNnew)
ret_crab_num

77897274-28626114
ret_crab_num/total_crab_num
ret_crab_num-total_crab_num

# #==males old and new shell total
#==from directed fishery
library(dplyr)
data<-filter(all_data,management_area=="bering_sea"&sex=="male")

LengthBins			<-seq(25,165,5)
tot_len_comp_old		<-rep(0,length(LengthBins))
 for(y in 1:(length(LengthBins)-1))
 {
   tot_len_comp_old[y]	<-sum(data$old[data$size>=LengthBins[y] & data$size<LengthBins[y+1]],
                         data$old[data$size>=LengthBins[y] & data$size<LengthBins[y+1]],na.rm=T)
 }

tot_len_comp_new		<-rep(0,length(LengthBins))
for(y in 1:(length(LengthBins)-1))
{
  tot_len_comp_new[y]	<-sum(data$new[data$size>=LengthBins[y] & data$size<LengthBins[y+1]],
                            data$new[data$size>=LengthBins[y] & data$size<LengthBins[y+1]],
                            data$new[data$size>=LengthBins[y] & data$size<LengthBins[y+1]],na.rm=T)
}


#==the raw numbers for both the observer and dockside length comps are not 'right', so we have to make sense of this
#==this requires getting them in a common currency
#==this will also be not necessary with GMACS

#==compositionify totals
tot_temp<-sum(tot_len_comp_new,tot_len_comp_old)
tot_comp_new<-tot_len_comp_new/tot_temp
tot_comp_old<-tot_len_comp_old/tot_temp

#==compositionify retained
ret_tot<-sum(TotalNnew,TotalNold)
ret_comp_new<-TotalNnew/ret_tot
ret_comp_old<-TotalNold/ret_tot

#==scale to the tot leng comps to total number of crab caught 
use_tot_new_n_len<-total_crab_num*tot_comp_new
use_tot_old_n_len<-total_crab_num*tot_comp_old

#==scale ret len comps to retained num caught
use_ret_new_n_len<-ret_crab_num*ret_comp_new
use_ret_old_n_len<-ret_crab_num*ret_comp_old

#==subtract retained from total
use_disc_new_len<-use_tot_new_n_len-use_ret_new_n_len
use_disc_old_len<-use_tot_old_n_len-use_ret_old_n_len

#==D'oh...there are negatives...'fix' them
use_disc_new_len[use_disc_new_len<0]<-0
use_disc_old_len[use_disc_old_len<0]<-0

plot(use_disc_new_len~LengthBins)
par(mfrow=c(3,1),mar=c(1,3,1,1))
barplot(use_tot_new_n_len,ylim=c(0,16000))
barplot(use_ret_new_n_len,ylim=c(0,16000))
barplot(use_disc_new_len,ylim=c(0,16000))

par(mfrow=c(4,1),mar=c(1,3,1,1))
barplot(use_tot_old_n_len,ylim=c(-1000,1600))
barplot(use_ret_old_n_len,ylim=c(-1000,1600))
barplot(use_disc_old_len,ylim=c(-1000,1600))
barplot(use_disc_old_len+use_ret_old_n_len,ylim=c(-1000,1600))

write.table(round(rbind(use_disc_new_len[1:22],use_disc_old_len[1:22]),1),"male_disc_comps_new_then_old.txt",row.names=FALSE,col.names=F)

# #==males old and new shell total
#==from directed fishery
library(dplyr)
data<-filter(all_data,management_area=="bering_sea"&sex=="female")
total_disc_num_f<-239332/1000 # from item 4
LengthBins			<-seq(25,165,5)
tot_len_comp_fem		<-rep(0,length(LengthBins))
for(y in 1:(length(LengthBins)-1))
{
  tot_len_comp_fem[y]	<-sum(data$total[data$size>=LengthBins[y] & data$size<LengthBins[y+1]],
                            data$total[data$size>=LengthBins[y] & data$size<LengthBins[y+1]],na.rm=T)
}

tot_f<-sum(tot_len_comp_fem)
tot_comp_f<-tot_len_comp_fem/tot_f
use_disc_n_len_f<-tot_comp_f*total_disc_num_f
plot(use_disc_n_len_f~LengthBins)
write.table(c(use_disc_n_len_f[1:22]),"female_disc_comps.txt",row.names=FALSE,col.names=F)


#=======================================================================
# bycatch length composition
# these data come from the "Observer data" tab on AKfin
# Enter 'NORPAC Length Report - Haul & Length"
# Download the data for snow crab
# Time period should be July 1 (last year) to June 30 (this year)
#=======================================================================

LenDatBig<-read.csv("data/norpac_length_report.csv")
names(LenDatBig) # check that you chop the top off (or you could skip over the lines in read.csv...)

LenDatBig$Haul.Offload.Date<-strptime(LenDatBig$Haul.Offload.Date,format="%d-%b-%y")
range(LenDatBig$Haul.Offload.Date)

#==CHECK DATES, EH!
LenDat<-LenDatBig[LenDatBig$Haul.Offload.Date >= "2019-07-01" & LenDatBig$Haul.Offload.Date <= "2020-06-30" & LenDatBig$Species.Name=="OPILIO TANNER CRAB",]
nrow(LenDat)
nrow(LenDatBig)

#==males old shell discard
LengthBins			<-seq(25,165,5)
BycatchFem			<-rep(0,length(LengthBins))
BycatchMale			<-rep(0,length(LengthBins))

 for(y in 1:(length(LengthBins)-1))
 {
  BycatchFem[y]	<-sum(LenDat$Frequency[LenDat$Length..cm.>=LengthBins[y] & LenDat$Length..cm.<LengthBins[y+1] & LenDat$Sex=="F"])
  BycatchMale[y]	<-sum(LenDat$Frequency[LenDat$Length..cm.>=LengthBins[y] & LenDat$Length..cm.<LengthBins[y+1] & LenDat$Sex=="M"])
 }

upperBnd							<-130
BycatchFem[which(LengthBins==upperBnd)]	<-sum(BycatchFem[(which(LengthBins==upperBnd)):length(LengthBins)])
BycatchFem						<-BycatchFem[1:which(LengthBins==upperBnd)]
BycatchMale[which(LengthBins==upperBnd)]	<-sum(BycatchMale[(which(LengthBins==upperBnd)):length(LengthBins)])
BycatchMale						<-BycatchMale[1:which(LengthBins==upperBnd)]

par(mfrow=c(1,2))
barplot(BycatchFem)
barplot(BycatchMale)

write.table(rbind(BycatchFem,BycatchMale),"bycatch_len_comps_f_then_m.txt",
 row.names=FALSE,col.names=F)

#=======================================================================
# bycatch numbers
# these data come from the "Observer data" tab on AKfin
# Enter 'NORPAC catch Report"
# Download the data for snow crab in BS of BSAI
# Time period should be July 1 (last year) to June 30 (this year)
#=======================================================================


#==WHAT IS THE POINT OF THIS CHUNK OF CODE?
bycatch_dat<-read.csv("data/norpac_catch_report.csv",skip=6)
temp<-strptime(bycatch_dat$Haul.Date,format="%d-%b-%y")
bycatch_dat$Haul.Date<-substr(temp,start=1,stop=10)

#==bycatch numbers total
  bycatchDat<-bycatch_dat[bycatch_dat$Haul.Date >= "2019-07-01" & bycatch_dat$Haul.Date <= "2020-06-30"& 
                            bycatch_dat$Species.Name=="OPILIO TANNER CRAB" & bycatch_dat$Gear.Description!="POT OR TRAP" ,]
  bycatch_num_tot<-sum(bycatchDat$Extrapolated.Number,na.rm=T)
  bycatch_wt_tot<-sum(bycatchDat$Extrapolated.Weight..kg.,na.rm=T)
  
  #==bycatch numbers total
  bycatchDat<-bycatch_dat[bycatch_dat$Haul.Date >= "2018-07-01" & bycatch_dat$Haul.Date <= "2019-06-30"& 
                            bycatch_dat$Species.Name=="OPILIO TANNER CRAB" & bycatch_dat$Gear.Description!="POT OR TRAP" ,]
  bycatch_num_tot<-sum(bycatchDat$Extrapolated.Number,na.rm=T)
  
  bycatch_wt_tot<-sum(bycatchDat$Extrapolated.Weight..kg.,na.rm=T)
#=========================================
# all years bycatch
#==================================

bycatch_dat_big<-read.csv("data/norpac_catch_report_big.csv",skip=6)
names(bycatch_dat_big)
temp<-strptime(bycatch_dat_big$Haul.Date,format="%d-%b-%y")
bycatch_dat_big$Haul.Date<-substr(temp,start=1,stop=10)
bycatch_dat_big$crab.year<-bycatch_dat_big$Year

unique(bycatch_dat_big$Gear.Description)

#==bycatch numbers total
bycatch_year    <-sort(unique(bycatch_dat_big$Year))
bycatch_num_tot		<-rep(0,length(bycatch_year))
bycatch_wt_tot		<-rep(0,length(bycatch_year))

for(y in 1:(length(bycatch_year)-1))
{  
  bycatchDat<-bycatch_dat_big[bycatch_dat_big$Haul.Date >= paste(bycatch_year[y],"-07-01",sep="") & bycatch_dat_big$Haul.Date <= paste(bycatch_year[y]+1,"-06-30",sep="") & 
                                bycatch_dat_big$Species.Name=="OPILIO TANNER CRAB" & bycatch_dat_big$Gear.Description!="POT OR TRAP" ,]
  bycatch_num_tot[y]<-sum(bycatchDat$Extrapolated.Number,na.rm=T)
  bycatch_wt_tot[y]<-sum(bycatchDat$Extrapolated.Weight..kg.,na.rm=T)  
}
write.table(cbind(bycatch_year,bycatch_num_tot),"bycatch_numbers_total.txt")
cbind(bycatch_year,bycatch_wt_tot/1000000)

plot(bycatch_wt_tot/1000000~bycatch_year,type='l')
write.table(cbind(bycatch_year,bycatch_wt_tot/1000000),"bycatch_wt_total.txt")

#==bycatch by gear type
library(dplyr)
library(ggplot2)
for(y in 1:(length(bycatch_year)-1))
  bycatch_dat_big$crab.year[bycatch_dat_big$Haul.Date >= paste(bycatch_year[y],"-07-01",sep="") & bycatch_dat_big$Haul.Date <= paste(bycatch_year[y]+1,"-06-30",sep="")]<-bycatch_year[y]

in_dat<-bycatch_dat_big[,-24]
temp<- in_dat %>%
  group_by(crab.year,Gear.Description) %>%
  summarise(Bycatch = sum(Extrapolated.Number))

temp$Year<-as.numeric(temp$Year)

p<-qplot(x=crab.year,y=Bycatch,col=Gear.Description,data=temp)
png("plots/bycatch.png",height=8,width=8,res=400,units='in')
p + geom_line() + .THEME
dev.off()



#=======================================================================
# Plot the proportion of bycatch in each fishery
#======================================================================
big_bycatch_09_17<-read.csv("2017_snowcrab_catch_data/Crab Bycatch Estimates_2009_2017.csv")
big_bycatch_91_09<-read.csv("2017_snowcrab_catch_data/Crab Bycatch Estimates_1991_2009.csv")

names(big_bycatch_91_09)

library(plyr)
library(ggplot2)
library(dplyr)

#==THIS IS THE ONE THAT GOES IN THE DAT FILE
tot_09_17<- big_bycatch_09_17 %>%
  group_by(Crab.Year) %>%
  summarise(Bycatch = sum(Estimate.Num))

temp<- big_bycatch_09_17 %>%
  group_by(Crab.Year,Agency.Gear.Code) %>%
  summarise(Bycatch = sum(Estimate.Num))
names(temp)[1]<-"Year"
names(temp)[2]<-"Gear"
temp$Year<-as.numeric(temp$Year)

temp2<- big_bycatch_91_09 %>%
  group_by(ï..Crab.Year,Gear) %>%
  summarise(Bycatch = sum(Estimate.Num..Sum.))


names(temp2)[1]<-"Year"
temp2$Year<-as.numeric(substr(temp2$Year,1,4))

stupid<-data.frame(rbind(as.matrix(temp),as.matrix(temp2)))
stupid$Bycatch<-as.numeric(as.character(stupid$Bycatch))
stupid$Year<-as.numeric(as.character(stupid$Year))
p<-qplot(x=Year,y=Bycatch,col=Gear,data=stupid)
p + geom_line()

p<-qplot(x=Year,y=Bycatch,col=Gear,data=stupid[stupid$Year>1992,])
p + geom_line()



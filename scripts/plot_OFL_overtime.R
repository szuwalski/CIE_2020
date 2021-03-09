library(ggplot2)
library(dplyr)
dat<-data.frame(read.csv("management_record.csv"))

#==plot theme
.THEME    = theme_bw(base_size = 12, base_family = "") +
  theme(strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(color="white",fill="white")) 
.COL =  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

taker<-c("ABC","OFL","total GHL","total TAC")
indat<-filter(dat,Quantity %in% taker)

both_proj<-ggplot() +
  geom_line(data=indat,aes(x=Year,y=Biomass,col=Quantity)) +
  geom_point(data=indat,aes(x=Year,y=Biomass,col=Quantity)) +
.THEME +
  theme(legend.position = c(.15,.85),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14)) + 
  expand_limits(y=0)
print(both_proj)



png("plots/projection_20.png",height=5,width=8,res=400,units='in')
print(both_proj)
dev.off()


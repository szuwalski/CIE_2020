
setwd("C:/snow_2020/retro/")

drop_dir<-list.dirs(recursive=FALSE)

for(x in 1:length(drop_dir))
{
 file.remove(paste(drop_dir[x],"/gmacs.exe",sep=""))
 file.copy("C:/gmacs_v2/src_rbar/gmacs.exe", drop_dir[x])
 
}



for(x in 1:length(drop_dir))
{
 file.copy("C:/snow_2020/19_gmacs/snow.CTL", drop_dir[x],overwrite=TRUE)
}
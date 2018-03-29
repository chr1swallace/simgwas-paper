#!/usr/bin/env Rscript
files <- list.files("timings",full=TRUE)
files
library(data.table)
library(ggplot2)
library(cowplot)
library(randomFunctions)
theme_set(theme_cowplot())


reader <- function(f) {
x <- fread(f,fill=TRUE,sep=" ",header=FALSE)
setnames(x, c("id","name","state","total","user"))
x[,min:=as.numeric(sub(":.*","",user))]
x[,sec:=as.numeric(sub(".*:","",user)) + 60 * min]
x[,method:=sub("-.*","",name)]
x[,method:=sub("hg","HapGen2 + SNPTEST2",method)]
x[,method:=sub("meth","simGWAS",method)]
x[,item:=as.numeric(sub(".*-","",name))]
x <- x[,.(n=sum(!is.na(sec)),mean=mean(sec,na.rm=TRUE),sd=sd(sec,na.rm=TRUE)),by=c("method","item")]
copy(x)
}

x <- reader(files[1])
p <- ggplot(x,aes(x=item,y=mean,ymin=mean-sd,ymax=mean+sd,col=method,group=method)) + geom_pointrange() + geom_smooth(method="lm") + labs(x="Number of causal variants",y="Time (seconds)")
nicecow(p)


x <- reader(files[2])
x[,mean:=mean/60]
x[,sd:=sd/60]
p <- ggplot(x,aes(x=item,y=mean,ymin=mean-sd,ymax=mean+sd,col=method,group=method)) + geom_pointrange() + geom_smooth(method="lm") + scale_x_log10("Number of simulations (log scale)") + scale_y_log10("Time (minutes, log scale)")
nicecow(p)

x <- reader(files[3])
x[,mean:=mean/60]
x[,sd:=sd/60]
p <- ggplot(x,aes(x=item,y=mean,ymin=mean-sd,ymax=mean+sd,col=method,group=method)) + geom_pointrange() + geom_smooth(method="lm") + scale_x_log10("Sample size (log scale)") + scale_y_log10("Time (minutes, log scale)")
nicecow(p)



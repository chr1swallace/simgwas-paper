#!/usr/bin/env Rscript
files <- list.files("timings",full=TRUE)
files
library(data.table)
library(ggplot2)
library(cowplot)
library(randomFunctions)
theme_set(theme_cowplot())

RESULTS <- vector("list",3)
reader <- function(f) {
x <- fread(f,fill=TRUE,sep=" ",header=FALSE)
setnames(x, c("id","name","state","total","user"))
x[,min:=as.numeric(sub(":.*","",user))]
x[,sec:=as.numeric(sub(".*:","",user)) + 60 * min]
x[,method:=sub("-.*","",name)]
x[,method:=sub("hg","HAPGEN+SNPTEST",method)]
x[,method:=sub("meth","simGWAS",method)]
x[,item:=as.numeric(sub(".*-","",name))]
x <- x[,.(nsim=sum(!is.na(sec)),mean=mean(sec,na.rm=TRUE),sd=sd(sec,na.rm=TRUE)),by=c("method","item")]
copy(x)
}

x <- reader(files[1])
x$item.what="n.causal.vars"
RESULTS[[1]] <- x
p1 <- ggplot(x,aes(x=item,y=mean,ymin=mean-sd,ymax=mean+sd,col=method,group=method)) + geom_pointrange() +  labs(x="Number of causal variants",y="Time (seconds)")

x <- reader(files[2])
x$item.what="n.replications"
RESULTS[[2]] <- x
x[,mean:=mean/60]
x[,sd:=sd/60]
p2 <- ggplot(x,aes(x=item,y=mean,ymin=mean-sd,ymax=mean+sd,col=method,group=method)) + geom_pointrange()  + scale_x_log10("Number of replications (log scale)",breaks=c(1,4,16,64,128)) + scale_y_log10("Time (minutes, log scale)")

x <- reader(files[3])
x$item.what="n.cases.controls"
RESULTS[[3]] <- x
x[,mean:=mean/60]
x[,sd:=sd/60]
p3 <- ggplot(x,aes(x=item,y=mean,ymin=mean-sd,ymax=mean+sd,col=method,group=method)) + geom_pointrange() + scale_x_log10("Sample size (log scale)",breaks=c(1000, 2000, 4000, 8000, 16000, 32000, 64000)) + scale_y_log10("Time (minutes, log scale)")

cols <- rgb(red=c(0,0,182,146,146,219)/256, green=c(146,109,109,0,73,209)/256, blue=c(146,219,255,0,0,0)/256)
nicer <- function(p) {
    nicecow(p) + geom_smooth(method="lm",linetype="dashed",se=FALSE) + theme(legend.position="none") +
      scale_colour_manual(values=c("HAPGEN+SNPTEST"=cols[2],"simGWAS"=cols[5]))
}
plot_grid(nicer(p1),
          nicer(p2) + theme(legend.position=c(0,0.9)),
          nicer(p3),
          labels=c("a","b","c"),nrow=3)

ggsave("timings.pdf",height=9,width=6)

write.table(do.call("rbind",RESULTS),file="timings.tsv",row.names=FALSE,sep="\t",quote=FALSE)

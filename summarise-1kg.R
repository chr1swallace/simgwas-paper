library(ggplot2)
library(cowplot)
library(ggridges)
library(magrittr)
library(data.table)
library(randomFunctions)
source("~/DIRS.txt") # sets SIMGWAS, REFDATA, BINDIR
source("~/Projects/simgwas/getspec.R")
d <- SIMGWAS
cols <- as.list(rgb(red=c(0,0,182,146,146,219)/256, green=c(146,109,109,0,73,209)/256, blue=c(146,219,255,0,0,0)/256,alpha=0.4))
names(cols)[c(1,2,5)] <- c("sG","HG","F")
cols <- list(sG="#Ca3542",HG="#27647B",F="#57575F")
cols <- list(sG="#F68930",HG="#697F90",F="#DA5526")
## read in haplotypes
dref <- file.path(REFDATA,"1000GP_Phase3")
fhg <- file.path(d,"input")#tempfile(tmpdir=dtmp)
h <- fread(paste0(fhg,".hap"))
h <- as.matrix(h)
snps <- fread(paste0(fhg,".leg"))

devtools::load_all("~/RP/simGWAS")
## library(simGWAS)
library(mvtnorm)
library(corpcor)
h1 <- h[, seq(1,ncol(h)-1,by=2)] #[,samples$GROUP=="AFR"]
h2 <- h[, seq(2,ncol(h),by=2)] #[,samples$GROUP=="AFR"]
freq <- as.data.frame(t(cbind(h1,h2))+1)
snps$rs <- make.names(snps$id)
colnames(freq) <- snps$rs
use <- snps$AFR>0.01 & snps$AFR < 0.99 & apply(freq,2,var)>0
freq <- freq[,use,drop=FALSE]
dfsnps <- snps[use,,drop=FALSE]
dfsnps$maf <- colMeans(freq-1)
LD <- cor(freq)
LD <- as.matrix(make.positive.definite(LD))
freq$Probability <- 1/nrow(freq)

getproxy <- function(x,lower,upper) {
    res <- apply(abs(x),2,function(.) which(.>lower & .<upper))
    l <- sapply(res,length)
    if(any(l>0))
        res[l > 0] <- lapply(res[ l > 0], sample, 1)
    if(any(l==0))
        res[l==0] <- lapply(res[l==0], function(.) NA)
    setdiff(unlist(res),NA)
}
get0 <- function(x,upper) {
    mx <- apply(abs(x),1,max)
    res <- which(mx<upper)
    if(length(res)>=ncol(x)) 
        return(sample(res,ncol(x)))
    setdiff(c(res,rep(NA,ncol(x)-length(res))),NA)
}

loadset <- function(f) {
    RESULTS <- vector("list",length(f))
    for(i in seq_along(f)) {
        (load(f[[i]]))
        results[,gsim:=paste(i,sim,sep="-")]
        RESULTS[[i]] <- results[,.(snp,sim,beta1,v,beta.hg,v.hg,gsim)]
    }
    rbindlist(RESULTS)
}

## read in results
## big number of sims on a small set of CVs
files <- list.files(d,pattern="sim",full=TRUE)
## split files by category
fcat <- basename(files)  %>%  sub("-[0-9a-f]+.RData","",.)
files <- split(files,fcat)
## to keep results figures readable, only use 1:4,6
#files <- files[ grep("spec5|spec7|spec8",names(files),invert=TRUE) ]
message("files found: ",length(unlist(files)))
pb <- txtProgressBar(min=1,max=length(files),style=3)
RESULTS <- vector("list",length(files))
for(i in seq_along(files)) { # i=33
    setTxtProgressBar(pb, i)
    results <- loadset(files[[i]])
    SPECIAL <- sub(".*spec","",names(files)[i])
    args <- list(SPECIAL=SPECIAL)
    source("~/Projects/simgwas/getspec.R")
    results[,z1:=beta1/sqrt(v)]
    ## results[,z2:=beta2/v]
    results[,z.hg:=beta.hg/sqrt(v.hg)]
    ## proxy <- getproxy(LD[,CV,drop=FALSE],sqrt(0.1), sqrt(0.4))
    unlinked <- get0(LD[,CV,drop=FALSE],0.4)
    ## asnps <- c(CV,setdiff(results$snp,CV)) #dfsnps[c(setdiff(proxy,NA),setdiff(unlinked,NA)),]$rs)
    asnps <- c(CV,dfsnps[c(#proxy,
               unlinked[1]),]$rs)
    lg1 <- c(log(g1),rep(0,length(asnps)-length(g1)))
    if(length(asnps)>1) {
        ld <- LD[asnps,asnps]
        beta.marg <- ld %*% lg1
    } else {
        beta.marg <- lg1
    }
    tmp2 <- data.table(snp=asnps,class=c(rep("causal",length(CV)),
                                         #rep("proxy",length(proxy))),
                                         rep("unlinked",1)), #length(unlinked))),
                       beta.marg=beta.marg[,1],
                       cat=names(files)[i])
    tmp2[,snplab:=paste(cat,seq_along(asnps),class,sep="-")]
    tmp2 <- tmp2[!duplicated(snp),]
    ## res <- results[snp %in% asnps , .(beta1=sort(beta1),
    ##                                   ## beta2=sort(beta2),
    ##                                   beta.hg=sort(beta.hg),
    ##                                   z1=sort(z1),
    ##                                   ## z2=sort(z2),
    ##                                   z.hg=sort(z.hg),
    ##                                   az1=sort(abs(z1)),
    ##                                   ## az2=sort(abs(z2)),
    ##                                   az.hg=sort(abs(z.hg))),
    ##                by=c("snp")]
    res <- results[snp %in% asnps , .(beta1=sort(beta1),
                                      ## beta2=sort(beta2),
                                      beta.hg=sort(beta.hg),
                                      z1=sort(z1),
                                      ## z2=sort(z2),
                                      z.hg=sort(z.hg),
                                      az1=sort(abs(z1)),
                                      ## az2=sort(abs(z2)),
                                      az.hg=sort(abs(z.hg))),
                   by=c("snp")]
    RESULTS[[i]] <- merge(results,tmp2,by="snp")
}

## forward sims
ffiles <- list.files(d,pattern="forward",full=TRUE)

################################################################################

## FIGURE 1
## manhattans

## FIRST 5000 snps in Africans only
dref <- file.path(REFDATA,"1000GP_Phase3")
fhg <- file.path(d,"input")#tempfile(tmpdir=dtmp)
library(data.table)
snps <- fread(paste0(fhg,".leg"))
snps[,snp:=make.names(id)]
makeman <- function(i) {
    stopifnot(i<=5 & i>=1)
    l <- getspec(as.character(i))
    man <- loadset(files[[5+i]])
    man[,z.sg:=-log10(pnorm(abs(beta1/sqrt(v)),lower.tail=FALSE)*2)]
    man[,z.hg:=-log10(pnorm(abs(beta.hg/sqrt(v.hg)),lower.tail=FALSE)*2)]
    summary(man)
    man <- man[!is.na(beta.hg),]
    man2 <- man[,.(lo.sg=quantile(z.sg,0.25),med.sg=median(z.sg),hi.sg=quantile(z.sg,0.75),
                   lo.hg=quantile(z.hg,0.25),med.hg=median(z.hg),hi.hg=quantile(z.hg,0.75)),
                by="snp"]
    man2[ ,CV:=snp %in% l$CV]
    merge(man2,snps[,.(snp,position)],by="snp")
}

plotman <- function(i=0,man=NULL,posl=536,posr=705) {
    if(i==0) {
        man2 <- man
    } else {
        man2 <- makeman(i)
    }
    man3 <- melt(man2,c("snp","position","CV"),list(lo=c("lo.sg","lo.hg"),
                                               med=c("med.sg","med.hg"),
                                               hi=c("hi.sg","hi.hg")))
    man3[,variable:=c("sG","HG")[variable]]
    man3 <- man3[order(position),]
    man3[,mx:=min(lo),by="position"]
    ## man3 <- man3[order(mx),]
    man3[,position:=as.numeric(position)]
    man3[,position:=as.numeric(1:.N),by="variable"]
    ## man3 <- man3[mx>1,]
    pmain <- ggplot(man3,aes(x=position,ymin=lo,ymax=hi,y=abs(med),pch=variable)) +
      geom_vline(aes(xintercept=position),data=man3[CV==TRUE,],linetype="dotted") +
      geom_linerange() +
      geom_point() +
      labs(x="SNP",y="-log10(p)") +
      ## scale_color_manual(values=c(sG=cols$F,HG=cols$HG)) +
      ## scale_linetype_manual("method",values=c(sG="solid",HG="dashed")) + 
      ## scale_shape_discrete("method") + 
      ## scale_color_manual("method",values=c(sG="gray",HG="black")) +
      theme(legend.position="none") + #c(0.1,0.9),legend.background = element_rect(fill="white",linetype=1,colour="black",size=0.5)) +
      facet_grid(variable ~ .)
    pmain
    
    ## man4 <- man3[CV==TRUE,]
    ## man4[,position:=as.numeric(factor(position))]
    ## man4[variable=="sG",position:=position+0.3]
    ## pl <- ggplot(man4,aes(x=position,ymin=lo,ymax=hi,y=abs(med),col=variable,pch=variable)) +
    ##   geom_vline(aes(xintercept=position),data=man4,linetype="dotted") +
    ##   geom_linerange() +
    ##   geom_point() +
    ##   labs(x="SNP",y="-log10(p)") +
    ##   scale_color_manual("method",values=c(sG=cols$F,HG=cols$HG)) +
    ##   ## scale_linetype_manual("method",values=c(sG="solid",HG="dashed")) + 
    ##   scale_shape_discrete("method") + 
    ##   ## scale_color_manual("method",values=c(sG="gray",HG="black")) +
    ##   theme(legend.position="none")
    ## pl
    
}
plotman(man=mans[[4]]) + xlim(450,750) + background_grid()

mans <- lapply(1:5,makeman)

w <-0.8 
pdf(file="manhattan.pdf",height=6*w,width=10*w)
plotman(man=mans[[4]]) + xlim(450,720) + background_grid()
dev.off()

## ggplot(man,aes(x=med.sg,y=med.hg)) + geom_point() + geom_abline()
## ggplot(man,aes(x=hi.sg-lo.sg,y=hi.sg-lo.hg)) + geom_point() + geom_abline()


##   man3 <- melt(man,c("snp","position","CV"),list(lo=c("lo.sg","lo.hg"),
##                                                med=c("med.sg","med.hg"),
##                                                hi=c("hi.sg","hi.hg")))
##     man3[,variable:=c("sG","HG")[variable]]
##     man3 <- man3[order(position),]
##     man3[,mx:=min(lo),by="position"]
##     ## man3 <- man3[order(mx),]
##     man3[,position:=as.numeric(position)]
## man3[,position:=as.numeric(1:.N),by="variable"]

##      ggplot(man,aes(x=med.sg,y=hi.sg-lo.sg,y=hi.sg-lo.hg)) + geom_point() + geom_abline()
## man3 <- man3[order(med),]
## ggplot(man3,aes(x=med,y=hi-lo,col=variable)) + geom_point() + geom_path()


################################################################################

## FIGURE 2
## qq plots

results <- rbindlist(RESULTS)
## ## relabel spec6 as spec5
## results[,cat:=sub("spec6","spec5",cat)]
table(results$cat)
results[,lab:=gsub("sim-c.*-spec|-causal|-unlinked|-proxy","",snplab)]
results[,n:=gsub("sim-c|-spec.*","",cat)]
table(paste(results$cat,results$lab),results$class)
## results[,lab:=sub("6-","5-",lab)]
results <- results[!is.na(beta.hg),] # one failure
results[,i:=sample(1:.N),by=c("lab","n")]
results <- results[i <= 1000,] # fix n sims = 1000 per scenario
table(paste(results$cat,results$lab),results$class)

## add forward sims
ffiles <- list.files(d,pattern="forward",full=TRUE)
fres <- lapply(ffiles, function(f) {
    eval(as.symbol(load(f)))
})  %>%  rbindlist()
table(fres$cat,fres$snp)

fres[,i:=1:.N,by=c("cat","snp")]
results <- merge(results,fres,by=c("cat","snp","i"),all.x=TRUE)

## table of summary stats
BY <- c("class","lab","n")

summ <- copy(results)
setnames(summ,c("beta1","z1"),c("beta.sg","z.sg"))
summ <- melt(summ,BY,list(beta=c("beta.sg","beta.hg","beta.f"),
                          z=c("z.sg","z.hg","z.f")))
summ[,method:=c("sG","HG","F")[variable]]
summ2 <- summ[,.(mn.b=mean(beta),
                mn.z=mean(z),
                sd.b=sd(beta),
                sd.z=sd(z)),by=c(BY,"method")]
summ2 <- melt(summ2,c(BY,"method"),list(mn=c("mn.b","mn.z"),sd=c("sd.b","sd.z")))
summ2[,stat:=c("beta","Z")[variable]]
summ2[,method:=factor(method,levels=c("sG","F","HG"))]
summ2[,x:=as.numeric(factor(lab)) + (as.numeric(as.factor(method))-2)/5]
lev <- levels(factor(summ2$lab))
br <- seq_along(lev)
summ2 <- summ2[!is.na(mn),]

ggplot(summ2[stat=="beta",],
       aes(x=x,y=mn,ymin=mn-sd,ymax=mn+sd,col=method)) +
  geom_pointrange() +
  geom_hline(yintercept=0,linetype="dotted") +
  scale_x_continuous("Scenario-SNP",breaks=br,labels=lev) +
  facet_wrap(~n) +
  background_grid()
  
results[class!="causal",beta.f:=0]
results[class=="causal",.(sg=mean(beta1),hg=mean(beta.hg),f=mean(beta.f)),by=c("class","lab","n")]
results[class=="causal",.(sg=sd(beta1),hg=sd(beta.hg),f=sd(beta.f)),by=c("class","lab","n")]
results[is.na(beta.f),beta.f:=0]
results[is.na(v.f),v.f:=0]
results[is.na(z.f),z.f:=0]
mean.tests.z <- results[, .(sg.hg=t.test(z1,z.hg)$p.value,
                               sg.f=t.test(z1,z.f)$p.value,
                            f.hg=t.test(z.hg,z.f)$p.value,
                            test="mean",
                            stat="Z"),
                           by=BY]
ks.tests.z <- results[, .(sg.hg=ks.test(z1,z.hg)$p.value,
                             sg.f=ks.test(z1,z.f)$p.value,
                             f.hg=ks.test(z.hg,z.f)$p.value,
                          test="KS",
                          stat="Z"),
                         by=BY]
mean.tests.beta <- results[, .(sg.hg=t.test(beta1,beta.hg)$p.value,
                               sg.f=t.test(beta1,beta.f)$p.value,
                               f.hg=t.test(beta.hg,beta.f)$p.value,
                               test="mean",
                               stat="beta"),
                           by=BY]
ks.tests.beta <- results[, .(sg.hg=ks.test(beta1,beta.hg)$p.value,
                             sg.f=ks.test(beta1,beta.f)$p.value,
                             f.hg=ks.test(beta.hg,beta.f)$p.value,
                             test="KS",
                             stat="beta"),
                         by=BY]
tests <- rbind(mean.tests.z,ks.tests.z,mean.tests.beta,ks.tests.beta)
tests[class!="causal",f.hg:=NA]
tests[class!="causal",sg.f:=NA]

stats.beta <- results[,.(sg.mean=mean(beta1), hg.mean=mean(beta.hg), f.mean=mean(beta.f),
                         sg.sd=sd(beta1),hg.sd=sd(beta.hg), f.sd=sd(beta.f),
                         stat="beta"),by=BY]
stats.z <- results[,.(sg.mean=mean(z1),hg.mean=mean(z.hg),f.mean=mean(z.f),
                      sg.sd=sd(z1),hg.sd=sd(z.hg), f.sd=sd(z.f),
                      stat="Z"),by=BY]
stats <- rbind(stats.beta,stats.z)
stats[class!="causal",f.mean:=NA]
stats[class!="causal",f.sd:=NA]

m <- merge(tests[test=="mean",],tests[test=="KS",],by=c(BY,"stat"),suffixes=c(".testmean",".testKS"))
m[,test.testmean:=NULL]
m[,test.testKS:=NULL]
m <- merge(stats,m,by=c(BY,"stat"))
head(m)


setnames(m,sub("f.","Forward.",names(m)))
setnames(m,sub("hg.","HapGen.",names(m)))
setnames(m,sub("sg.","simGWAS.",names(m)))
m[class=="unlinked",lab:=sub("-.*","-0",lab)]
head(m)
head(m[class=="unlinked",])
m <- m[order(lab,n,stat),]
m
fwrite(m,file="~/simgwas-supptab.csv")

## library(Hmisc)



## p.mean.beta <- results[, .(p.mean=t.test(beta1,beta.hg)$p.value),by=BY]
## p.mean.z <- results[, .(p.mean=t.test(z1,z.hg)$p.value),by=BY]
## ks.beta <- results[!is.na(beta.hg),ks.test(beta1,beta.hg)[c("statistic","p.value")],by=BY]
## ks.z <- results[!is.na(beta.hg),ks.test(z1,z.hg)[c("statistic","p.value")],by=BY]
## mr.beta <- results[!is.na(beta.hg),.(md.sg=median(beta1),md.hg=median(beta.hg),iqr.sg=quantile(beta1,.75)-quantile(beta1,.25),iqr.hg=quantile(beta.hg,0.75)-quantile(beta.hg,0.25)),by=BY]
## mr.z <- results[!is.na(z.hg),.(md.sg=median(z1),md.hg=median(z.hg),iqr.sg=quantile(z1,.75)-quantile(z1,.25),iqr.hg=quantile(z.hg,0.75)-quantile(z.hg,0.25)),by=BY]
## n <- results[!is.na(beta.hg),.(nsim=.N),by=BY]
## m1 <- merge(ms.beta,ks.beta,by=BY)  %>%  merge(.,p.mean.beta)
## m2 <- merge(ms.z,ks.z,by=BY)  %>% merge(.,p.mean.z)
## m1$n  %<>% as.numeric()
## m2$n  %<>% as.numeric()
## m1[class=="unlinked",lab:=paste0(sub("-.","",lab),"-0")]
## m2[class=="unlinked",lab:=paste0(sub("-.","",lab),"-0")]
## m1 <- m1[order(m1$lab),]
## m2 <- m2[order(m2$lab),]

## m1[class=="causal" & n==5000,.(lab,n,mn.sg,mn.hg,sd.sg,sd.hg)][lab %in% c("5-1","5-2"),]
## m2[class=="causal" & n==5000,.(lab,n,mn.sg,mn.hg,sd.sg,sd.hg)][lab %in% c("5-1","5-2"),]

## latex(m1[,.(lab,n,mn.sg,mn.hg,p.mean,sd.sg,sd.hg,p.value)],
##       digits=2,rowname=NULL,
##       cgroup=c("Scenario","n","Mean","Std. Dev.","K-S"),
##       n.cgroup=c(1,1,3,2,1),
##       colheads=c("","","sG","HG","p","sG","HG","p"),
##       file="beta.tex",
##       table.env=FALSE)

## latex(m2[,.(lab,n,mn.sg,mn.hg,p.mean,sd.sg,sd.hg,p.value)],
##       digits=2,rowname=NULL,
##       cgroup=c("Scenario","n","Mean","Std. Dev.","K-S"),
##       n.cgroup=c(1,1,3,2,1),
##       colheads=c("","","sG","HG","p","sG","HG","p"),
##       file="z.tex",
##       table.env=FALSE)

## cols <- as.list(rgb(red=c(0,0,182,146,146,219)/256, green=c(146,109,109,0,73,209)/256, blue=c(146,219,255,0,0,0)/256,alpha=0.4))

nicer <- function(p) {
    nicecow(p) +
      scale_fill_manual(values=c("F"=cols$F,"HG"=cols$HG,"sG"=cols$sG)) +
      scale_x_continuous(breaks=scales:::extended_breaks(n=3)) +
      theme(panel.spacing.x = unit(1, "lines"))
}
niceqq <- function(p) {
    nicecow(p) +
      scale_col_manual(values=c("F"=cols[1],"HG"=cols[5],"sG"=cols[2])) +
      scale_x_continuous(breaks=scales:::extended_breaks(n=3)) +
      theme(panel.spacing.x = unit(1, "lines"))
}
## nicer(pz1)



## save(results,file="results.RData")
## causal - main

pz1 <- plotz(results[class=="causal" & n==1000,]) + ggtitle("a: Z score, n=1000")
pz5 <- plotz(results[class=="causal" & n==5000,]) + ggtitle("b: Z score, n=5000")
pb1 <- plotb(results[class=="causal" & n==1000,]) + ggtitle("c: log OR, n=1000")
pb5 <- plotb(results[class=="causal" & n==5000,]) + ggtitle("d: log OR, n=5000")

myqq <- function(x,y,z) {
    if(!all(is.na(y))) {
        d1 <- as.data.frame(qqplot(x, y, plot.it=FALSE));
        d1$what="y"
    }
    if(!all(is.na(z))) {
        d2 <- as.data.frame(qqplot(x, z, plot.it=FALSE));
        d2$what="z"
    }
    if(all(is.na(y)))
        return(d2)
    if(all(is.na(z)))
        return(d1)
    rbind(d1,d2)
}

prepqq <- function(what) {
    v <- paste0(what,c("1",".f",".hg"))
    ret <- split(results,by=BY)  %>%  lapply(., function(dt) {
        tmp <- as.data.table(myqq(dt[[v[1]]], dt[[v[2]]], dt[[v[3]]]))
        tmp[what=="y",what:=v[[2]]]
        tmp[what=="z",what:=v[[3]]]
        cbind(tmp,dt[1,BY,with=FALSE])
    })  %>%
      rbindlist()  
    ret[what==v[2],variable:="F"]
    ret[what==v[3],variable:="HG"]
}

bqq <- prepqq("beta")
bqq$stat <- "beta"
bqz <- prepqq("z")
bqz$stat <- "Z"
bqq <- rbind(bqq,bqz)

plots <- list()
for(istat in c("beta","Z")) {
    for(i.n in c("1000","5000")) {
        plots <- c(plots,
                   list(ggplot(bqq[class=="causal" & stat==istat & n==i.n,],
                          aes(x=x,y=y,col=variable))  + geom_point() + geom_abline() +
                   facet_wrap(~lab,scales="free") +
                   labs(x="simGWAS beta",y="forward/HapGen beta")))
    }
}
        
ggplot(bqq[class=="causal" ,],
       aes(x=x,y=y,col=variable,pch=n))  + geom_point() + geom_abline() +
  facet_wrap(stat~lab,scales="free",nrow=4) +
  labs(x="simGWAS beta",y="forward/HapGen beta")
        
tmp <- bqq[class=="causal" & stat=="beta" & n=="5000",]
p <- ggplot(tmp) + geom_abline(linetype="dotted") +
  geom_point(aes(x=x,y=y,col=variable,pch=variable)) + 
  facet_wrap(~lab,scales="free") +
  labs(x="simGWAS beta",y="forward/HapGen beta") +
  ## scale_colour_manual("Method",values=c(HG=cols$HG,F=cols$F)) +
  scale_colour_manual("Method",values=c(HG="black",F="gray")) +
  scale_shape_discrete("Method") + 
  theme(legend.position="bottom")
ks.tests <- results[class=="causal", .(sg.hg=ks.test(beta1,beta.hg)$p.value,
                                       sg.f=ks.test(beta1,beta.f)$p.value,
                                       f.hg=ks.test(beta.hg,beta.f)$p.value),
                    by=BY]
tmp2 <- tmp[stat=="beta",.(minx=min(x),maxy=max(y)),by=BY]
ks.tests <- merge(ks.tests,tmp2,by=BY)
ks.tests[,label:=paste0("sG-HG: ",format.pval(sg.hg,digits=2),"\n",
                        " F-sG: ",format.pval(sg.f,digits=2),"\n",
                        " F-HG: ",format.pval(f.hg,digits=2))]
p <- p + geom_label(aes(x=minx,y=maxy,label=label),data=ks.tests,
              vjust="top",hjust="left",family="mono",size=3)
p
pdf("qqplots.pdf",height=8,width=8)
print(p)
dev.off()

## ## Fig 2 - unlinked
## tmp <- bqq[class!="causal" & n=="5000" & what!="beta.f" & what!="z.f",]
## tmp[,what:=sub(".hg","",what)]
## tmp[,scen:=sub("-.*","",lab)]
## p <- ggplot(tmp) + geom_point(aes(x=x,y=y,col=stat)) + geom_abline() +
##   facet_wrap(stat~scen,scales="free",nrow=2) +
##   labs(x="simGWAS beta",y="HapGen beta") +
##   scale_colour_manual("Method",values=c(beta=cols$HG,Z=cols$F)) +
##   theme(legend.position="bottom")
## ks.tests <- results[class=="causal", .(sg.hg=ks.test(beta1,beta.hg)$p.value),
##                     by=BY]
## tmp2 <- tmp[stat=="beta",.(minx=min(x),maxy=max(y)),by=BY]
## ks.tests <- merge(ks.tests,tmp2,by=BY)
## ks.tests[,label:=paste0("sG-HG: ",format.pval(sg.hg,digits=2),"\n")]
## p <- p + geom_label(aes(x=minx,y=maxy,label=label),data=ks.tests,
##               vjust="top",hjust="left",family="mono",size=3)
## p
## pdf("qqplots.pdf",height=8,width=8)
## print(p)
## dev.off()

## bqq <- prepqq("z")
## ggplot(bqq[class=="causal" & n==5000,],
##        aes(x=x,y=y,col=variable))  + geom_point() + geom_abline() + facet_wrap(~lab,scales="free") +
##   labs(x="simGWAS Z",y="forward/HapGen Z")

## z <- with(results[class=="causal" & n==5000,], myqq(beta1,beta.f,beta.hg))

##     ggplot(d) + geom_point(aes(x=x, y=y))

## ggplot(results[class=="causal" & n==5000,],aes(x=beta1,y=beta.f)) + geom_qq()


## w <- 1.2
## library(ggridges)
## pdf(file="ridge.pdf",height=11*w,width=9*w)
## plot_grid(pz1,pz5,pb1,pb5,nrow=2)
## dev.off()

## unlinked - supp
## results[class!="causal",lab:=sub("-.*","",lab)]
## pz1 <- plotz(results[class!="causal" & n==1000,]) + ggtitle("a: Z score, n=1000")
## pz5 <- plotz(results[class!="causal" & n==5000,]) + ggtitle("b: Z score, n=5000")
## pb1 <- plotb(results[class!="causal" & n==1000,]) + ggtitle("c: log OR, n=1000")
## pb5 <- plotb(results[class!="causal" & n==5000,]) + ggtitle("d: log OR, n=5000")


## Figure s2
plotz <- function(results) {
    tmp <- melt(results[,.(n,lab,z1,z.hg,z.f)],c("n","lab"))
    tmp[variable=="z1",variable:="sG"]
    tmp[variable=="z.hg",variable:="HG"]
    tmp[variable=="z.f",variable:="F"]
    tmp[,variable:=factor(variable,levels=c("F","sG","HG"))]
    setnames(tmp,"n","Sample.size")
    pz <- ggplot(tmp,aes(x=value,y=variable,fill=variable)) +
      geom_density_ridges() +
      facet_wrap(~lab,nrow=3) + #geom_vline(xintercept=0,linetype="dashed") +
      xlab("Z score") + ylab("") +
      theme(legend.position="none")
    nicer(pz)
## pz
}
plotb <- function(results) {
    results[,beta.mn:=mean(beta.hg),by=c("lab","class")]
    tmp <- melt(results[,.(n,lab,beta1,beta.hg,beta.f,beta.mn)],c("n","lab","beta.mn"))
    tmp[variable=="beta1",variable:="sG"]
    tmp[variable=="beta.hg",variable:="HG"]
    tmp[variable=="beta.f",variable:="F"]
    tmp[,variable:=factor(variable,levels=c("F","sG","HG"))]
    setnames(tmp,"n","Sample.size")
    pb <- ggplot(tmp,aes(x=value,y=variable,fill=variable)) + geom_density_ridges() + facet_wrap(~lab,nrow=3) + #geom_vline(aes(xintercept=beta.mn),linetype="dashed") +
      xlab("log odds ratio") + ylab("") +
  theme(legend.position="none")
    nicer(pb)
}
f <- function(s) {
    results[,scen:=sub("-.*","",lab)]
tmp <- results[scen==s & n==5000,]
tmp <- tmp[class=="unlinked",beta.f:=NA]
tmp <- tmp[class=="unlinked",z.f:=NA]
tmp <- tmp[,lab:=sub("-3","-0",lab)]
    list(z=plotz(tmp),b= plotb(tmp))
}

plots <- c(f("4"),f("5")) 
p <- plot_grid(plotlist=plots,ncol=4)
ggsave(p,file="ridge-4-5.pdf",height=9,width=9)

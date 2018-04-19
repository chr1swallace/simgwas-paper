library(ggplot2)
library(cowplot)
library(magrittr)
library(data.table)
library(randomFunctions)
d <- "/rds/user/cew54/hpc-work/simgwas"

## read in haplotypes
dref <- "/home/cew54/newscratch/Data/reference/1000GP_Phase3"
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
## XX <- new("SnpMatrix", as.matrix(freq))
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

## cv <- results[snp %in% dfsnps[unlinked,]$rs,]
## cv <- results[snp %in% CV[[1]],]
## G1 <- g1[1]
## with(cv,qqplot(z.hg,z1)); abline(0,1)
## with(cv,qqplot(z.hg,z2)); abline(0,1)
## with(cv,qqplot(beta.hg,beta1)); abline(0,1); abline(v=log(G1))
## with(cv,mean(beta1-G1))
## with(cv,hist(beta.hg)); abline(v=log(G1))
## with(cv,qqplot(beta.hg,beta2)); abline(0,1)
## with(cv,qqplot(beta1,be); abline(0,1)
## with(cv,qqplot(z1,z2)); abline(0,1)

## with(results,qqplot(beta.hg,beta1)); abline(0,1)
## with(results,qqplot(beta.hg,beta2)); abline(0,1)
## with(results,qqplot(z.hg,z1)); abline(0,1)
## with(results,qqplot(z.hg,z2)); abline(0,1)

## tmp <- results[!is.na(z.hg),ks.test(beta.hg,beta2),by="snp"]
## with(cv,ks.test(z.hg,z2)$p.value)

loadset <- function(f) {
    RESULTS <- vector("list",length(f))
    for(i in seq_along(f)) {
        (load(f[[i]]))
        results[,gsim:=paste(i,sim,sep="-")]
        RESULTS[[i]] <- results
    }
    rbindlist(RESULTS)
}

## read in results
## big number of sims on a small set of CVs
files <- list.files(d,pattern="spec",full=TRUE)
# files <- list.files(d,pattern="s100",full=TRUE)
## split files by category
fcat <- basename(files)  %>%  sub("-[0-9a-f]+.RData","",.)
files <- split(files,fcat)
files <- files[ grep("spec7|spec8",names(files),invert=TRUE) ]
message("files found: ",length(unlist(files)))
pb <- txtProgressBar(min=1,max=length(files),style=3)
RESULTS <- vector("list",length(files))
for(i in seq_along(files)) { # i=33
    setTxtProgressBar(pb, i)
    results <- loadset(files[[i]])
    SPECIAL <- sub(".*spec","",names(files)[i])
    CV <- switch(SPECIAL,
                 "1"="rs367629917.9654916.T.C", # common snp, small effect
                 "2"="rs376308868.9828660.G.C", # maf=0.02, larger effect
                 "3"=c("rs77603406.9413839.C.T", "rs377252712.9827703.T.C", "rs75434219.9671019.T.C"), # 3 unlinked
                 "4"= c("rs75112728.9830024.A.C", "X21.9680193.G.GA"), # 2 with r=-0.15
                 "5"=c("rs71247672.9723463.A.G", "rs79178122.9724174.G.A"), # 2 with r=0.5, similar maf
                 "6"=c("rs371462627.9695792.A.G", "X21.9695816.A.G"), # 2 with rs=0.8
                 "7"="rs367629917.9654916.T.C", # common snp, big effect
                 "8"=c("rs71247672.9723463.A.G", "rs79178122.9724174.G.A")) # 2 with rs=0.5
    g1 <- switch(SPECIAL,
                 "1"=1.1,
                 "2"=1.5,
                 "3"=c(1.3,1.2,1.1),
                 "4"=c(1.2,1/1.2),
                 "5"=c(1.2,1/1.2),
                 "6"=c(1.2,1/1.2),
                 "7"=2,
                 "8"=c(2,2))
##     if(length(CV)>1) {
##         ld <- LD[CV,CV]
##         diag(ld) <- 0
##         print(max(abs(ld)))
##     }
## }
    results[,z1:=beta1/v]
    results[,z2:=beta2/v]
    results[,z.hg:=beta.hg/v.hg]
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
    res <- results[snp %in% asnps , .(beta1=sort(beta1),beta2=sort(beta2),beta.hg=sort(beta.hg),
                                      z1=sort(z1),z2=sort(z2),z.hg=sort(z.hg),
                                      az1=sort(abs(z1)),az2=sort(abs(z2)),az.hg=sort(abs(z.hg))),
                   by=c("snp")]
    RESULTS[[i]] <- merge(res,tmp2,by="snp")
}

results <- rbindlist(RESULTS)
table(results$cat)
results[,lab:=gsub("sim-c.*-spec|-causal|-unlinked|-proxy","",snplab)]
results[,n:=gsub("sim-c|-spec.*","",cat)]
table(paste(results$cat,results$lab),results$class)

cols <- rgb(red=c(0,0,182,146,146,219)/256, green=c(146,109,109,0,73,209)/256, blue=c(146,219,255,0,0,0)/256,alpha=0.4)
library(ggridges)
nicer <- function(p) nicecow(p) + scale_fill_manual(values=c("1000"=cols[2],"5000"=cols[5]))

tmp <- melt(results[class=="causal" & lab!="7-1" & lab!="8-1" & lab!="8-2",.(n,lab,z1,z.hg)],c("n","lab"))
tmp[variable=="z1",variable:="simGWAS"]
tmp[variable=="z.hg",variable:="HAPGEN"]
setnames(tmp,"n","Sample.size")
pz <- ggplot(tmp,aes(x=value,y=variable,fill=Sample.size)) + geom_density_ridges() + facet_wrap(~lab,scales="free_x",nrow=4) + geom_vline(xintercept=0,linetype="dashed") + xlab("Z score") + ylab("") +
  theme(legend.position=c(0.75,0.1))
nicer(pz)
ggsave(nicer(pz),file="ridge-z.pdf",height=8,width=6)

results[,beta.mn:=mean(beta.hg),by=c("lab","class")]
tmp <- melt(results[class=="causal" & lab!="7-1" & lab!="8-1" & lab!="8-2",.(n,lab,beta1,beta.hg,beta.mn)],c("n","lab","beta.mn"))
tmp[variable=="beta1",variable:="simGWAS"]
tmp[variable=="beta.hg",variable:="HAPGEN"]
setnames(tmp,"n","Sample.size")
pz <- ggplot(tmp,aes(x=value,y=variable,fill=Sample.size)) + geom_density_ridges() + facet_wrap(~lab,scales="free_x",nrow=4) + geom_vline(aes(xintercept=beta.mn),linetype="dashed") + xlab("log odds ratio") + ylab("") +
  theme(legend.position=c(0.75,0.1))
ggsave(nicer(pz),file="ridge-beta.pdf",height=8,width=6)

tmp <- melt(results[class!="causal" & lab!="7-1" & lab!="8-1" & lab!="8-2",.(n,lab,z1,z.hg)],c("n","lab"))
tmp[variable=="z1",variable:="simGWAS"]
tmp[variable=="z.hg",variable:="HAPGEN"]
setnames(tmp,"n","Sample.size")
pz <- ggplot(tmp,aes(x=value,y=variable,fill=Sample.size)) + geom_density_ridges() + facet_wrap(~lab,scales="free_x",nrow=4) + geom_vline(xintercept=0,linetype="dashed") + xlab("Z score") + ylab("") +
  theme(legend.position="below")
nicer(pz)
ggsave(nicer(pz),file="ridge-z-unlinked.pdf",height=8,width=6)

results[,beta.mn:=mean(beta.hg),by=c("lab","class")]
tmp <- melt(results[class!="causal" & lab!="7-1" & lab!="8-1" & lab!="8-2",.(n,lab,beta1,beta.hg,beta.mn)],c("n","lab","beta.mn"))
tmp[variable=="beta1",variable:="simGWAS"]
tmp[variable=="beta.hg",variable:="HAPGEN"]
setnames(tmp,"n","Sample.size")
pz <- ggplot(tmp,aes(x=value,y=variable,fill=Sample.size)) + geom_density_ridges() + facet_wrap(~lab,scales="free_x",nrow=4) + geom_vline(aes(xintercept=beta.mn),linetype="dashed") + xlab("log odds ratio") + ylab("") +
  theme(legend.position="below")
ggsave(nicer(pz),file="ridge-beta-unlinked.pdf",height=8,width=6)



## pz1 <- ggplot(results[grepl("c1000",cat) & class=="causal",],aes(x=z.hg,y=z1)) + geom_point() +
##   facet_wrap(~lab,scales="free") + geom_abline() 
## nicer(pz1)

## pz5 <- ggplot(results[grepl("c5000",cat) & class=="causal",],aes(x=z.hg,y=z1)) + geom_point() +
##   facet_wrap(~lab,scales="free") + geom_abline() 
## nicecow(pz5)

## pz <- ggplot(results[class=="causal",],aes(x=z.hg,y=z2,col=n)) + geom_point() +
##   facet_wrap(~lab,scales="free") + geom_abline() 
## nicer(pz)

## pz <- ggplot(results[class!="causal",],aes(x=z.hg,y=z1,col=n)) + geom_point() +
##   facet_wrap(~lab,scales="free") + geom_abline() 
## nicer(pz)


## pz <- ggplot(results[class!="causal",],aes(x=z.hg,y=z1,col=n)) + geom_point() +
##   facet_wrap(~lab,scales="free") + geom_abline() 
## nicer(pz)


## pz5 <- ggplot(results[grep("c5000",cat),],aes(x=z.hg,y=z1,col=class)) + geom_point() +
##   facet_wrap(~snplab,scales="free") + geom_abline() 
## nicecow(pz5)

## pb1 <- ggplot(results[grepl("c1000",cat) & class=="causal",],aes(x=beta.hg,y=beta1)) + geom_point() +
##   facet_wrap(~snplab,scales="free") + geom_abline() + geom_vline(aes(xintercept=beta.marg)) +
##   geom_hline(aes(yintercept=beta.marg))
## nicecow(pb1)

## pb <- ggplot(results[class=="causal",],aes(x=beta.hg,y=beta1,col=n)) + geom_point() +
##   facet_wrap(~lab,scales="free") + geom_abline() + geom_vline(aes(xintercept=beta.marg)) +
##   geom_hline(aes(yintercept=beta.marg))
## nicecow(pb)

## pb5 <- ggplot(results[grep("c5000",cat),],aes(x=beta.hg,y=beta1,col=class)) + geom_point() +
##   facet_wrap(~snplab,scales="free") + geom_abline() + geom_vline(aes(xintercept=beta.marg)) +
##   geom_hline(aes(yintercept=beta.marg))
## nicecow(pb5)

## ggsave(nicecow(pz1),file="z1000.pdf",height=6,width=6)
## ggsave(nicecow(pb1),file="beta1000.pdf",height=6,width=6)
## ggsave(nicecow(pz5),file="z5000.pdf",height=6,width=6)
## ggsave(nicecow(pb5),file="beta5000.pdf",height=6,width=6)

## ##     with(res,qqplot(beta1,beta2)); abline(0,1,col="red")
## ##     with(res,qqplot(beta.hg,beta1)); abline(0,1,col="red")
## ##     with(res,qqplot(z.hg,z1)); abline(0,1,col="red")
## ##     ## res <- results[snp %in% CV,]
## ##     with(res[snp==CV[[1]],], qqplot(z.hg,z1)); abline(0,1,col="red")
## ##     #with(res[snp==CV[[1]],], qqplot(z.hg,z2)); abline(0,1,col="red")
## ##     with(res[snp==CV[[1]],], qqplot(beta.hg,beta1)); abline(0,1,col="red")

## ##     unlinked <- names(unlinked)
## ##     res <- results[snp %in% unlinked,]
## ##     with(res[snp==unlinked[[1]],], qqplot(z.hg,z1)); abline(0,1,col="red")
## ##     #with(res[snp==unlinked[[1]],], qqplot(z.hg,z2)); abline(0,1,col="red")
## ##     with(res[snp==unlinked[[1]],], qqplot(beta.hg,beta1)); abline(0,1,col="red")


## ##     tmp <- melt(results[snp %in% c(CV,proxy,unlinked),.(snp,beta1,beta2,beta.hg)],"snp")
## ##     tmp <- merge(tmp,tmp2,by="snp")

## ##     ggplot(tmp,aes(x=value-beta.marg))+ geom_histogram() + facet_grid(variable ~ snp)
    
## ##     mn <- results[snp %in% asnps , lapply(.SD,mean,na.rm=TRUE),by="snp"]
## ##     effects <- data.table(snp=asnps,beta.true=lg1,beta.marg=beta.marg,
## ##                           class=c(rep("cv",length(CV)), rep("proxy",sum(!is.na(proxy))),
## ##                                   rep("unlinked",sum(!is.na(unlinked)))))
## ##     effects <- merge(effects,dfsnps[,.(rs,maf)],by.x="snp",by.y="rs")
## ##     mn <- merge(mn,effects,by="snp")
## ##     mn[,stat:="mean"]
## ##     var <- results[snp %in% asnps, lapply(.SD,var),by="snp"]
## ##     var <- merge(var,effects,by="snp")
## ##     var[,stat:="var"]
## ##     RESULTS[[i]] <- rbind(mn,var)
## ##     RESULTS[[i]][,ncv:=length(CV)]
## ##     RESULTS[[i]][,sim:=i]
## ##     RESULTS[[i]][,ncse:=args$NCSE]
## ## }
## ## close(pb)



## ## files <- list.files(d,pattern="RData",full=TRUE)

## ## message("files found: ",length(files))
## ## pb <- txtProgressBar(min=1,max=length(files),style=3)
## ## RESULTS <- vector("list",length(files))
## ## for(i in seq_along(files)) { # i=33
## ##     setTxtProgressBar(pb, i)
## ##     (load(files[[i]]))
## ##     message(i, "\t",length(CV))
## ## ##     if(length(CV)>1) {
## ## ##         ld <- LD[CV,CV]
## ## ##         print(ld)
## ## ##     }
## ## ## }
## ##     head(results)
## ##     results[,z1:=beta1/v]
## ##     results[,z2:=beta2/v]
## ##     results[,z.hg:=beta.hg/v.hg]
## ##     results[,sim:=i]
## ##     proxy <- getproxy(LD[,CV,drop=FALSE],sqrt(0.4), sqrt(0.8))
## ##     unlinked <- get0(LD[,CV,drop=FALSE],0.1)
## ##     ## asnps <- c(CV,setdiff(results$snp,CV)) #dfsnps[c(setdiff(proxy,NA),setdiff(unlinked,NA)),]$rs)
## ##     asnps <- c(CV,dfsnps[c(setdiff(proxy,NA),setdiff(unlinked,NA)),]$rs)
## ##     lg1 <- c(log(g1),rep(0,length(asnps)-length(g1)))
## ##     if(length(asnps)>1) {
## ##         ld <- LD[asnps,asnps]
## ##         beta.marg <- ld %*% lg1
## ##     } else {
## ##         beta.marg <- lg1
## ##     }
## ##     mn <- results[snp %in% asnps , lapply(.SD,mean,na.rm=TRUE),by="snp"]
## ##     effects <- data.table(snp=asnps,beta.true=lg1,beta.marg=beta.marg,
## ##                           class=c(rep("cv",length(CV)), rep("proxy",sum(!is.na(proxy))),
## ##                                   rep("unlinked",sum(!is.na(unlinked)))))
## ##     effects <- merge(effects,dfsnps[,.(rs,maf)],by.x="snp",by.y="rs")
## ##     mn <- merge(mn,effects,by="snp")
## ##     mn[,stat:="mean"]
## ##     var <- results[snp %in% asnps, lapply(.SD,var),by="snp"]
## ##     var <- merge(var,effects,by="snp")
## ##     var[,stat:="var"]
## ##     RESULTS[[i]] <- rbind(mn,var)
## ##     RESULTS[[i]][,ncv:=length(CV)]
## ##     RESULTS[[i]][,sim:=i]
## ##     RESULTS[[i]][,ncse:=args$NCSE]
## ## }
## ## close(pb)

## ## data <- rbindlist(RESULTS)
## ## head(data)
## ## summary(data)
## ## setnames(data,"beta.marg.V1","beta.marg")

## ## mn <- data[stat=="mean",]
## ## mn <- melt(mn[!is.na(beta.hg),.(sim,ncse,beta.true,beta.marg,ncv,class,beta1,beta2,beta.hg)],c("sim","ncse","beta.true","beta.marg","ncv","class"))

## ## THR <- 0
## ## message("about given values")
## ## tmp <- mn[(beta.true-beta.marg)>=THR,
## ##           .(mn.beta=mean(value - beta.true,na.rm=TRUE),
## ##              var.beta=var(value-beta.true,na.rm=TRUE),
## ##              n=sum(!is.na(value))),by=c("variable","ncv","ncse","class")]
## ## ## head(tmp)
## ## tmp[order(class,ncv,ncse,variable),]

## ## message("about marg values")
## ## tmp <- mn[(beta.true-beta.marg)>=THR,
## ##           .(mn.beta=mean(value - beta.marg,na.rm=TRUE),
## ##             var.beta=var(value-beta.marg,na.rm=TRUE),
## ##             n=sum(!is.na(value))),by=c("variable","ncv","ncse","class")]
## ## ## head(tmp)
## ## tmp[order(class,ncv,ncse,variable),]

## ## message("VAR about marg values")
## ## v <- data[stat=="var" & !is.na(v.hg),.(ncv,ncse,v,v.hg)]
## ## v[,lapply(.SD,mean,na.rm=TRUE),by=c("ncv","ncse")]
## ## v[,lapply(.SD,median,na.rm=TRUE),by=c("ncv","ncse")]

## ## mn <- data[stat=="mean",]
## ## mn <- melt(mn[!is.na(beta.hg) & abs(beta.marg -beta.true)<0.05,
## ##               .(sim,ncse,beta.marg,ncv,beta1,beta2,beta.hg)],c("sim","ncse","beta.marg","ncv"))

## ## library(ggplot2)
## ## ggplot(mn,aes(x=value - beta.marg)) + geom_histogram() + facet_wrap(variable~ncv) + geom_vline(xintercept=0,col="red",linetype="dashed")
## ## mn[,.N,by=c("variable","ncv")]

## ## with(data[stat=="mean" & !is.na(beta.hg),],qqplot(beta1,beta.hg)); abline(0,1)
## ## with(data[stat=="mean" & !is.na(beta.hg),],qqplot(beta2,beta.hg)); abline(0,1)
## ## with(data[stat=="mean" & !is.na(beta.hg),],qqplot(z2,z.hg)); abline(0,1)

## ## with(mn[ncv==3 & variable=="beta.hg",], {
## ##     qqnorm(value - beta.marg)
## ##     qqline(value - beta.marg) })

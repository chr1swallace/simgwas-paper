library(ggplot2)
library(cowplot)
library(magrittr)
library(data.table)
library(randomFunctions)
d <- "/rds/user/cew54/hpc-work/simgwas"

## read in haplotypes
dref <- "/home/cew54/share/Data/reference/1000GP_Phase3"
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
files <- files[ grep("spec5|spec7|spec8",names(files),invert=TRUE) ]
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
## relabel spec6 as spec5
results[,cat:=sub("spec6","spec5",cat)]
table(results$cat)
results[,lab:=gsub("sim-c.*-spec|-causal|-unlinked|-proxy","",snplab)]
results[,n:=gsub("sim-c|-spec.*","",cat)]
table(paste(results$cat,results$lab),results$class)
results[,lab:=sub("6-","5-",lab)]

cols <- rgb(red=c(0,0,182,146,146,219)/256, green=c(146,109,109,0,73,209)/256, blue=c(146,219,255,0,0,0)/256,alpha=0.4)
library(ggridges)

nicer <- function(p) {
    nicecow(p) +
      scale_fill_manual(values=c("HG"=cols[5],"sG"=cols[2])) +
      scale_x_continuous(breaks=scales:::extended_breaks(n=3)) +
      theme(panel.spacing.x = unit(1, "lines"))
}
## nicer(pz1)


plotz <- function(results) {
tmp <- melt(results[,.(n,lab,z1,z.hg)],c("n","lab"))
tmp[variable=="z1",variable:="sG"]
tmp[variable=="z.hg",variable:="HG"]
setnames(tmp,"n","Sample.size")
pz <- ggplot(tmp,aes(x=value,y=variable,fill=variable)) + geom_density_ridges() + facet_wrap(~lab,scales="free_x",nrow=3) + #geom_vline(xintercept=0,linetype="dashed") +
  xlab("Z score") + ylab("") +
  theme(legend.position="none")
nicer(pz)
## pz
}
plotb <- function(results) {
    results[,beta.mn:=mean(beta.hg),by=c("lab","class")]
    tmp <- melt(results[,.(n,lab,beta1,beta.hg,beta.mn)],c("n","lab","beta.mn"))
    tmp[variable=="beta1",variable:="sG"]
    tmp[variable=="beta.hg",variable:="HG"]
    setnames(tmp,"n","Sample.size")
    pb <- ggplot(tmp,aes(x=value,y=variable,fill=variable)) + geom_density_ridges() + facet_wrap(~lab,scales="free_x",nrow=3) + #geom_vline(aes(xintercept=beta.mn),linetype="dashed") +
      xlab("log odds ratio") + ylab("") +
  theme(legend.position="none")
    nicer(pb)
}

## save(results,file="results.RData")
## causal - main

pz1 <- plotz(results[class=="causal" & n==1000,]) + ggtitle("a: Z score, n=1000")
pz5 <- plotz(results[class=="causal" & n==5000,]) + ggtitle("b: Z score, n=5000")
pb1 <- plotb(results[class=="causal" & n==1000,]) + ggtitle("c: log OR, n=1000")
pb5 <- plotb(results[class=="causal" & n==5000,]) + ggtitle("d: log OR, n=5000")

w <- 1.2
pdf(file="ridge.pdf",height=11*w,width=9*w)
plot_grid(pz1,pz5,pb1,pb5,nrow=2)
dev.off()

## unlinked - supp
results[class!="causal",lab:=sub("-.*","",lab)]
pz1 <- plotz(results[class!="causal" & n==1000,]) + ggtitle("a: Z score, n=1000")
pz5 <- plotz(results[class!="causal" & n==5000,]) + ggtitle("b: Z score, n=5000")
pb1 <- plotb(results[class!="causal" & n==1000,]) + ggtitle("c: log OR, n=1000")
pb5 <- plotb(results[class!="causal" & n==5000,]) + ggtitle("d: log OR, n=5000")

ggsave(plot_grid(pz1,pz5,pb1,pb5,nrow=2),file="ridge-unlinked.pdf",height=11,width=9)

library(data.table)
d <- "/rds/user/cew54/hpc-work/simgwas"

## read in haplotypes
library(data.table)
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
    unlist(res)
}
get0 <- function(x,upper) {
    mx <- apply(abs(x),1,max)
    res <- which(mx<upper)
    if(length(res)>=ncol(x)) 
        return(sample(res,ncol(x)))
    c(res,rep(NA,ncol(x)-length(res)))
}

cv <- results[snp %in% dfsnps[unlinked,]$rs,]
cv <- results[snp %in% CV[[1]],]
G1 <- g1[1]
with(cv,qqplot(z.hg,z1)); abline(0,1)
with(cv,qqplot(z.hg,z2)); abline(0,1)
with(cv,qqplot(beta.hg,beta1)); abline(0,1); abline(v=log(G1))
with(cv,mean(beta1-G1))
with(cv,hist(beta.hg)); abline(v=log(G1))
with(cv,qqplot(beta.hg,beta2)); abline(0,1)
with(cv,qqplot(beta1,be); abline(0,1)
with(cv,qqplot(z1,z2)); abline(0,1)

with(results,qqplot(beta.hg,beta1)); abline(0,1)
with(results,qqplot(beta.hg,beta2)); abline(0,1)
with(results,qqplot(z.hg,z1)); abline(0,1)
with(results,qqplot(z.hg,z2)); abline(0,1)

tmp <- results[!is.na(z.hg),ks.test(beta.hg,beta2),by="snp"]
with(cv,ks.test(z.hg,z2)$p.value)

## read in results
 files <- list.files(d,pattern="s1000",full=TRUE)
files <- list.files(d,pattern="RData",full=TRUE)
message("files found: ",length(files))
pb <- txtProgressBar(min=1,max=length(files),style=3)
RESULTS <- vector("list",length(files))
for(i in seq_along(files)) { # i=33
    setTxtProgressBar(pb, i)
    (load(files[[i]]))
    message(i, "\t",length(CV))
    if(length(CV)>1) {
        ld <- LD[CV,CV]
        print(ld)
    }
}

    head(results)
    results[,z1:=beta1/v]
    results[,z2:=beta2/v]
    results[,z.hg:=beta.hg/v.hg]
    results[,sim:=i]
    proxy <- getproxy(LD[,CV,drop=FALSE],sqrt(0.4), sqrt(0.8))
    unlinked <- get0(LD[,CV,drop=FALSE],0.1)
    ## asnps <- c(CV,setdiff(results$snp,CV)) #dfsnps[c(setdiff(proxy,NA),setdiff(unlinked,NA)),]$rs)
    asnps <- c(CV,dfsnps[c(setdiff(proxy,NA),setdiff(unlinked,NA)),]$rs)
    lg1 <- c(log(g1),rep(0,length(asnps)-length(g1)))
    if(length(asnps)>1) {
        ld <- LD[asnps,asnps]
        beta.marg <- ld %*% lg1
    } else {
        beta.marg <- lg1
    }
    mn <- results[snp %in% asnps , lapply(.SD,mean,na.rm=TRUE),by="snp"]
    effects <- data.table(snp=asnps,beta.true=lg1,beta.marg=beta.marg,
                          class=c(rep("cv",length(CV)), rep("proxy",sum(!is.na(proxy))),
                                  rep("unlinked",sum(!is.na(unlinked)))))
    effects <- merge(effects,dfsnps[,.(rs,maf)],by.x="snp",by.y="rs")
    mn <- merge(mn,effects,by="snp")
    mn[,stat:="mean"]
    var <- results[snp %in% asnps, lapply(.SD,var),by="snp"]
    var <- merge(var,effects,by="snp")
    var[,stat:="var"]
    RESULTS[[i]] <- rbind(mn,var)
    RESULTS[[i]][,ncv:=length(CV)]
    RESULTS[[i]][,sim:=i]
    RESULTS[[i]][,ncse:=args$NCSE]
}
close(pb)

data <- rbindlist(RESULTS)
head(data)
summary(data)
setnames(data,"beta.marg.V1","beta.marg")

mn <- data[stat=="mean",]
mn <- melt(mn[!is.na(beta.hg),.(sim,ncse,beta.true,beta.marg,ncv,class,beta1,beta2,beta.hg)],c("sim","ncse","beta.true","beta.marg","ncv","class"))

THR <- 0
message("about given values")
tmp <- mn[(beta.true-beta.marg)>=THR,
          .(mn.beta=mean(value - beta.true,na.rm=TRUE),
             var.beta=var(value-beta.true,na.rm=TRUE),
             n=sum(!is.na(value))),by=c("variable","ncv","ncse","class")]
## head(tmp)
tmp[order(class,ncv,ncse,variable),]

message("about marg values")
tmp <- mn[(beta.true-beta.marg)>=THR,
          .(mn.beta=mean(value - beta.marg,na.rm=TRUE),
            var.beta=var(value-beta.marg,na.rm=TRUE),
            n=sum(!is.na(value))),by=c("variable","ncv","ncse","class")]
## head(tmp)
tmp[order(class,ncv,ncse,variable),]

message("VAR about marg values")
v <- data[stat=="var" & !is.na(v.hg),.(ncv,ncse,v,v.hg)]
v[,lapply(.SD,mean,na.rm=TRUE),by=c("ncv","ncse")]
v[,lapply(.SD,median,na.rm=TRUE),by=c("ncv","ncse")]

mn <- data[stat=="mean",]
mn <- melt(mn[!is.na(beta.hg) & abs(beta.marg -beta.true)<0.05,
              .(sim,ncse,beta.marg,ncv,beta1,beta2,beta.hg)],c("sim","ncse","beta.marg","ncv"))

library(ggplot2)
ggplot(mn,aes(x=value - beta.marg)) + geom_histogram() + facet_wrap(variable~ncv) + geom_vline(xintercept=0,col="red",linetype="dashed")
mn[,.N,by=c("variable","ncv")]

with(data[stat=="mean" & !is.na(beta.hg),],qqplot(beta1,beta.hg)); abline(0,1)
with(data[stat=="mean" & !is.na(beta.hg),],qqplot(beta2,beta.hg)); abline(0,1)
with(data[stat=="mean" & !is.na(beta.hg),],qqplot(z2,z.hg)); abline(0,1)

with(mn[ncv==3 & variable=="beta.hg",], {
    qqnorm(value - beta.marg)
    qqline(value - beta.marg) })

#!/usr/bin/env Rscript
library(data.table)
library(randomFunctions)
library(simGWAS)
library(mvtnorm)
library(corpcor)
## working directory
d <- "/rds/user/cew54/hpc-work/simgwas"

args <- getArgs(defaults=list(N=1000,NEFF=1,NSIM=100,
                              file.ldd="/home/cew54/newscratch/Data/reference/lddetect/EUR/fourier_ls-chr22.bed",
                              file.vcf="/home/cew54/newscratch/Data/reference/UK10K/chr22.bcf.gz",
                              file.cv=file.path(d,"chrcv.csv")),
                numeric=c("N","NEFF","NSIM"))
                
cv <- fread(args$file.cv)

## ldblocks
ldd <- fread(args$file.ldd)

## split bcf by ldblocks
tmp <- tempfile()
ldd[,blocknum:=1:.N]

ldd[,comm:=paste0("/home/cew54/localc/bin/bcftools view ",args$file.vcf,
                  " --min-af 0.01:minor --max-alleles 2 --min-alleles 2 ",
                  " -r chr",22,":",start,"-",stop," -Ov ")] # -o ",tmp)]

gethap <- function(i) {
    y=fread(ldd$comm[i])
    ha <- simGWAS:::vcf2haps(as.matrix(y[,-c(1:9)]))
    rownames(ha) <- paste0("pos",y$POS)
    t(ha)
}
cor2 <- function (x) {
    1/(NROW(x) - 1) * crossprod( scale(x, TRUE, TRUE) )
}

message("simulating for ",nrow(ldd),"blocks")
simv <- simz <- simbeta <- vector("list",nrow(ldd))
for(i in 1:nrow(ldd)) {

    cat(i,"\r")
    h <- gethap(i)

    use <- apply(h,2,var)>0
    if(any(!use))
        h <- h[,use,drop=FALSE]
    freq <- as.data.frame(h+1)
    freq$Probability <- 1/nrow(freq)
    LD <- cor2(h)

    pcv <- cv[block==i,]$pos
    if(!length(pcv)) {
        pcv <- colnames(h)[1]
        g1 <- 0
    } else {
        pcv <- paste0("pos",pcv)
        g1 <- cv[block==i,]$eff
    }
    FP <- make_GenoProbList(snps=colnames(h),W=pcv,freq=freq)

## method 1 - simulate Z scores and adjust by simulated variance to get beta
    EZ <- est_statistic(N0=args$N, # number of controls
                            N1=args$N, # number of cases
                            snps=colnames(h), # column names in freq of SNPs for which Z scores should be generated
                            W=pcv, # causal variants, subset of snps
                            gamma1=g1, # odds ratios
                            freq=freq, # reference haplotypes
                            GenoProbList=FP) # FP above
    simz[[i]] <- t(rmvnorm(n = args$NSIM, mean = EZ, sigma = LD))
    tmpsimv <- sim_vbeta(N0=args$N, # number of controls
                         N1=args$N, # number of cases
                         snps=colnames(h), # column names in freq of SNPs for which Z scores should be generated
                         W=pcv, # causal variants, subset of snps
                         gamma1=g1, # odds ratios
                         freq=freq, # reference haplotypes
                         GenoProbList=FP,
                         nsim=args$NSIM)
    simv[[i]] <- 1/do.call("cbind",tmpsimv)
    simbeta[[i]] <- simz[[i]] * sqrt(simv[[i]])
}

simv <- do.call("rbind",simv)
simz <- do.call("rbind",simz)
simbeta <- do.call("rbind",simbeta)

## now do whatever with the output

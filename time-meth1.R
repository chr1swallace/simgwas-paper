#!/usr/bin/env Rscript
library(randomFunctions)
args <- getArgs(defaults=list(N=1000,NEFF=1,NSIM=100),numeric=c("N","NEFF","100"))

CV <- c("X21.9557710.C.CAA", "rs372513062.9429998.GA.G", "X21.9579507.T.C", 
        "rs75467879.9828538.G.A", "rs376042200.9507863.A.C", "X21.9476405.C.A")[args$NEFF]
g1 <- c(1.5, 1.8, 1.2, 1.8, 1.2, 1.5)[args$NEFF]

d <- "/rds/user/cew54/hpc-work/simgwas"

## FIRST 5000 snps in Africans only
dref <- "/home/cew54/newscratch/Data/reference/1000GP_Phase3"
fhg <- file.path(d,"input")#tempfile(tmpdir=dtmp)

## read in haplotypes
library(data.table)
h <- fread(paste0(fhg,".hap"))
h <- as.matrix(h)
snps <- fread(paste0(fhg,".leg"))

## devtools::load_all("~/RP/simGWAS")
library(simGWAS)
library(mvtnorm)
library(corpcor)
h1 <- h[, seq(1,ncol(h)-1,by=2)] #[,samples$GROUP=="AFR"]
h2 <- h[, seq(2,ncol(h),by=2)] #[,samples$GROUP=="AFR"]
freq <- as.data.frame(t(cbind(h1,h2))+1)
snps$rs <- make.names(snps$id)
colnames(freq) <- snps$rs
use <- apply(freq,2,var)>0
freq <- freq[,use,drop=FALSE]
dfsnps <- snps[use,,drop=FALSE]
## dfsnps$maf <- colMeans(freq-1)
LD <- cor(freq)
LD <- as.matrix(make.positive.definite(LD))
freq$Probability <- 1/nrow(freq)
FP <- make_GenoProbList(snps=dfsnps$rs,W=CV,freq=freq)

## method 1 - simulate Z scores and adjust by simulated variance to get beta
EZ <- est_statistic(N0=args$N, # number of controls
                    N1=args$N, # number of cases
                    snps=dfsnps$rs, # column names in freq of SNPs for which Z scores should be generated
                    W=CV, # causal variants, subset of snps
                    gamma1=log(g1), # odds ratios
                    freq=freq, # reference haplotypes
                    GenoProbList=FP) # FP above
simz <- t(rmvnorm(n = args$NSIM, mean = EZ, sigma = LD))
simv <- sim_vbeta(N0=args$N, # number of controls
                N1=args$N, # number of cases
                snps=dfsnps$rs, # column names in freq of SNPs for which Z scores should be generated
                W=CV, # causal variants, subset of snps
                gamma1=log(g1), # odds ratios
                freq=freq, # reference haplotypes
                GenoProbList=FP,
                nsim=args$NSIM)
simv <- 1/do.call("cbind",simv)
simbeta <- simz* sqrt(simv)


## ## method 2 - calculate an expected beta and simulate about that,
## valt <- est_vbeta(N0=args$NCTL, # number of controls
##                   N1=args$NCSE, # number of cases
##                   snps=dfsnps$rs, # column names in freq of SNPs for which Z scores should be generated
##                   W=dfsnps$rs[CV], # causal variants, subset of snps
##                   gamma1=log(g1), # odds ratios
##                   freq=freq, # reference haplotypes
##                   GenoProbList=FP) # FP above
## Ebeta <- EZ * valt
## simbeta2 <- lapply(1:args$NSIM, function(i) {
##     vv <- tcrossprod(sqrt(simv[,i]))
##     t(rmvnorm(n=1,mean=Ebeta,sigma=vv * LD))
## })
## simbeta2 <- do.call("cbind",simbeta2)

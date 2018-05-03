#!/usr/bin/env Rscript
library(data.table)
library(randomFunctions)
library(simGWAS)
library(mvtnorm)
library(corpcor)
source("~/DIRS.txt") # SIMGWAS, REFDATA, BINDIR locations
args <- getArgs(defaults=list(N=1000,NEFF=1,NSIM=100,chr=22),
                numeric=c("N","NEFF","NSIM"))
file.ldd=file.path(REFDATA,paste0("lddetect/EUR/fourier_ls-chr",args$chr,".bed"))
file.vcf=file.path(REFDATA,paste0("UK10K/chr",args$chr,".bcf.gz"))
                

g1 <- c(1.5, 1.8, 1.2, 1.8, 1.2, 1.5)[args$NEFF]

## working directory
d <- SIMGWAS

## ldblocks
ldd <- fread(file.ldd)

## split bcf by ldblocks
tmp <- tempfile()
ldd[,blocknum:=1:.N]
ldd[,comm:=paste0(BINDIR,"/bcftools view ",file.vcf,
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
simv <- simz <- simbeta <- vector("list",nrow(ldd))
for(i in 1:nrow(ldd)) {

    h <- gethap(i)

    use <- apply(h,2,var)>0
    h <- h[,use,drop=FALSE]
    freq <- as.data.frame(h+1)
    freq$Probability <- 1/nrow(freq)
    LD <- readRDS(file.path(d,"uk10k",paste0("ld-",args$chr,"-",i,".rds"))) #cor2(h)
    CV <- sample(colnames(h),1)
    
    FP <- make_GenoProbList(snps=colnames(h),W=CV,freq=freq)

    ## method 1 - simulate Z scores and adjust by simulated variance to get beta
    EZ <- simGWAS:::est_statistic(N0=args$N, # number of controls
                        N1=args$N, # number of cases
                        snps=colnames(h), # column names in freq of SNPs for which Z scores should be generated
                        W=CV, # causal variants, subset of snps
                        gamma.W=log(g1), # odds ratios
                        freq=freq, # reference haplotypes
                        GenoProbList=FP) # FP above
    simz[[i]] <- rmvnorm(n = args$NSIM, mean = EZ, sigma = LD)
    tmpsimv <- simulated_vbeta(N0=args$N, # number of controls
                         N1=args$N, # number of cases
                         snps=colnames(h), # column names in freq of SNPs for which Z scores should be generated
                         W=CV, # causal variants, subset of snps
                         gamma.W=log(g1), # odds ratios
                         freq=freq, # reference haplotypes
                         GenoProbList=FP,
                         nrep=args$NSIM)
    simv[[i]] <- 1/tmpsimv
    simbeta[[i]] <- simz[[i]] * sqrt(simv[[i]])
}


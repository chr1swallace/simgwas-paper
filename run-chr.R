#!/usr/bin/env Rscript
library(data.table)
library(randomFunctions)
library(simGWAS)
library(mvtnorm)
library(corpcor)
args <- getArgs(defaults=list(N=1000,NEFF=1,NSIM=100,
                              file.ldd="/home/cew54/newscratch/Data/reference/lddetect/EUR/fourier_ls-chr22.bed",
                              file.vcf="/home/cew54/newscratch/Data/reference/UK10K/chr22.bcf.gz"),
                numeric=c("N","NEFF","NSIM"))
                

g1 <- c(1.5, 1.8, 1.2, 1.8, 1.2, 1.5)[args$NEFF]

## working directory
d <- "/rds/user/cew54/hpc-work/simgwas"

## ldblocks
ldd <- fread(args$file.ldd)

## split bcf by ldblocks
tmp <- tempfile()
ldd[,blocknum:=1:.N]
## ldd[,hapfile:=paste0(tmp,"_",blocknum,".hap.gz")]
## ldd[,comm:=paste0("/home/cew54/localc/bin/bcftools view ",args$file.vcf," --min-af 0.01:minor --max-alleles 2 --min-alleles 2 -r chr",22,":",start,"-",stop," -Ou |/home/cew54/localc/bin/bcftools convert --haplegendsample ",tmp,"_",blocknum)]
## system.time({
##     system(ldd$comm[i])
##     x=as.matrix(fread(paste("zcat",ldd$hapfile[i],"| sed 's/\\*//g'")))
## })

ldd[,comm:=paste0("/home/cew54/localc/bin/bcftools view ",args$file.vcf,
                  " --min-af 0.01:minor --max-alleles 2 --min-alleles 2 ",
                  " -r chr",22,":",start,"-",stop," -Ov ")] # -o ",tmp)]
## system.time({
##     y=fread(ldd$comm[i])
##     yl <- y[,1:9]
##     y <- as.matrix(y[,-c(1:9)])
##     y1 <- matrix(as(substr(y,1,1),"numeric"),nrow(y))
##     y2 <- matrix(as(substr(y,3,3),"numeric"),nrow(y))
##     h <- cbind(y1,y2)
## })

## system.time({
##     y=fread(ldd$comm[i])
##     y <- as.matrix(y[,-c(1:9)])
##     ha <- vcf2haps(y[,-c(1:9)])
## })

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
    LD <- cor2(h)
    ## system.time(LD1 <- corpcor::make.positive.definite(LD))
    ## system.time(LD2 <- Matrix::nearPD(LD,corr=TRUE,ensureSymmetry=TRUE))
    CV <- sample(colnames(h),1)
    
    FP <- make_GenoProbList(snps=colnames(h),W=CV,freq=freq)

    ## method 1 - simulate Z scores and adjust by simulated variance to get beta
    EZ <- est_statistic(N0=args$N, # number of controls
                        N1=args$N, # number of cases
                        snps=colnames(h), # column names in freq of SNPs for which Z scores should be generated
                        W=CV, # causal variants, subset of snps
                        gamma1=log(g1), # odds ratios
                        freq=freq, # reference haplotypes
                        GenoProbList=FP) # FP above
    simz[[i]] <- t(rmvnorm(n = args$NSIM, mean = EZ, sigma = LD))
    tmpsimv <- sim_vbeta(N0=args$N, # number of controls
                         N1=args$N, # number of cases
                         snps=colnames(h), # column names in freq of SNPs for which Z scores should be generated
                         W=CV, # causal variants, subset of snps
                         gamma1=log(g1), # odds ratios
                         freq=freq, # reference haplotypes
                         GenoProbList=FP,
                         nsim=args$NSIM)
    simv[[i]] <- 1/do.call("cbind",tmpsimv)
    simbeta[[i]] <- simz[[i]] * sqrt(simv[[i]])
}


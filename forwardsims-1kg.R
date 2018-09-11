#!/usr/bin/env Rscript
library(randomFunctions)
library(magrittr)
library(data.table)
args <- getArgs(defaults=list(NCSE=5000,NCTL=5000,NSIM=1000,NCV=1,SPECIAL=2),numeric=c("NCSE","NCTL","NSIM","NCV"))
source("~/DIRS.txt")

d <- SIMGWAS

## FIRST 5000 snps in Africans only
dref <- file.path(REFDATA,"1000GP_Phase3")
fhg <- file.path(d,"input")#tempfile(tmpdir=dtmp)
## crange <- paste0("7-",ncol(h)+6)
## system(paste("zcat",
##             file.path(dref,"chr21_afr.hap.gz"),
##             "|head -n 5000 | cut -d' ' -f",crange,
##             "> ", paste0(fhg,".hap")))
## system(paste("zcat",
##              file.path(dref,"1000GP_Phase3_chr21.legend.gz"),
##              "| head -n 5001 > ", paste0(fhg,".leg")))


## read in haplotypes
## h <- fread(file.path(d,"example/ex.haps"))[1:1000,]
## snps <- fread(file.path(d,"example/ex.leg"))
## map <- fread(file.path(d,"example/ex.map"))
## wh <- which(c(0,diff(map[["Genetic_Map(cM)"]]))>0.1)
h <- fread(paste0(fhg,".hap"))
h <- as.matrix(h)
snps <- fread(paste0(fhg,".leg"))
samples <- fread(file.path(dref,"1000GP_Phase3.sample"))
## map <- fread(file.path(d,"example/ex.map"))
## wh <- which(c(0,diff(map[["Genetic_Map(cM)"]]))>0.1)

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
## XX <- new("SnpMatrix", as.matrix(freq))S""
freq$Probability <- 1/nrow(freq)
sum(freq$Probability)

#' @title inverse logit function
#' @param x value to evaluate inverse logit function at
#' @return value of inverse logit function at x
inv.logit.fn <-function(x) return(exp(x)/(1+exp(x)))

if(args$SPECIAL=="0") {
    CV=sample(which(dfsnps$AFR > 0.2 & dfsnps$AFR < 0.8),args$NCV)
#c(474) #sample(1:nrow(snps),2)
## g1 <- rep(2,args$NCV) #sample(c(1.2,1.5,1.8),args$NCV,replace=TRUE)
    g1 <- sample(c(1.1,1.2,1.3),args$NCV,replace=TRUE)
} else {
    source("~/Projects/simgwas/getspec.R")
    sCV <- CV
    CV <- match(CV,dfsnps$rs) # make numeric
}

## probs that someone is a case, depending on genotype at CV
loga <- log(0.5) # common disease
a <- 0.5 # common disease

## repeatedly sample two rows from h, keeping only the CV, and assign individual as case or control depending on genotype
N <- 20000 # simulate more than we need and throw away
N.target <- args$NCSE
hCV <- freq[,CV,drop=FALSE]-1
nh <- nrow(hCV)

fun <- function() {
    X <- as.matrix(hCV[ sample(1:nh,N,replace=TRUE),,drop=FALSE ] +
                   hCV[ sample(1:nh,N,replace=TRUE),,drop=FALSE])
    X <- cbind(1,X)
    p <- inv.logit.fn(X %*% matrix(log(c(a,g1)),ncol=1))
    Y <- runif(N) < p
    use <- c(which(Y==TRUE)[1:N.target],
             which(Y==FALSE)[1:N.target])
    X <- X[use,,drop=FALSE]
    Y <- as.numeric(Y[use])
    m <- lapply(seq_along(CV), function(i)  {
        m <- glm(Y ~ X[,i+1],family="binomial")
        c(beta=coefficients(m)[2],v=vcov(m)[2,2])
    })  %>%  do.call("rbind",.)  %>%  as.data.table()
    setnames(m,c("beta.f","v.f"))
    m[,snp:=sCV]
    ## names(m) <- paste(c("beta","z"), rep(seq_along(CV),each=2),sep=".f")
    m
}

fresults <- replicate(args$NSIM,fun(),simplify=FALSE)  %>%  do.call("rbind",.)
fresults[,cat:=paste0("sim-c",args$NCSE,"-spec",args$SPECIAL)]
fresults[,z.f:=beta.f/sqrt(v.f)]
save(fresults,file=file.path(d,paste0("forward-c",args$NCSE,"-spec",args$SPECIAL,".RData")))



library(randomFunctions)
NCSE <- 1000
NCTL <- 1000

## read in haplotypes
library(data.table)
h <- fread("/home/cew54/example/ex.haps")[1:1000,]
snps <- fread("~/example/ex.leg")
map <- fread("~/example/ex.map")
wh <- which(c(0,diff(map[["Genetic_Map(cM)"]]))>0.1)

library(simGWAS)
library(mvtnorm)
library(corpcor)
freq <- as.data.frame(t(h)+1)
snps$maf <- colMeans(freq-1)
colnames(freq) <- snps$rs
use <- apply(freq,2,var)>0
freq <- freq[,use,drop=FALSE]
dfsnps <- snps[use,,drop=FALSE]
LD <- cor(freq)
LD <- as.matrix(make.positive.definite(LD))
## XX <- new("SnpMatrix", as.matrix(freq))
freq$Probability <- 1/nrow(freq)
sum(freq$Probability)
CV=c(474) #sample(1:nrow(snps),2)
g1 <- c(1.5)
FP <- make_GenoProbList(snps=dfsnps$rs,W=dfsnps$rs[CV],freq=freq)

## method 1 - simulate Z scores and adjust by expected variance to get beta
z <- est_statistic(N0=NCTL, # number of controls
              N1=NCSE, # number of cases
              snps=dfsnps$rs, # column names in freq of SNPs for which Z scores should be generated
              W=dfsnps$rs[CV], # causal variants, subset of snps
              gamma1=log(g1), # odds ratios
              freq=freq, # reference haplotypes
              GenoProbList=FP) # FP above
## cs <- col.summary(XX)
## LD <- snpStats::ld(XX, XX, stat = "R", symmetric = TRUE)
sim_z_score <- t(rmvnorm(n = 10000, mean = z, sigma = LD))
p <- 2*pnorm(abs(sim_z_score),lower.tail=FALSE)
Var.data.cc <- function(f, N, s) {
  1 / (2 * N * f * (1 - f) * s * (1 - s))
}
v <- Var.data.cc(dfsnps$maf,NCSE+NCTL,NCSE/(NCSE+NCTL))
beta <- sim_z_score * sqrt(v)

## method 2 - calculate an expected beta and simulate about that,
## using a variance-covariance matrix of V * LD
expbeta <- z * sqrt(v)
V <- sqrt(v) %*% t(sqrt(v))
simbeta <- t(rmvnorm(n=1000,mean=expbeta,sigma=V * LD))
simz <- simbeta/sqrt(v)

## method 3 - method 1 but with alternate estimate of E(v(beta))
valt <- est_vbeta(N0=NCTL, # number of controls
              N1=NCSE, # number of cases
              snps=dfsnps$rs, # column names in freq of SNPs for which Z scores should be generated
              W=dfsnps$rs[CV], # causal variants, subset of snps
              gamma1=log(g1), # odds ratios
              freq=freq, # reference haplotypes
              GenoProbList=FP) # FP above
betaalt <- sim_z_score * valt

## method 4 - method 2 but with alternate estimate of E(v(beta))
expbetaalt <- z * valt
V <- valt %*% t(valt)
simbetaalt <- t(rmvnorm(n=1000,mean=expbetaalt,sigma=V * LD))

## method 5 - method 1 but with randomly sampled v(beta)
devtools::load_all("~/RP/simGWAS")
VB <- sim_vbeta(N0=NCTL, # number of controls
              N1=NCSE, # number of cases
              snps=dfsnps$rs, # column names in freq of SNPs for which Z scores should be generated
              W=dfsnps$rs[CV], # causal variants, subset of snps
              gamma1=log(g1), # odds ratios
              freq=freq, # reference haplotypes
              GenoProbList=FP,
              nsim=10000)
vb <- do.call("cbind",VB)
betab <- sim_z_score * sqrt(1/vb)
head(dv <- data.frame(maf=pmin(dfsnps$maf,1-dfsnps$maf),oneoverv=1/v,v=v,valt=valt,vb1=VB[[1]],vb2=VB[[2]],hg1=hg.se[,1],hg2=hg.se[,2]))
ggplot(dv,aes(x=sqrt(v),y=1/sqrt(vb1),col=maf)) + geom_point() + geom_smooth() + geom_abline()
    

## method 6 - method 2 but with randomly sampled v(beta)
simbetab <- lapply(1:100, function(i) {
    vv <- sqrt(1/VB[[i]]) %*% t(sqrt(1/VB[[i]]))
    t(rmvnorm(n=1,mean=expbetaalt,sigma=vv * LD))
})
simbetab <- do.call("cbind",simbetab)
    
## do the same for hapgen
fhg <- function() {
    tmp <- tempfile(tmpdir=".")
    tmp2 <- tempfile(tmpdir=".")
    system(paste("~/localc/bin/hapgen2 -m ~/example/ex.map -l ~/example/ex.leg -h ~/example/ex.haps -o",
                 tmp,
                 "-dl", dfsnps$position[CV[1]], 1, g1[1], g1[1]^2,
                 dfsnps$position[CV[2]], 1, g1[2], g1[2]^2,
                 "-n",NCTL,NCSE,
                 " -t ~/example/ex.tags -no_haps_output"))
    system(paste("~/localc/bin/snptest -data ",paste0(tmp,".controls.gen"),paste0(tmp,".controls.sample"),paste0(tmp,".cases.gen"),paste0(tmp,".cases.sample"),
                 "-o",tmp2,
                 "-frequentist 1 -method score -pheno pheno"))
    ret <- fread(tmp2)[rsid %in% colnames(freq),.(rsid,frequentist_add_pvalue,frequentist_add_beta_1,frequentist_add_se_1)]
    ## cleanup
    unlink(list.files(".",pattern=tmp,full=TRUE))
    unlink(list.files(".",pattern=tmp2,full=TRUE))
    return(ret)
}

## library(parallel)
results <- replicate(10,fhg(),simplify=FALSE)
#results <- mclapply(1:20,function(i) fhg(),mc.cores=3)

hg.beta <- sapply(results,function(x) x$frequentist_add_beta_1)
hg.se <- sapply(results,function(x) x$frequentist_add_se_1)
hg.p <- sapply(results,function(x) x$frequentist_add_pvalue)
hg <- results[[1]]

## compare var(beta)
head(dv <- data.frame(maf=pmin(dfsnps$maf,1-dfsnps$maf),oneoverv=1/v,v=v,valt=valt,vb1=VB[[1]],vb2=VB[[2]],hg1=hg.se[,1],hg2=hg.se[,2]))
library(ggplot2)
ggplot(dv,aes(x=sqrt(v),y=sqrt(1/vb1),col=maf)) + geom_point() + geom_smooth() + geom_abline()
ggplot(dv,aes(x=valt,y=sqrt(1/vb1),col=maf)) + geom_point() + geom_smooth() + geom_abline()

ggplot(dv,aes(x=sqrt(v),y=valt,col=maf)) + geom_point() + geom_smooth() + geom_abline()
ggplot(dv,aes(x=sqrt(v),y=sqrt(vb2),col=maf)) + geom_point() + geom_smooth() + geom_abline()
ggplot(dv,aes(x=sqrt(v),y=hg1,col=maf)) + geom_point() + geom_smooth() + geom_abline()
ggplot(dv,aes(x=hg1,y=1/sqrt(vb1),col=maf)) + geom_point() + geom_smooth() + geom_abline()

## compare

## beta at causal var
qqplot(betab[CV[[1]],], hg.beta[CV[[1]],]); abline(0,1)
qqplot(simbetab[CV[[1]],], hg.beta[CV[[1]],]); abline(0,1)
qqplot(betab[CV[[1]],], simbetab[CV[[1]],]); abline(0,1)


df <- data.frame(maf=pmin(dfsnps$maf,1-dfsnps$maf),
                 hgz=apply(hg.beta/hg.se,1,mean),
                 hgzmn=apply(hg.beta/hg.se,1,min),
                 hgzmx=apply(hg.beta/hg.se,1,max),
                 hgse=apply(hg.se,1,mean),
                 hgsemn=apply(hg.se,1,min),
                 hgsemx=apply(hg.se,1,max),
                 sgz=rowMeans(sim_z_score),
                 sgzmn=apply(sim_z_score,1,min),
                 sgzmx=apply(sim_z_score,1,max),
                 sgse=sqrt(v),
                 sgsealt=valt,
                 hg=rowMeans(hg.beta),
                 hgmn=apply(hg.beta,1,min),
                 hgmx=apply(hg.beta,1,max),
                 sg=rowMeans(beta),
                 sgmn=apply(beta,1,min),
                 sgmx=apply(beta,1,max),
                 sgalt=rowMeans(betaalt),
                 sgmnalt=apply(betaalt,1,min),
                 sgmxalt=apply(betaalt,1,max))

library(ggplot2)
## z + crosshairs
ggplot(df,aes(x=hgz,y=sgz,ymin=sgzmn,ymax=sgzmx,xmin=hgzmn,xmax=hgzmx,col=maf)) + geom_pointrange() + geom_errorbarh() + geom_smooth(method="lm") + geom_abline()

## beta + crosshairs
ggplot(df,aes(x=hg,y=sg,ymin=sgmn,ymax=sgmx,xmin=hgmn,xmax=hgmx,col=maf)) + geom_pointrange() + geom_errorbarh() + geom_smooth(method="lm") + geom_abline()

## se + errorbars
library(cowplot)
p1 <- ggplot(df,aes(x=hgse,y=sgse,xmin=hgsemn,xmax=hgsemx,col=maf)) + geom_point() + geom_errorbarh() + geom_smooth(method="lm") + geom_abline()
p2 <- ggplot(df,aes(x=hgse,y=sgsealt,xmin=hgsemn,xmax=hgsemx,col=maf)) + geom_point() + geom_errorbarh() + geom_smooth(method="lm") + geom_abline()
plot_grid(p1,p2)

sum((hg.se-sqrt(v))^2,na.rm=TRUE)
sum((hg.se-valt)^2,na.rm=TRUE)

sum((df$hgse-df$sgse)^2,na.rm=TRUE)
sum((df$hgse-df$sgsealt)^2,na.rm=TRUE)
sum((df$hg-df$sg)^2,na.rm=TRUE)
sum((df$hg-df$sgalt)^2,na.rm=TRUE)

sum((hg.beta - expbeta)^2,na.rm=TRUE)
sum((hg.beta - expbeta)^2,na.rm=TRUE)
sum((hg.beta - expbetaalt)^2,na.rm=TRUE)

ggplot(df,aes(x=maf,y=hgse,ymin=hgsemn,ymax=hgsemx,col=maf)) + geom_pointrange() + geom_smooth(method="lm") 

ggplot(df,aes(x=maf,y=hgz-sgz,col=maf)) + geom_point() #+ geom_smooth(method="lm") 






identical(snps$rs,hg$rs)
plot(snps$position,-log10(p[,1]))
points(snps$position,-log10(hg$frequentist_add_pvalue),col="red")
points(snps$position,-log10(p[,2]),col="grey")
points(snps$position,-log10(p[,3]),col="grey")



summary(apply(beta,2,min))
summary(beta[,1])
summary(hg$frequentist_add_beta)
plot(hg$frequentist_add_beta_1,beta[,1])
points(hg$frequentist_add_beta_1,apply(beta,1,min),col="grey")

cbind(head(sqrt(v)),head(hg.se))
cbind(head(beta[,1:4]),head(hg.beta[,1:4]))
summary(-log10(apply(p,2,min)))
summary(-log10(hg$frequentist_add_pvalue))


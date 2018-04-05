#!/usr/bin/env Rscript
library(data.table)
library(randomFunctions)
library(simGWAS)

## working directory
d <- "/rds/user/cew54/hpc-work/simgwas"

## PEFF = prob there is at least one CV in an ldblock
## EXPN = expected number of CV per ld block, given there is at least one CV
args <- getArgs(defaults=list(EXPN=1.5,PEFF=0.5,
                              CHR=22,
                              file.ldd="/home/cew54/newscratch/Data/reference/lddetect/EUR/fourier_ls-chr22.bed",
                              file.vcf="/home/cew54/newscratch/Data/reference/UK10K/chr22.bcf.gz",
                              file.cv=file.path(d,"chrcv.csv")),
                numeric=c("PEFF","EXPN"))

## ldblocks
ldd <- fread(args$file.ldd)

## split bcf by ldblocks
tmp <- tempfile()
ldd[,blocknum:=1:.N]

## number of cvs per block
ldd[,any.cv:=rbinom(nrow(ldd),1,args$PEFF)]
ldd[,ncv:=as.integer(0)]
ldd[any.cv>0,ncv:=rpois(sum(any.cv>0),args$EXPN)]

ldd[,comm:=paste0("/home/cew54/localc/bin/bcftools view ",args$file.vcf,
                  " --min-af 0.05:minor --max-alleles 2 --min-alleles 2 ",
                  " -r chr",args$CHR,":",start,"-",stop," -Ov ")] # -o ",tmp)]
gethap <- function(i) {
    y=fread(ldd$comm[i])
    ha <- simGWAS:::vcf2haps(as.matrix(y[,-c(1:9)]))
    rownames(ha) <- paste0("pos",y$POS)
    t(ha)
}
cor2 <- function (x) {
    1/(NROW(x) - 1) * crossprod( scale(x, TRUE, TRUE) )
}

block <- cv <- eff <- vector("list",nrow(ldd))
for(i in which(ldd$ncv>0)) {
    cat(i,"\r")
    h <- gethap(i)
    block[[i]] <- rep(i,ldd$ncv[i])
    cv[[i]] <- sample(colnames(h),ldd$ncv[i])
    eff[[i]] <- sample(c(1.1,1.2,1.3),ldd$ncv[i],replace=TRUE)
}

cv <- data.table(block=unlist(block),
                 pos=as.numeric(sub("pos","",unlist(cv))),
                 eff=log(unlist(eff)))

fwrite(cv,file=args$file.cv)

## make gstr for hapgen2
cv[,str:=paste(pos,1,eff,eff^2,sep=" ")]
cvstr <- paste(cv$str,collapse=" ")
cat(cvstr,file=sub(".csv",".txt",args$file.cv))


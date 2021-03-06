getspec <- function(x) {
    CV <- switch(x,
             "1"="rs367629917.9654916.T.C", # common snp, small effect
             "2"="rs376308868.9828660.G.C", # maf=0.02, larger effect
             "3"=c("rs77603406.9413839.C.T", "rs377252712.9827703.T.C", "rs75434219.9671019.T.C"), # 3 unlinked
             "4"= c("rs75112728.9830024.A.C", "X21.9680193.G.GA"), # 2 with r=-0.15
             ## "5"=c("rs71247672.9723463.A.G", "rs79178122.9724174.G.A"), # 2 with r=0.5, similar maf
             "5"=c("rs371462627.9695792.A.G", "X21.9695816.A.G"), # 2 with rs=0.8
             ## "7"="rs367629917.9654916.T.C", # common snp, big effect
             ## "8"=c("rs71247672.9723463.A.G", "rs79178122.9724174.G.A")) # 2 with r=0.5, similar maf
             )
    g1 <- switch(x,
             "1"=1.1,
             "2"=1.5,
             "3"=c(1.3,1.2,1.1),
             "4"=c(1.2,1/1.2),
             ## "5"=c(1.2,1/1.2),
             "5"=c(1.2,1/1.2),
             ## "7"=2,
             ## "8"=c(2,2))
             )
    return(list(CV=CV,g1=g1))
}

l <- getspec(args$SPECIAL)
CV <- l$CV
g1 <- l$g1

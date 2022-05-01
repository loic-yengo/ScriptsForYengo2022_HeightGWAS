library(meta)
library(parallel)
args   <- commandArgs(trailingOnly = T)
chunck <- args[1]

cat(paste0("\tChunk #",chunck,"...\n"))
cat("\tLoading data...\n")

setwd("/home/l.yengo/loic/GIANT2020/trans-ethnic-analyses")
cc   <- c(rep("character",3),rep("numeric",5))
afr  <- read.table("meta-analyses/AA.ma", header=TRUE,stringsAsFactors=FALSE,colClasses=cc)
eas  <- read.table("meta-analyses/EAS.ma",header=TRUE,stringsAsFactors=FALSE,colClasses=cc)
sas  <- read.table("meta-analyses/SAS.ma",header=TRUE,stringsAsFactors=FALSE,colClasses=cc)
his  <- read.table("meta-analyses/HIS.ma",header=TRUE,stringsAsFactors=FALSE,colClasses=cc)
eur  <- read.table("meta-analyses/EUR.ma",header=TRUE,stringsAsFactors=FALSE,colClasses=cc)

Bmat <- cbind(EUR=eur[,"beta"],AFR=afr[,"beta"],EAS=eas[,"beta"],SAS=sas[,"beta"],HIS=his[,"beta"])
Smat <- cbind(EUR=eur[,"se"],AFR=afr[,"se"],EAS=eas[,"se"],SAS=sas[,"se"],HIS=his[,"se"])
Pmat <- cbind(EUR=eur[,"p"],AFR=afr[,"p"],EAS=eas[,"p"],SAS=sas[,"p"],HIS=his[,"p"])
Fmat <- cbind(EUR=eur[,"freq"],AFR=afr[,"freq"],EAS=eas[,"freq"],SAS=sas[,"freq"],HIS=his[,"freq"])
Nmat <- cbind(EUR=eur[,"N"],AFR=afr[,"N"],EAS=eas[,"N"],SAS=sas[,"N"],HIS=his[,"N"])
pops <- colnames(Bmat)

M    <- nrow(Bmat)
nblc <- 100
grps <- as.numeric(cut(1:M,breaks = quantile(1:M,probs=seq(0,1,len=nblc+1)),include.lowest=TRUE))

cn <- c("nstudies","ancestry","freq","n",
        "beta_fixed","se_fixed","p_fixed",
        "beta_random","se_random","p_random",
        "Isq","Q","pv_Q")

runMeta <- function(i){
  b <- Bmat[i,]
  s <- Smat[i,]
  p <- Pmat[i,]
  f <- Fmat[i,]
  N <- Nmat[i,]
  l <- which(!is.na(b) & !is.na(s))
  n <- length(l)
  if(n==0){
    rs <- cbind.data.frame(nstudies=0,anc=NA,freq=NA,n=NA,beta_fixed=NA,se_fixed=NA,p_fixed=NA,beta_random=NA,se_random=NA,p_random=NA,Isq=NA,Q=NA,pv_Q=NA)
  }else{
    b  <- b[l]
    s  <- s[l]
    p  <- p[l]
    f  <- f[l]
    N  <- N[l]
    anc<- paste(sort(pops[l]),collapse="+")
    if(n==1){
      rs <- cbind.data.frame(nstudies=n,anc=anc,freq=f,n=N,beta_fixed=b,se_fixed=s,p_fixed=p,beta_random=b,se_random=s,p_random=p,Isq=0,Q=0,pv_Q=1)
    }else{
      m  <- metagen(TE=b,seTE=s)
      Nx <- sum(N)
      fx <- sum(f*N)/Nx
      rs <- cbind.data.frame(nstudies=n,anc=anc,freq=fx,n=Nx,beta_fixed=m$TE.fixed,se_fixed=m$seTE.fixed,p_fixed=m$pval.fixed,beta_random=m$TE.random,se_random=m$seTE.random,p_random=m$pval.random,Isq=m$I2,Q=m$Q,pv_Q=m$pval.Q)
    }
  }
  colnames(rs) <- cn
  return(rs)
}

index <- which(grps==chunck)
cat("\tRunning meta-analysis...\n")
rmeta <- do.call("rbind",mclapply(index,runMeta,mc.cores=10))
dt    <- cbind.data.frame(afr[index,c("SNP","A1","A2")],rmeta)

cat("\tWriting results...\n")
save(list="dt",file=paste0("Rmeta-results/grp",chunck,".RData"))
cat("\t[done].\n")



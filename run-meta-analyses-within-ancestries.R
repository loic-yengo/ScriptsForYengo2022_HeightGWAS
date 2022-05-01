setwd("/home/l.yengo/loic/GIANT2020/trans-ethnic-analyses")
M  <- 1385132
cc <- c(rep("character",3),rep("numeric",5))
fe_meta_analysis <- function(studies,outputname){
  nstudies <- length(studies)
  Bmat     <- matrix(NA,nrow=M,ncol=nstudies)
  SEmat    <- matrix(NA,nrow=M,ncol=nstudies)
  Nmat     <- matrix(NA,nrow=M,ncol=nstudies)
  Fmat     <- matrix(NA,nrow=M,ncol=nstudies)
  var_y    <- rep(NA,nstudies)
  for(i in 1:nstudies){
    cat(paste0("\tReading ",studies[i],"..."))
    tmp       <- read.table(studies[i],header=TRUE,stringsAsFactors=FALSE,colClasses=cc)
    Bmat[,i]  <- tmp[,"beta"]
    SEmat[,i] <- tmp[,"se"]
    Nmat[,i]  <- tmp[,"N"]
    Fmat[,i]  <- tmp[,"freq"]
    var_y[i]  <- median( tmp[,"N"] * 2 * tmp[,"freq"] * (1-tmp[,"freq"]) * tmp[,"se"] * tmp[,"se"], na.rm=TRUE)
    cat(paste0("E[var(Y)] = ",round(var_y[i],4),".\n"))
  }
  
  cat("\tCorr of frequencies\n")
  print(cor(Fmat,use="pairwise"))
  
  cat("\tCorr of BETA's\n")
  print(cor(Bmat,use="pairwise"))

  slopeMat <- matrix(NA,nrow=nstudies,ncol=nstudies)
  for(i in 1:nstudies){
    for(j in 1:nstudies){
      slopeMat[i,j] <- ifelse(i==j,1,coef(lm(Bmat[,i]~Bmat[,j]))[2])
    }
  }
  cat("\tSlope of BETA's\n")
  print(slopeMat)

  wgtr  <- 1/(SEmat^2)
  wgts  <- wgtr / rowSums(wgtr,na.rm=TRUE)
  BETA  <- rowSums(Bmat * wgts,na.rm=TRUE)
  SE    <- 1/sqrt(rowSums(wgtr,na.rm=TRUE))
  Z     <- BETA/SE
  CHISQ <- Z*Z
  P     <- pchisq(q=CHISQ,df=1,lower.tail=FALSE)
  N     <- rowSums(Nmat,na.rm=TRUE)
  FREQ  <- rowSums(Nmat*Fmat,na.rm=TRUE)/N
  
  ## filters
  f <- which(SE==0 | SE==Inf | abs(BETA)==Inf | N==0 | N==Inf | is.na(N) | is.na(BETA) | is.na(SE) | is.na(P) | is.na(FREQ) | FREQ==0 | FREQ==1)
  if(length(f)>0){
    print("Issue here!")
    BETA[f] <- SE[f] <- P[f] <- N[f] <- FREQ[f] <- NA
  }
  tmp[,"freq"] <- FREQ
  tmp[,"beta"] <- BETA
  tmp[,"se"]   <- SE
  tmp[,"p"]    <- P
  tmp[,"N"]    <- N
  write.table(tmp,outputname,row.names=FALSE,quote=FALSE)
}

sourceDir <- "/gpfs/gpfs01/polaris/Q0286/loic/GIANT2020/HM3/"

eur <- paste0(sourceDir,c("23andMe_EUR_int_height_female.HM3_v2.txt",
                          "23andMe_EUR_int_height_male.HM3_v2.txt",
                          "HEIGHT_EUR_ALL_FINAL_allchr_HM3_v2.txt"))

eas <- paste0(sourceDir,c("23andMe_EAS_int_height_female.HM3_v2.txt",
                          "23andMe_EAS_int_height_male.HM3_v2.txt",
                          "HEIGHT_EAS_1KGP3_withOut_UKB_FINAL_allchr_HM3_v2.txt"))

sas <- paste0(sourceDir,c("23andMe_SAS_int_height_female.HM3_v2.txt",
                          "23andMe_SAS_int_height_male.HM3_v2.txt",
                          "HEIGHT_SAS_1KGP3_withOUTUKB_FINAL_allchr_HM3_v2.txt"))

his  <- paste0(sourceDir,c("23andMe_HA_int_height_female.HM3_v2.txt",
                           "23andMe_HA_int_height_male.HM3_v2.txt",
                           "HEIGHT_HIS_1KGP3_FINAL_allchr_HM3_v2.txt"))

aa   <- paste0(sourceDir,c("23andMe_AA_int_height_female.HM3_v2.txt",
                           "23andMe_AA_int_height_male.HM3_v2.txt",
                           "HEIGHT_AA_1KGP3_FINAL_allchr_HM3_v2.txt"))

eur2 <- paste0(sourceDir,c("23andMe_EUR_int_height_female.HM3_v2.txt",
                           "23andMe_EUR_int_height_male.HM3_v2.txt",
                           "HEIGHT_EUR_ALL_withoutUKB_FINAL_allchr_HM3_v2.txt",
                           "HEIGHT_UKB_v3EUR_maf0001_nosibs_LMM_inf_HM3_v2.txt"))

giant2020 <- paste0(sourceDir,c("HEIGHT_EUR_ALL_withoutUKB_FINAL_allchr_HM3_v2.txt",
                                "HEIGHT_UKB_v3EUR_maf0001_withsibs_LMM_inf_HM3_v2.txt"))

print("Analysis 1")
#fe_meta_analysis(eur,"meta-analyses/EUR.ma")

print("Analysis 2")
#fe_meta_analysis(eas,"meta-analyses/EAS.ma")

print("Analysis 3")
#fe_meta_analysis(sas,"meta-analyses/SAS.ma")

print("Analysis 4")
#fe_meta_analysis(his,"meta-analyses/HIS.ma")

print("Analysis 5")
#fe_meta_analysis(aa,"meta-analyses/AA.ma")

#print("Analysis 6 -- naive")
#fe_meta_analysis(paste0("meta-analyses/",dir("meta-analyses")),"meta-analyses/ALL.ma")

print("Analysis 7 -- no UKB sibs")
#fe_meta_analysis(eur2,"meta-analyses/EUR_noUKBsibs.ma")

print("Analysis 8 - nosibs + UKB")
#fe_meta_analysis(giant2020,"meta-analyses/GIANT2020_inclUKB_v2.ma")

### New analyses >23APR2020
eas <- c(paste0(sourceDir,c("23andMe_EAS_int_height_female.HM3_v2.txt",
                            "23andMe_EAS_int_height_male.HM3_v2.txt")),
         "../download/formatted/GIANT2020noUKB_EAS.ma")

afr <- c(paste0(sourceDir,c("23andMe_AA_int_height_female.HM3_v2.txt",
                            "23andMe_AA_int_height_male.HM3_v2.txt")),
         "../download/formatted/GIANT2020noUKB_AFR.ma")

sas <- c(paste0(sourceDir,c("23andMe_SAS_int_height_female.HM3_v2.txt",
                            "23andMe_SAS_int_height_male.HM3_v2.txt")),
         "../download/formatted/GIANT2020noUKB_SAS.ma")

his <- c(paste0(sourceDir,c("23andMe_HA_int_height_female.HM3_v2.txt",
                            "23andMe_HA_int_height_male.HM3_v2.txt")),
         "../download/formatted/GIANT2020noUKB_HIS.ma")

eur <- c(paste0(sourceDir,c("23andMe_EUR_int_height_female.HM3_v2.txt",
                            "23andMe_EUR_int_height_male.HM3_v2.txt",
                            "HEIGHT_UKB_v3EUR_maf0001_nosibs_v2_LMM_inf_HM3_v2.txt")),
         "../download/formatted/GIANT2020noUKB_EUR.ma")
all <- c(eas,sas,afr,his,eur)

fe_meta_analysis(eas,"meta-analyses/EAS.ma")
fe_meta_analysis(sas,"meta-analyses/SAS.ma")
fe_meta_analysis(afr,"meta-analyses/AFR.ma")
fe_meta_analysis(his,"meta-analyses/HIS.ma")
#fe_meta_analysis(eur,"meta-analyses/EUR.ma")
#fe_meta_analysis(all,"meta-analyses/ALL.ma")




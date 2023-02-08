### Winner's curse correction
WC_correction <- function(BETA,        # Effect size 
                        SE,          # Standard Error
                        alpha=5e-8){ # Significance threshold
  Q     <- qchisq(alpha,df=1,lower.tail=FALSE)
  cs    <- sqrt(Q)
  
  bias  <- function(betaTrue,betaObs,se){
    z   <- betaTrue/se
    num <- dnorm(z-cs)-dnorm(-z-cs)
    # changed here from "-" to "+" 
    # Thanks Prof. Yoichiro Kamatani for pointing this out
    den <- pnorm(z-cs) + pnorm(-z-cs) 
    return(betaObs-betaTrue + se * num/den)
  }
  
  solveBetaTrue <- function(betaObs,se){
    md <- uniroot(f=function(b) bias(b,betaObs,se),lower=-100,upper=100)
    return(md$root)
  }
  
  ## Height
  BETA_corrected  <- sapply(1:length(BETA),function(i) solveBetaTrue(BETA[i],SE[i]))
  return(BETA_corrected)
}

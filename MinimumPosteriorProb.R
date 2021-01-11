# Transform pvalue in minimum posterior probability (probability of the H0 of being true) . Sellke (2001) and Held (2010)

# Minimum Bayes Factor of H0 to H1 given the p-value. 
# BF = -e*p*log(p) when p<1/e and 1 otherwise.
bf <- function(p, NaN.2.zero=T){
  
  bf <- (p<1/exp(1)&p>0)*(-exp(1)*p*log(p)) + (p>1/exp(1))*1
  
  if(NaN.2.zero){ bf[is.na(bf)] <- 0 }
  return(bf)
}

# Minimum Posterior Probability given the Bayes Factor:
mmp <- function(bf, PH0=0.5, selfloop=F) {
  mmp_matrix <- (1+(bf*PH0/(1-PH0))^-1)^-1
  
  if(is.matrix(bf)){
    if(!selfloop) {diag(mmp_matrix) <- 0}
  }
  return(mmp_matrix)
  }
  
# Functions to compute different types of random matrices.

# Sample 0-1 values from a given distribution
sample01 <- function(N , distr , param , normalize01=F ){
  if( !(distr %in% c('rand.adj','rand.adj.gnp','omog','unif','beta','norm','lnorm','poisson'))){stop("Non ammissible distribution. Choose among 'rand.adj','rand.adj.gnp','omog','unif', beta','norm','lnorm','poisson'")}
  
  switch(distr,
         
         rand.adj = {         
           x = sample(c(0,1) , N , replace = T)
         },
         
         rand.adj.gnp = {
           if(length(param)!=1){stop("Wrong number of parameter: Insert threshold p")}
           y = runif( N )
           x = numeric(N)
           x[which(y < param)] = 1
         },
         
         omog = {
           if(length(param)!=1){stop("The parameter must be 0<p<1")}
           x = rep(param , N)
         },
         
         unif = {
           if(length(param)!=2){stop("Wrong number of parameter (insert min and max bounds of Uniform Distribution")}
           x = runif( N )
         },
         
         beta = {
           if(length(param)==1){stop("Wrong number of parameter")}
           x = rbeta( N , param[1] , param[2])
         },
         
         norm = {
           if(length(param)==1){stop("Wrong number of parameter")}
           x = rnorm(N,param[1],param[2])
         },
         
         lnorm = {
           if(length(param)==1){stop("Wrong number of parameter")}
           x = rlnorm(N,param[1],param[2])
         },
         
         poisson = {
           if(length(param)!=1){stop("Only one parameter is required")}
           
           x = rpois(N,param)
           
         }
         
  )
  
  p <- x
  
  if(normalize01==T & (distr=='norm'|distr=='lnorm'|distr=='poisson')){p <- (x-min(x))/(max(x)-min(x))} #ci sarà 0 e 1...
  
  return( p )
  
}

#Generate Random Matrices (sym==symmetric , loop==diag!=0 , wigner==gaussian normalized matrix)
#INPUT: N = numel 
#       distr = selected distribution among 'rand.adj','rand.adj.gnp','omog','unif', beta','norm','lnorm'"
#       param = parameter of the distribution

#library(algstat)
upper <- function(v){
  n <- (1+sqrt(1+8*length(v)))/2
  mat <- matrix(0,n,n)
  mat[upper.tri(mat, diag=FALSE)] <- v
  return(mat)
}

lower <- function(v){
  n <- (1+sqrt(1+8*length(v)))/2
  mat <- matrix(0,n,n)
  mat[lower.tri(mat, diag=FALSE)] <- v
  return(mat)
}

symat <- function(v){
  n <- (1+sqrt(1+8*length(v)))/2
  mat <- matrix(0,n,n)
  mat[lower.tri(mat, diag=FALSE)] <- v
  mat[upper.tri(mat)] <-  t(mat)[upper.tri(mat)]
  return(mat)
}

random_sym_mat <- function(N,distr,param,normalize01=F){
  A <- matrix(0, N, N)
  upperN <- (N^2-N)/2

  x <- sample01(upperN,distr,param,normalize01)
  A <- symat(x)
  
  return(A)
}

random_sym_given_mat <- function(N,vector){
  A <- matrix(0, N, N)
  upperN <- (N^2-N)/2
  
  x <- sample(vector,upperN,replace = T)
  A <- symat(x)
  
  return(A)
}


random_sym_loop_mat <- function(N,distr,param,normalize01=F){
  A <- matrix(0, N, N)
  upperN <- (N^2-N)/2
  
  x <- sample01(upperN,distr,param,normalize01)
  diagA <- sample01(N,distr,param,normalize01)
  A <- upper(x) + lower(x) + diag(diagA)
  
  return(A)
}

randomat <- function(N,distr,param,normalize01=F){
  #each p_ij one with probability p sampled from U(0,1)
  A <- sample01(N^2,distr,param,normalize01) %>% matrix(N,N)
  
  return(A)
}

randomat_noloop <- function(N,distr,param,normalize01=F){
  #each p_ij one with probability p sampled from U(0,1)
  A <- sample01(N^2,distr,param,normalize01) %>% matrix(N,N)
  diag(A)=0
  
  return(A)
}

randomat_rect <- function(N,M,distr,param,normalize01=F){
  #each p_ij one with probability p sampled from U(0,1)
  A <- sample01(N*M,distr,param,normalize01) %>% matrix(N,M)
  return(A)
}

random_diff_beta <- function(N,a.vec,b.vec){
  
  upperN <- (N^2-N)/2
  A <- matrix(0, N, N)
  param.list = cbind( a.vec , b.vec )
  
  x <- lapply(1:upperN, function(t) rbeta( 1 , param.list[t,1] , param.list[t,2]) ) %>% unlist 
  
  A <- symat(x)
  Ex <- a.vec/(a.vec+b.vec)
  
  return(list(A=A,Ex=Ex))

  }

random_diff_norm <- function(N,param){
  
  upperN <- (N^2-N)/2
  A <- matrix(0, N, N)
  a <- numeric()
  b <- numeric()
  x <- numeric()
  Ex <- numeric()
    
  for(i in 1:upperN){
    
    a[i] = sample(param[,1],1)
    b[i] = sample(param[,2],1)
    x[i] <- rnorm(1,a[i],b[i])
    Ex[i] <- a[i]
  }
  A <- x %>% matrix(.,N,N)
  
  return(list(Ex=Ex,A=A))
}

mathresh <- function(W,threshold=0.75, type='mag', diag.value=NULL){
  
  if(type=='mag'){wnew <- (W > threshold)*1}
  if(type=='min'){wnew <- (W < threshold)*1}
  if(!is.null(diag.value)){ diag(wnew) <- diag.value }
  return(wnew)
  
}  




# This script contains the functions to compute the fuzzy descriptors for the nodes of an undirected network
# INPUT  --> P = matrix of probabilities of existence
#        --> analytic = T if the fuzzy degree must be from Poisson-Distribution rather than sampling FuzzyProb
#        --> Nsamples = number of random sampling from the FuzzyProb (if analytic==F)
# OUTPUT --> Pdnode = Dataframe with probabilities of a node having a certain value of the descriptor

library(poibin)

####################################################################################################################################
#~~~~~~~~~~~~~~~~~ ANCILLARY FUNCTIONS ~~~~~~~~~~~~~~~~~#
####################################################################################################################################.

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

# Thresholding probability/weighted adjacency matrix
mathresh <- function(W,threshold=0.75, type='mag', diag.value=NULL){
  
  if(type=='mag'){wnew <- (W > threshold)*1}
  if(type=='min'){wnew <- (W < threshold)*1}
  if(!is.null(diag.value)){ diag(wnew) <- diag.value }
  return(wnew)
  
}  

# Sampling a probability matrix, obtaining a list of binary adjacency matrix
sampling.P <- function(P, Nsamples){
  if(dim(P)[1] != dim(P)[2] | is.null(dim(P)[1])) {stop('The probability matrix is not square')}
  N <- dim(P)[1]
  lapply(1:Nsamples, function(x) (random_sym_mat(N,'unif', c(0,1) ) < P)*1 )
}

####################################################################################################################################.
#~~~~~~~~~~~~~~~~~ FUZZY FRAMEWORK FUNCTIONS ~~~~~~~~~~~~~~~~~#
####################################################################################################################################.

# Fuzzy Degree --------
# This function compute the fuzzy degree for the nodes of a network
# INPUT  --> P = matrix of probabilities of existence
#        --> analytic = T if the fuzzy degree must be from Poisson-Distribution rather than sampling FuzzyProb
#        --> Asample = List of sampled adjacency matrices (N x N each)
# OUTPUT --> Pdnode = Dataframe with probabilities of a node having a certain degree (degrees on rows and nodes on columns)

get_fuzzy_Degree <- function(P, analytic=T, Asample=NULL){
  
  if(analytic) {
    N <- dim(P)[1]
  } else{
    N <- dim(Asample[[1]])[1]
  }
  
  Pdnode <- matrix(NA,N,N)
  
  if(analytic) {
    # ANALYITIC (with poibin library, more efficient!)
    for (i in seq(1, N)) {
      for (k in seq(1, N)) {
        Pdnode[k, i] = dpoibin(kk = k - 1 , pp = P[i, ]) # TODO : ATTENZIONE! Probabilità che il nodo i abbia grado NON DIRETTO uguale a k!!! Non valido per CCM e TE che hanno matrici P asimmetriche
        
      }
    }
    
  } else{
    # NUMERIC (without poibin library)
    cat('Sampling Done!\n')
    simuls.node <- lapply(Asample, function(x) rowSums(x)) # sum the "successes" for each node (degree)
    n.links <- matrix(unlist(simuls.node) , Nsam , N , byrow = T) # matrix of degree for each node (column) for each sample (row)
    
    # Pdnode_table <- apply(n.links, 2, function(x) table(x)/Nsam)
    # Pdnode <- t( plyr::ldply(Pdnode_table, rbind) )
    hist.node.raw <- apply(n.links , 2 , function(x) hist(x, plot = F, breaks = seq(-1, N - 1))) # histogram of degrees for each node (along columns). Done instead of the previous two lines in order to have get results also for deterministic networks
    Pdnode <- matrix(unlist(lapply(hist.node.raw, "[[" , 3)) , N , N , byrow = T) # link distribution in N classes
    Pdnode <- t(Pdnode)
  }
  
  Pdnode <- as.data.frame(Pdnode)
  names(Pdnode) <- 1:N
  Pdnode <- cbind('degree' = seq(0, N-1), Pdnode)
  Pdnode[Pdnode < 1e-8] <- 0
  
  return(Pdnode)
}

# Fuzzy Clustering Coefficient ------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Given a network with N nodes, this function computes for each node, the probability that the (local) clustering coefficient is 
# equal to the value in the support. 
# This function takes the sampled adjacency matrices and computes the transitivity for each node of each sampled network. 
# Then, it compute the frequency of the transitivity values for every node.

# INPUT: 
# Asample = list of sampled adjacency matrix (N x N each)
# OUTPUT:
# sim.cl.pr = Probability of a node having clustering coefficent equal to values in the support

library(igraph)
get_fuzzy_Clustering_numerical <- function(Asample){
  
  N <- dim(Asample[[1]])[1]
  Nsim <- length(Asample)
  sim.cl <- lapply( Asample , function(x) transitivity( graph_from_adjacency_matrix( x, mode = 'undirected'), 'local') )
  sim.cl.df <- data.frame( matrix(sim.cl %>% unlist, ncol=N, byrow = T) ); names(sim.cl.df) <- 1:N #paste('n',seq(1,N),sep='')
  
  #apply(sim.cl.df,2,function(x) table(x,useNA = 'ifany') )
  sim.cl.values <- sort(unique(unlist(sim.cl.df)), na.last = F)
  
  # In the following line, the NA value (no possibility for triangles) is separated from the rest. sim.cl.pr compute the frequence of all possible clustering coefficient values for every node.
  sim.cl.pr <- lapply( 1:N , function(n) c( sum(is.na(sim.cl.df[,n]))/Nsim, ( lapply(sim.cl.values, function(x) sum(sim.cl.df[,n]==x , na.rm = T)/Nsim)  %>% unlist)[-1] ) )
  sim.cl.pr <- data.frame( matrix(sim.cl.pr %>% unlist, ncol=N, byrow = F) ) 
  names(sim.cl.pr) <- 1:N
  sim.cl.pr <- cbind('clustering'=sim.cl.values, sim.cl.pr)
  
  # in the following sim.cl.values[3] corresponds to the first clustering coefficient value, different from NA and zero. So the x axis is drawn symmetrically
  # if(length(sim.cl.values)>2 ){
  #   matplot( c(-sim.cl.values[3], sim.cl.values[-1]) , t(sim.cl.pr) , pch = 16 , type='b', lty = 1,  col = rainbow(N) , main = 'Local Clustering Distributions\nNumeric' , xlab = "Clustering Coefficient" , ylab = "Probability" , xaxt='n' )#, ylim=c(0,1));
  #   #legend('topright' , legend=paste('Node-',1:N,sep=''),col=rainbow(N), lty=1, cex=0.6 , bty ='n' , ncol = ceiling(N/3) )
  #   axis(side=1, at=c(-sim.cl.values[3], sim.cl.values[-1]) , label=sim.cl.values %>% round(3))
  #}
  
  return(sim.cl.pr)
  
}

# Theoretical support of the clustering coefficient probaility distribution (and comparison with the simulated one)
clustering_theoretical_support <- function(N){
  CC <- 0L
  i <- 1L
  for(ki in 2:(N-1)){
    for(triangles in 0:(ki*(ki-1))){   #ki*ki-1 is the maximum number of triangles attached to the node i
      
      CC[i] <- triangles/(ki*(ki-1))
      i <- i+1
    }
  }
  return( sort(unique(CC)) ) 
}

# theSup <- c(NaN, clustering_theoretical_support(N))
# length(theSup)
# length(sim.cl.values)
# round(sim.cl.values, 9) %in% round(c(NaN,theSup),9) %>% sum/length(sim.cl.values)

# Fuzzy Connectivity ------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Given an undirected network with N nodes, this function compute NUMERICALLY for each node, the probability that the network is connected with a certain number of edges. 
# This function takes the sampled adjacency matrices and computes the number of components for each node of each sampled network. 
# Then, it compute the frequency of the networks actually connected.
# NB at the end, the probabilities does not sum (necessarily) to one! (It does if also the other components are considered)
#
# INPUT  --> Asample = List of sampled adjacency matrices (N x N each)
# OUTPUT --> simPcc = Dataframe with probabilities of the network of being connected with the corresponding number of edges

get_fuzzy_Connectivity_numerical <- function(Asample, doPlot=F){
  
  dp <- # number n.CC of components connected with n.edges number of edges
    lapply( Asample , function(x){
      g <- graph_from_adjacency_matrix(x, mode = 'undirected')
      n.edges <- length(E(g))
      n.CC <- count_components(g) 
      return( cbind(n.edges, n.CC) )
    }) %>% 
    do.call(rbind,.) %>% 
    as.data.frame 
  
  simPcc <- 
    as.data.frame(table(dp)/length(Asample)) %>% 
    dplyr::filter(n.CC==1) %>% 
    dplyr::select(n.edges, Freq)
  
  if (dim(simPcc)[1] > 0) {
    simPcc$n.edges <- as.numeric(levels(simPcc$n.edges))
    return(simPcc)
  } else{
    return(data.frame(n.edges = NA, Freq = NA))
  }
  
  # if( length(simPcc)>2 & doPlot){
  #   plot( as.numeric(names(simPcc)), simPcc , type='b', pch=1, col='firebrick', xlim=c(0,m), main='Probability of existence of the Giant Component\n(Numerical)', sub = paste('N =',N), xlab = 'Number of edges', ylab='Probability')#,xaxt='n')
  #   legend('topright' , legend=c('Analytical','Numerical'), col=c('darkblue','firebrick'), lty=1, cex=0.6 , bty ='n' )
  # }
  
}




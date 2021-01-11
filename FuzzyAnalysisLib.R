# Fuzzy Degree ------------------------------------------------------------------------------------------------------------------------------------------------------------------

# This function compute the fuzzy degree for the nodes of a network
# INPUT --> FuzzyProb = matrix of probabilities of existence
#       --> analytic = T if the fuzzy degree must be from Poisson-Distribution rather than sampling FuzzyProb
#       --> Nsam --> numeber of random sampling from the FuzzyProb (if analytic==F)
# OUTPUT --> Pdnode = Dataframe with probabilities of a node having a certain degree (degrees on rows and nodes on columns)

library(poibin)
source('./Functions/randomat.R')

get_fuzzy_Degree <- function(P, analytic=T, Wsample=NULL){
  
  if(analytic) {
    N <- dim(P)[1]
  } else{
    N <- dim(Wsample[[1]])[1]
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
    simuls.node <- lapply(Wsample, function(x) rowSums(x)) # sum the "successes" for each node (degree)
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

# Given a network with N nodes, this function compute NUMERICALLY for each node, the probability that the (local) clustering coefficient is equal to a certain value. 
# This function takes the sampled adjacency matrices and computes the transitivity for each node of each sampled network. 
# Then, it compute the frequency of the transitivity values for every node.

# INPUT: 
# Wsample = list of sampled adjacency matrix (NxN each)

library(igraph)
get_fuzzy_Clustering_numerical <- function(Wsample){
  
  N <- dim(Wsample[[1]])[1]
  Nsim <- length(Wsample)
  sim.cl <- lapply( Wsample , function(x) transitivity( graph_from_adjacency_matrix( x, mode = 'undirected'), 'local') )
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

# Given a network with N nodes, this function compute NUMERICALLY for each node, the probability that the network is connected with a certain number of edges. 
# This function takes the sampled adjacency matrices and computes the number of components for each node of each sampled network. 
# Then, it compute the frequency of the networks actually connected.
# NB at the end, the probabilities does not sum (necessarely) to one! (It does if also the other components are considered)

get_fuzzy_Connectivity_numerical <- function(Wsample, doPlot=T){
  
  N <- dim(Wsample[[1]])[1]
  Nsim <- length(Wsample)
  m <- choose(N,2)
  sim_nedges <- lapply(Wsample, function(x) x[upper.tri(x)] %>% sum ) %>% unlist # number o edges for each sample
  sim_nCC <- lapply( Wsample , function(x) count_components( graph_from_adjacency_matrix(x, mode = 'undirected')) ) %>% unlist # number of components connected with sim_nedges number of edges
  
  dp <- data.frame(sim_nedges, sim_nCC)
  simPccTot <- as.matrix(table(dp)/Nsim)
  simPcc <- simPccTot[,1]
  simPcc <- cbind('nedges'=as.numeric(names(simPcc)),'probability'=simPcc)
  
  # if( length(simPcc)>2 & doPlot){
  #   plot( as.numeric(names(simPcc)), simPcc , type='b', pch=1, col='firebrick', xlim=c(0,m), main='Probability of existence of the Giant Component\n(Numerical)', sub = paste('N =',N), xlab = 'Number of edges', ylab='Probability')#,xaxt='n')
  #   legend('topright' , legend=c('Analytical','Numerical'), col=c('darkblue','firebrick'), lty=1, cex=0.6 , bty ='n' )
  # }
  
  return(simPcc)
  
}




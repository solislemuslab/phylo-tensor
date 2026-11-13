library(SiPhyNetwork)

##Simulation parameters
set.seed(47176870) # random seed. BB(5) 
simulation_dir <-"./sims/nets/"
n_phy <- 10
ntips <- 15
lambda <- 0.9 #speciation
mu <- 0.0 #no extinction
nu <- 0.01 #hybridization. Set to 0 for trees
hybprops <- c(0.0, 0.0, 1.0) #only lineage neutral  
hyb.inher.fxn <- make.beta.draw(10, 10) #generally inheritance proportions around 0.5
scalar<- 1 ##scale the branch lengths by this amount. Higher means less ILS. Closer to 0 means more ILS


##A function to define what a 'valid' simulation is
##You can write any conditions here or comment out any irrelevant conditions
is_valid <- function(net){
  
  #check that the simulation return a phylogeny
  #only really relevant with incomplete sampling or extinction
  if(!('evonet' %in% class(net))){
    return(FALSE)
  }
  
  #check that the simulation generated a network (more than 1 reticulation)
  nrets <- nrow(net$reticulation)
  if(nrets==0){ #could be changed to require exactly one reticulation, etc.
    return(FALSE)
  }
  
  # #Commented out -  Only relevant if we have lineage degenerative hybridization or extinction
  # #check that the root is the Last stable ancestor (LSA).
  # ltt <- ltt.network(net)
  # if(any(ltt$n_lineages == 1)){ #Root not LSA if we ever come back to one lineage
  # 
  # }
  
  #check network level
  nlevel <-getNetworkLevel(net)
  if(nlevel != 1){#reject if not level-1
    return(FALSE)
  }
  
  return(TRUE)
}



dir.create(simulation_dir,showWarnings=F,recursive=T)

n_success <- 0 
while(n_success < n_phy){
  #simulate network
  net <- sim.bdh.taxa.ssa(n = ntips, numbsim = 1, lambda = lambda, mu = mu, nu = nu,
                   hybprops = hybprops, hyb.inher.fxn = hyb.inher.fxn,
                   frac = 1, twolineages = TRUE, complete = TRUE,
                   stochsampling = FALSE, hyb.rate.fxn = NULL, trait.model = NULL)[[1]]
  #print(paste("Reticulations:", nrow(net$reticulation),
            #"Level:", getNetworkLevel(net)))
  if(!is_valid(net)){#Check validity of simulation. See function definition above
    next #retry
  }
  
  #we have a valid simulation
  
  #scale branch lengths by scalar
  net$edge.length<- net$edge.length*scalar
  
  #save result
  n_success<-n_success+1
  net_dir <- paste0(simulation_dir,'net_',n_success,'/')
  dir.create(net_dir,showWarnings=FALSE,recursive = T)
  write.net(net,paste0(net_dir,'net.tre')) #save the network
  print(n_success)
}


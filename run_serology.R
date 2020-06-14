# Code for simulations in paper: Potential biases arising from epidemic dynamics in observational seroprotection studies
# Code adapted from Matt Hitchings: https://github.com/mhitchings/Code

require(NetSurv)
require(Matrix)
require(Rlab)
require(igraph)
require(deSolve)
require(reshape2)
require(ggplot2)
require(caTools)
require(ggpubr)
library(survival)
library(coxme)
require(tidyverse)
require(splitstackshape)
require(dplyr)
require(coxphf)


# code to output multiple items from function
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}


# Make network, run epidemic -----
make_network <- function(ave_community_size, community_size_range, num_communities, rate_within, rate_between) {
  # Function to make a network of closely-connected communities that are more sparsely
  # connected with each other. I use a stochastic block model.
  # Inputs:
  # Ave_community_size is the average number of people in a community
  # Community_size_range is the range of community sizes. 
  # is being chosen from a uniform distribution on (ave-range/2, ave+range/2)
  # Num_communities is the number of communities in the study population
  # rate_within is the probability of an edge between any two nodes within the same community
  # rate_between is the probability of an edge between any two nodes in "different" communities
  
  require(NetSurv)
  require(Matrix)
  require(Rlab)
  require(igraph)
  require(deSolve)
  require(reshape2)
  require(ggplot2)
  
  # Create the network, and assign all members a community number
  community_sizes <- ave_community_size + round(runif(num_communities,-community_size_range/2,community_size_range/2))
  studypop_size <- sum(community_sizes)
  # Currently all communities have the same connectedness, and all communities are equally
  # connected to each other
  within_rates <- diag(nrow=num_communities,ncol=num_communities,x=rate_within)
  between_rates <- matrix(rate_between,nrow=num_communities,ncol=num_communities) -
    diag(nrow=num_communities,ncol=num_communities,x=rate_between)
  rates<-within_rates+between_rates
  
  g <- sample_sbm(studypop_size,rates,community_sizes)
  # Give the nodes a name so that igraph remembers them
  V(g)$name<-1:studypop_size
  V(g)$community<-rep(1:num_communities,community_sizes)

  return(g)
  
}

network_epidemic<-function(g,beta,SP,control_day,beta_control,num_introductions,num_timesteps,
                           incperiod_shape,incperiod_rate,infperiod_shape,infperiod_rate,sim) {
  # Inputs:
  # g - the graph to run the epidemic on
  # beta - Every infectious individual contacts all their neighbours in a time step
  # and infects each susceptible with hazard beta. So beta represents the hazard of infection from one
  # contact between an infectious and a susceptible.
  # num_introductions - how many separate introductions we expect on average from the main epidemic. This is used
  # to calibrate the external force of infection
  
  
  list <- structure(NA,class="result")
  "[<-.result" <- function(x,...,value) {
    args <- as.list(match.call())
    args <- args[-c(1:2,length(args))]
    length(value) <- length(args)
    for(i in seq(along=args)) {
      a <- args[[i]]
      if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
    }
    x
  }
  
  # Recover and spread functions
  recover<-function(e_nodes,i_nodes,r_nodes,infperiod_shape,infperiod_rate) {
    # Input is a list of the exposed nodes, 
    # with number of days since infection and total incubation/latent
    # period, and equivalently for the infectious nodes.
    # For each of these nodes, we will add it to newinfectious if the number of days exposed has
    # reached the total length of the incubation period, and equivalently for the infectious nodes.
    
    # Advance infectious nodes first, otherwise you will doubly advance any nodes switching from
    # exposed to infectious at this time step
    indices_to_remove <- i_nodes[2,]>=i_nodes[3,]
    newremoved<-as.vector(i_nodes[1,indices_to_remove])
    
    # Add one day to the length of each infected individual's time infected
    i_nodes[2,] <- i_nodes[2,]+1
    
    # Remove any recovered from i_nodes and add them to recovered (note this is the second S in the SEIS')
    i_nodes <- i_nodes[,!(i_nodes[1,] %in% newremoved),drop=FALSE]
    r_nodes <- c(r_nodes,newremoved)
    
    # Now advance exposed nodes
    indices_to_remove <- e_nodes[2,]>=e_nodes[3,]
    newinfectious<-as.vector(e_nodes[1,indices_to_remove])
    
    # Add one day to the length of each infected individual's time infected
    e_nodes[2,] <- e_nodes[2,]+1
    
    # Remove any progressing from e_nodes and add to i_nodes
    e_nodes <- e_nodes[,!(e_nodes[1,] %in% newinfectious),drop=FALSE]
    inf_periods <- rgamma(length(newinfectious),infperiod_shape,infperiod_rate)
    i_nodes <- cbind(i_nodes,rbind(newinfectious,rep(0,length(newinfectious)),inf_periods))
    
    list(e_nodes, i_nodes, r_nodes, sort(newinfectious))
  }
  
  spread<-function(g, s_nodes, e_nodes, i_nodes,r_nodes,t,Susceptible,
                   beta,SP,control_day,beta_control,
                   incperiod_shape, incperiod_rate,
                   connected_nodes,external_inf_F,source_num_inf){
    
    # Spread will create new infected nodes from two sources: infectious nodes within the the study
    # population, and external pressure from the source population
    # Inputs:
    # g is the graph, used to find neighbours of infected nodes
    # s_nodes, e_nodes, i_nodes, r_nodes are susceptible, exposed, infected and susceptible' nodes
    # beta is the hazard of infection for one contact
    # incperiod_shape and rate are used to assign each newly exposed node a latent/incubation period
    # length, currently drawn from a gamma distribution
    # connected_nodes is a list of nodes that are connected the the source population
    # external_inf_F is a constant of proportionality that defines infectious pressure from source population 
    # to an individual
    # source_num_inf is the number of infectious individuals in the source population
    
    # Process: go through list of i_nodes, and choose a random number of its susceptible neighbours to
    # be infected, according to beta and choose a random number of its susceptible seroprotected neighbours to
    # be infected, according to beta and SP (seroprotection)
    # Then go through list of nodes that are connected to source population and infect each susceptible
    # one with probability 1-exp(-FI), where I is the number/proportion of infectious, and F is a constant
    
    if (ncol(i_nodes)>0) {
      
      # note if control_day has occurred; if so reduce beta
      if (t >= control_day){
        beta <- beta*beta_control
      }
      #print(beta)
      
      # Make a beta vector
      betavec <- rep(beta,ncol(i_nodes))
      
      # Get a list of all susceptible neighbours of all infected nodes
      potential_contacts<-lapply(i_nodes[1,],function(x) neighbors(g,x))
      susc_contacts<-lapply(potential_contacts,function(x,susceptibles) intersect(x,susceptibles),susceptibles=s_nodes)
      num_neighbours_susc<-rapply(susc_contacts,length)
      # Sample from each group of neighbours in turn
      # First choose how many neighbours each node infects
      num_contacts_susc<-rbinom(length(num_neighbours_susc),num_neighbours_susc,1-exp(-betavec))
      # Then sample from the neighbours
      # If one node gets picked twice by "different" nodes, just discard the duplicate.
      # In the very rare case that each i_nodes makes a number of new infectees equal to the number
      # of infectious nodes, mapply will make a matrix and unlist won't work. Therefore, the c() around
      # it ensures we turn the matrix into a vector. Unique then removes duplicates.
      infectees_s<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=susc_contacts,y=num_contacts_susc)))
      infectees_s<-unique(infectees_s)
      
      # Make a beta vector
      betavec_R <- rep(beta*SP,ncol(i_nodes)) # if SP is 1, beta is same; if 0, Rs cannot get infected
      
      # Get a list of all recovered neighbours of all infected nodes
      
      if  (length(r_nodes)>0){
        # Get a list of all neighbours of all infected nodes
        potential_contacts<-lapply(i_nodes[1,],function(x) neighbors(g,x))
        r_contacts<-lapply(potential_contacts,function(x,susceptibles) intersect(x,susceptibles),susceptibles=r_nodes)
        num_neighbours_r<-rapply(r_contacts,length)
        # Sample from each group of neighbours in turn
        # First choose how many neighbours each node infects
        num_contacts_r<-rbinom(length(num_neighbours_r),num_neighbours_r,1-exp(-betavec_R))
        # Then sample from the neighbours
        # If one node gets picked twice by "different" nodes, just discard the duplicate.
        # In the very rare case that each i_nodes makes a number of new infectees equal to the number
        # of infectious nodes, mapply will make a matrix and unlist won't work. Therefore, the c() around
        # it ensures we turn the matrix into a vector. Unique then removes duplicates.
        infectees_r<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=r_contacts,y=num_contacts_r)))
        infectees_r<-unique(infectees_r)
      } else{
        infectees_r <- c()
      }
      infectees <- unique(c(infectees_s,infectees_r))
      
    } else {
      infectees <- c()
    }
    # Pick out the nodes connected to the source that are still susceptible 
    # and haven't just been infected
    target_snodes <- setdiff(intersect(connected_nodes,s_nodes),infectees) 
    if (length(r_nodes)>0){
      target_rnodes <- sample(r_nodes,length(r_nodes)*SP,replace=FALSE)
    } else{
      target_rnodes <- c()
    }
    target_cnodes <- c(target_snodes,target_rnodes)
    
    # Make a vector to represent external infection hazard for each individual
    communities_s <- V(g)[target_cnodes]$community
    comm_sizes_s <- sapply(1:num_communities,function(x) sum(communities_s==x))
    
    # Hazard of infection
    extFs_s<-rep(extF,comm_sizes_s)
    
    # Probability of infection
    prob_inf_fromsource <- 1 - exp(-mean(extFs_s)*source_num_inf)
    
    # Choose a number of individuals to be infected, then sample those individuals
    if (length(target_cnodes)>0) {
      num_conn_inf <- rbinom(1,length(target_cnodes),prob_inf_fromsource)
      conn_inf <- target_cnodes[sample.int(length(target_cnodes),num_conn_inf,prob=extFs_s)]
    } else {
      conn_inf <- c()
    }
    
    
    newinfected <- c(infectees,conn_inf)
    newinfected <- unique(newinfected)
    
    if (length(newinfected)>0) {
      
      # Give each newly exposed node an incubation/latent period
      inc_periods <- rgamma(length(newinfected),incperiod_shape,incperiod_rate)
      # Add them to e_nodes and remove from s_nodes
      e_nodes <- cbind(e_nodes,rbind(newinfected,rep(0,length(newinfected)),inc_periods))
      s_nodes<-setdiff(s_nodes,newinfected)
      
    }
    
    Susceptible <- rbind(Susceptible,cbind(t,length(s_nodes)))
    
    list(s_nodes, e_nodes,newinfected,Susceptible)
  }
  
  
  #### RUN THE EPIDEMIC IN THE SOURCE POPULATION ####
  # This is to define external infectious pressure to the network
  
  model <- function(t, y, parms) {
    with(as.list(c(y,parms)), {
      beta <- betahat * exp(-a1 * (t - atau))
      dS <- -beta * S * I / (S+E+I+R)
      dE <- beta * S * I/  (S+E+I+R) - sigma *E
      dI <- sigma * E - gamma * I
      dR <- gamma * I
      list(c(dS,dE,dI,dR))
    })
  }
  
  # number of people in source population
  NS <- 1000000
  # starting values
  y<- c(S=NS-1,E=0,I=1,R=0)
  times<-seq(0,50,1)
  parms<-c(betahat=0.95,a1=0.04,a2=0.2,atau=20,sigma=0.2,gamma=0.2)
  out<-as.data.frame(lsoda(y,times,model,parms))
  out[51:num_timesteps,4] <- 0 #only allow introductions to occur until day 50
  #plot(out$I)

    #### RUN THE EPIDEMIC IN THE STUDY POPULATION ####
  
  # Define how the study population is linked to the source population
  # Connect all individuals to source population at same hazard
  # Constant of proportionality varies by community
  studypop_size<-length(V(g))
  connected_to_source <- V(g)$name
  
  # Calibrate extF to the number of introductions, given the progression of the epidemic in the source population
  num_communities <- max(V(g)$community)
  comm_sizes <- sapply(1:num_communities,function(x) length(V(g)[community==x]))
  sumsqrt <- sum(sqrt(comm_sizes))
  extF <- -log(1-num_introductions/(sqrt(comm_sizes)*sumsqrt))/trapz(times,out$I)
  
  # Initialize the S, E, I, and R nodes. I seed the epidemic from an SEIR curve in a source population,
  # so initially all nodes in the study population are susceptible
  # i_nodes and e_nodes are matrices. The first row is the  identity of the node. The second row
  # is the number of days since infection/infectiousness. The third row is the total incubation/infectious period, 
  # drawn from a distribution when it becomes infected/infectious.
  e_nodes<-matrix(nrow=3,ncol=0)
  i_nodes<-matrix(nrow=3,ncol=0)
  s_nodes <- as.vector(V(g))
  r_nodes <- c()
  Susceptible <- NULL
  
  
  # Initialize results.
  results<-data.frame("SimulationNumber"=rep(NA,studypop_size*12),
                      "InfectedNode"=rep(NA,studypop_size*12),
                      "DayInfected"=rep(NA,studypop_size*12),
                      "Community"=rep(NA,studypop_size*12),
                      "Seroconvert"=rep(NA,studypop_size*12))
  numinfectious<-0
  
  for (t in 1:num_timesteps) {
    cat(t)
    # I'm recovering first, so I need to ensure that everyone has at least one chance to infect.
    # I do this by initializing an infectious node with 0 days since infection, seeing whether they
    # recover, then advancing them one day along their infectious period.
    
    # Only need to recover if there are any infected or exposed
    if ((ncol(i_nodes)>0)||(ncol(e_nodes)>0)) {
      list[e_nodes,i_nodes,r_nodes,newinfectious]<-
        recover(e_nodes,i_nodes,r_nodes,infperiod_shape,infperiod_rate)
      
    } else {
      newinfectious <- c()
    }
    
    if (ncol(e_nodes)==0) {
      newinfected <- 0
    }
    
    Susceptible <- as.data.frame(Susceptible)
    
    list[s_nodes,e_nodes,newinfected,Susceptible]<-
      spread(g,s_nodes,e_nodes,i_nodes,r_nodes,t,Susceptible,
             beta,SP,control_day,beta_control,
             incperiod_shape,incperiod_rate,
             connected_to_source,extF,out$I[t])
    
    
    numnewinfectious<-length(newinfectious)
    if (numnewinfectious>0) {
      
      newcommunities <- V(g)[newinfectious]$community
      
      # Update results
      results$SimulationNumber[(numinfectious+1):(numinfectious+numnewinfectious)]<-rep(sim,numnewinfectious)
      results$InfectedNode[(numinfectious+1):(numinfectious+numnewinfectious)]<-newinfectious
      results$DayInfected[(numinfectious+1):(numinfectious+numnewinfectious)]<-rep(t,numnewinfectious)
      results$Community[(numinfectious+1):(numinfectious+numnewinfectious)]<-newcommunities
      # seroconversion occurs 7 days after become infectious
      results$Seroconvert[(numinfectious+1):(numinfectious+numnewinfectious)]<-rep((t+7),numnewinfectious)

      numinfectious <- numinfectious+numnewinfectious
      
    }
  }
  
  # Tidy up results
  if (numinfectious>0) {
    results<-results[1:numinfectious,]
  } else {
    results<-results[1,]
    results$SimulationNumber[1]<-sim
  }
  
  # track susceptibles remaining
  Susceptible$newexposed <- lag(Susceptible$V2) - Susceptible$V2
  Susceptible$newexposed[is.na(Susceptible$newexposed)] <- 0
  Susceptible$hazard <- Susceptible$newexposed/Susceptible$V2
  
  list(results,Susceptible)
  
}


analyze <- function(results,start_followup1,start_followup2,start_followup3,num_timesteps,match_analysis) {
  
  # for people infected multiple times, change Seroconversion date to that of first time infected
  results_agg <- aggregate(results$SimulationNumber,by=list(results$InfectedNode),length)
  results_agg_mult <- results_agg[results_agg$x>1,]
  if (nrow(results_agg_mult)>0){
    for (n in 1:nrow(results_agg_mult)){
      results$Seroconvert[results$InfectedNode==results_agg_mult$Group.1[n]] <- 
        min(results$Seroconvert[results$InfectedNode==results_agg_mult$Group.1[n]])
    }
  }
  
  # make results list of those infected
  results_analysis<-results[results$DayInfected<=num_timesteps,]
  
  # Assign them eventstatus=1 for the Cox analysis
  results_analysis$eventstatus<-rep(1,nrow(results_analysis))
  
  # Get a list of nodes that were enrolled in the trial but never infected
  noninf<-setdiff(V(g)$name,results$InfectedNode)
  
  # Make data frame for those who were never infected (i.e. censored by end of study)
  noninfdf<-data.frame(InfectedNode=noninf,
                       DayInfected=rep(num_timesteps,length(noninf)),
                       Community=V(g)$community[V(g)$name %in% noninf],
                       Seroconvert=rep(NA,length(noninf)),
                       eventstatus=rep(0,length(noninf)))
  
  # Remove column with simulation number so the columns match up
  results_analysis$SimulationNumber<-NULL
  results_analysis<-rbind(results_analysis,noninfdf)
  
  # Randomly assign day of enrollment and serostatus relative to day of enrollment
  Day_enrollment <- cbind("InfectedNode" = V(g)$name, "Day_enrollment" = sample(c(start_followup1,start_followup2,start_followup3),length(V(g)),replace=TRUE))
  results_analysis <- merge(Day_enrollment,results_analysis,by="InfectedNode")
  results_analysis$serostatus <- ifelse(results_analysis$Seroconvert<=results_analysis$Day_enrollment,1,0)
  results_analysis$serostatus[is.na(results_analysis$serostatus)] <- 0
  
  # Remove those that were infected before day of enrollment from results_analysis
  results_analysis_pre <- results_analysis[results_analysis$DayInfected <= results_analysis$Day_enrollment,]
  results_analysis_pre <- results_analysis_pre[!is.na(results_analysis_pre$InfectedNode),]
  results_analysis <- results_analysis[results_analysis$DayInfected > results_analysis$Day_enrollment,]
  results_analysis <- results_analysis[!is.na(results_analysis$InfectedNode),]
  
  # make event status 0 and day infected end of trial for those infected before enrollment
  results_analysis_pre$eventstatus <- 0
  results_analysis_pre$DayInfected <- num_timesteps
  # add back in those from pre that weren't reinfected
  infected_once <- setdiff(results_analysis_pre$InfectedNode,results_analysis$InfectedNode)
  results_analysis <- rbind(results_analysis,unique(results_analysis_pre[results_analysis_pre$InfectedNode %in% infected_once,]))
  
  # make day infected relative to day of enrollment
  results_analysis$DayInfected <- results_analysis$DayInfected - results_analysis$Day_enrollment
  # for those infected multiple times after enrollment, just keep the first time
  results_analysis <- results_analysis[order(results_analysis$DayInfected),]
  results_analysis <- distinct(results_analysis, InfectedNode, .keep_all= TRUE)
  
  # if analysis is matched, match on day of enrollment and community
  if (match_analysis==1){
    results_analysis_sample <- data.frame()
    communities <- unique(results$Community)
    for (com in 1:length(communities)){
      results_analysis_comm <- results_analysis[results_analysis$Community==communities[com],]
      results_analysis_1_pos <- results_analysis_comm[results_analysis_comm$Day_enrollment==start_followup1 & 
                                                        results_analysis_comm$serostatus==1,]
      results_analysis_2_pos <- results_analysis_comm[results_analysis_comm$Day_enrollment==start_followup2 & 
                                                        results_analysis_comm$serostatus==1,]
      results_analysis_3_pos <- results_analysis_comm[results_analysis_comm$Day_enrollment==start_followup3 & 
                                                        results_analysis_comm$serostatus==1,]
      results_analysis_1_neg <- results_analysis_comm[results_analysis_comm$Day_enrollment==start_followup1 & 
                                                        results_analysis_comm$serostatus==0,]
      results_analysis_2_neg <- results_analysis_comm[results_analysis_comm$Day_enrollment==start_followup2 & 
                                                        results_analysis_comm$serostatus==0,]
      results_analysis_3_neg <- results_analysis_comm[results_analysis_comm$Day_enrollment==start_followup3 & 
                                                        results_analysis_comm$serostatus==0,]
      
      # match 1:1 (so need to see which group is larger - seropositive or seronegative)
      if(nrow(results_analysis_1_pos) > nrow(results_analysis_1_neg)){
        larger_1 <- results_analysis_1_pos
        smaller_1 <- results_analysis_1_neg
      } else{
        smaller_1 <- results_analysis_1_pos
        larger_1 <- results_analysis_1_neg
      }
      
      if(nrow(results_analysis_2_pos) > nrow(results_analysis_2_neg)){
        larger_2 <- results_analysis_2_pos
        smaller_2 <- results_analysis_2_neg
      } else{
        smaller_2 <- results_analysis_2_pos
        larger_2 <- results_analysis_2_neg
      }
      
      if(nrow(results_analysis_3_pos) > nrow(results_analysis_3_neg)){
        larger_3 <- results_analysis_3_pos
        smaller_3 <- results_analysis_3_neg
      } else{
        smaller_3 <- results_analysis_3_pos
        larger_3 <- results_analysis_3_neg
      }
      
      frac_1 <- nrow(smaller_1)/nrow(larger_1)
      frac_2 <- nrow(smaller_2)/nrow(larger_2)
      frac_3 <- nrow(smaller_3)/nrow(larger_3)
      
      results_analysis_sample <- rbind(results_analysis_sample,
                                       smaller_1,
                                       dplyr::sample_frac(larger_1,frac_1,replace=FALSE),
                                       smaller_2,
                                       dplyr::sample_frac(larger_2,frac_2,replace=FALSE),
                                       smaller_3,
                                       dplyr::sample_frac(larger_3,frac_3,replace=FALSE)
      )
    } 
  } else{
    # if analysis not matched take a random sample 50% 
    results_analysis_sample <- dplyr::sample_frac(results_analysis,0.5,replace=FALSE)
  }
  
  # do an unstratified cox model
  survmodel<-try(coxph(Surv(DayInfected,eventstatus)~serostatus,results_analysis_sample),silent=T)
  
  coxph(Surv(Day_enrollment,Day_enrollment+DayInfected,eventstatus)~serostatus,results_analysis_sample)
  # PH test
  #PH <- cox.zph(survmodel)
  #PH_pval <- PH$table[3]
  usesurvmod <- !inherits(survmodel, 'try-error')
  
  if (usesurvmod && vcov(survmodel)>=0){
    # If no error was thrown and the variance is positive, use the results of the model
    
    serEffEst <- exp(survmodel$coefficient + c(0, -1.96, 1.96)*sqrt(survmodel$var))
    zval <- survmodel$coefficient/sqrt(survmodel$var)
    
  } else {
    
    serEffEst<-c(NA,NA,NA)
    
  }
  
  # analysis stratified by Community and day enrollment
  survmodel<-try(coxph(Surv(DayInfected,eventstatus)~serostatus+strata(Day_enrollment)+strata(Community),results_analysis_sample),silent=T)
  
  #PH_strat <- cox.zph(survmodel)
  #PH_pval_strat <- PH_strat$table[3]
  usesurvmod <- !inherits(survmodel, 'try-error')
  
  if (usesurvmod && vcov(survmodel)>=0){
    # If no error was thrown and the variance is positive, use the results of the model
    
    serEffEst_strat <- exp(survmodel$coefficient + c(0, -1.96, 1.96)*sqrt(survmodel$var))
    zval_strat <- survmodel$coefficient/sqrt(survmodel$var)
    
  } else {
    
    serEffEst_strat<-c(NA,NA,NA)
    
  }
  
  # do an unstratified cox model with left truncation
  survmodel<-try(coxph(Surv(Day_enrollment,Day_enrollment+DayInfected,eventstatus)~serostatus,results_analysis_sample),silent=T)
  
  # PH test
  #PH <- cox.zph(survmodel)
  #PH_pval <- PH$table[3]
  usesurvmod <- !inherits(survmodel, 'try-error')
  
  if (usesurvmod && vcov(survmodel)>=0){
    # If no error was thrown and the variance is positive, use the results of the model
    
    serEffEst_LT <- exp(survmodel$coefficient + c(0, -1.96, 1.96)*sqrt(survmodel$var))
    zval <- survmodel$coefficient/sqrt(survmodel$var)
    
  } else {
    
    serEffEst_LT<-c(NA,NA,NA)
    
  }
  
  # analysis with left truncation stratified by Community and day enrollment
  survmodel<-try(coxph(Surv(Day_enrollment,Day_enrollment+DayInfected,eventstatus)~serostatus+strata(Day_enrollment)+strata(Community),results_analysis_sample),silent=T)
  
  #PH_strat <- cox.zph(survmodel)
  #PH_pval_strat <- PH_strat$table[3]
  usesurvmod <- !inherits(survmodel, 'try-error')
  
  if (usesurvmod && vcov(survmodel)>=0){
    # If no error was thrown and the variance is positive, use the results of the model
    
    serEffEst_strat_LT <- exp(survmodel$coefficient + c(0, -1.96, 1.96)*sqrt(survmodel$var))
    zval_strat <- survmodel$coefficient/sqrt(survmodel$var)
    
  } else {
    
    serEffEst_strat_LT<-c(NA,NA,NA)
    
  }
  
  # if only one arm has a case, set HR to either 0 or 1
  
  results_analysis_sample %>%
    group_by(serostatus,eventstatus) %>%
    filter(eventstatus==1)%>%
    summarise(n=n())-> results_analysis_sample_summary
  
  serostatus_sample <- unique(results_analysis_sample_summary$serostatus)
  if (length(serostatus_sample)==1){
    if(serostatus_sample==0){
      serEffEst <- c(0,0,0)
      serEffEst_strat <- c(0,0,0)
      serEffEst_LT <- c(0,0,0)
      serEffEst_strat_LT <- c(0,0,0)
    } else if(serostatus_sample==1){
      serEffEst <- c(Inf,Inf,Inf)
      serEffEst_strat <- c(Inf,Inf,Inf)
      serEffEst_LT <- c(Inf,Inf,Inf)
      serEffEst_strat_LT <- c(Inf,Inf,Inf)
    }
  }
  
  
  analysis_results <- cbind(serEffEst[1],serEffEst[2],serEffEst[3],
                            serEffEst_strat[1],serEffEst_strat[2],serEffEst_strat[3],
                            serEffEst_LT[1],serEffEst_LT[2],serEffEst_LT[3],
                            serEffEst_strat_LT[1],serEffEst_strat_LT[2],serEffEst_strat_LT[3],sim,nrow(results))
  
  list(results_analysis,analysis_results)
  
}


#Constant parameters
nsim<-10
community_size_range<-0
rate_between<-0
num_introductions<-10 

# Gamma-distribution parameters of incubation and infectious period
incperiod_shape<-5
incperiod_rate<-0.9
infperiod_shape<-1.13
infperiod_rate<-0.226
ave_inc_period <- ceiling(incperiod_shape/incperiod_rate)

# length of run
num_timesteps <- 200
beta_control <- 0.5

start_followup1 <- 50
start_followup2 <- 100
start_followup3 <- 150


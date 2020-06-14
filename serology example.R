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


source('./run_serology.R')

#1: 1 comm null well mixed
ave_community_size <- 10000
num_communities <- 1
SP <- 1
rate_within <- 1
beta <- 0.00002
control_day <- 1000

results_master_1 <- data.frame()
analysis_master_1 <- data.frame()
Susceptible_master <- data.frame()

for (sim in 1:nsim){
  cat(sim)
  
  g<-make_network(ave_community_size, community_size_range, num_communities,rate_within, rate_between)
  list[results,Susceptible]<-
    network_epidemic(g,beta,SP,control_day,beta_control,num_introductions,num_timesteps,
                     incperiod_shape,incperiod_rate,infperiod_shape,infperiod_rate,sim)
  results_master_1 <- rbind(results_master_1,cbind(results,rep(sim,nrow(results))))
  Susceptible_master <- rbind(Susceptible_master,cbind(Susceptible,rep(sim,nrow(Susceptible))))
  
  list[results_analysis,analysis_results] <-
    analyze(results,start_followup2,start_followup2,start_followup2,num_timesteps,0)
  analysis_master_1 <- rbind(analysis_master_1,cbind(analysis_results,"Num_comm"=num_communities,"SP"=SP,"follow-up"="same",match="0","rate_within"=rate_within,"control_day"=control_day))
  list[results_analysis,analysis_results] <-
    analyze(results,start_followup1,start_followup2,start_followup3,num_timesteps,0)
  analysis_master_1 <- rbind(analysis_master_1,cbind(analysis_results,"Num_comm"=num_communities,"SP"=SP,"follow-up"="different",match="0","rate_within"=rate_within,"control_day"=control_day))
  list[results_analysis,analysis_results] <-
    analyze(results,start_followup1,start_followup2,start_followup3,num_timesteps,1)
  analysis_master_1 <- rbind(analysis_master_1,cbind(analysis_results,"Num_comm"=num_communities,"SP"=SP,"follow-up"="different",match="1","rate_within"=rate_within,"control_day"=control_day))
  
  write.csv(Susceptible_master,paste0('/n/holyscratch01/lipsitch_lab/rkahn/serology/',j,"_",num_communities,"_",SP,"_",rate_within,"_",beta,"_",control_day,"_","Susceptibility_master.csv"))
  write.csv(results_master_1,paste0('/n/holyscratch01/lipsitch_lab/rkahn/serology/',j,"_",num_communities,"_",SP,"_",rate_within,"_",beta,"_",control_day,"_","results_master_1.csv"))
  write.csv(analysis_master_1,paste0('/n/holyscratch01/lipsitch_lab/rkahn/serology/',j,"_",num_communities,"_",SP,"_",rate_within,"_",beta,"_",control_day,"_","analysis_master_1.csv"))
}

names(analysis_master_1)[1:14] <- c("serEffEst","serEffEst_low","serEffEst_high",
                                    "serEffEst_strat","serEffEst_low_strat","serEffEst_high_strat",
                                    "serEffEst_LT","serEffEst_low_LT","serEffEst_high_LT",
                                    "serEffEst_strat_LT","serEffEst_low_strat_LT","serEffEst_high_strat_LT","sim","num_events")

names(Susceptible_master) <- c("t","newly_exposed","num_susceptible","sim")


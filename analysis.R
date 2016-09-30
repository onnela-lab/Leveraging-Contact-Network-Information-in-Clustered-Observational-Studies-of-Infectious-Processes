rm(list=ls())
library(batch,lib='/home/pcs18/R/x86_64-unknown-linux-gnu-library/3.1')
packages = c("lmtest", "reshape", "data.table", "nlme", "lme4", "dplyr","splines", "stats")
use <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
use(packages)

nc = 15 # number of covariates

setups_num = cbind(
  rep(0:1,each=2**5,times=2**0),
  rep(0:1,each=2**4,times=2**1),
  rep(0:1,each=2**3,times=2**2),
  rep(0:1,each=2**2,times=2**3),
  rep(0:1,each=2**1,times=2**4),
  rep(0:1,each=2**0,times=2**5)
)
setups=apply(setups_num,1,paste,collapse="")
improvements       <- matrix(NA,2**6,nc)
biases             <- matrix(NA,2**6,nc)
powers             <- matrix(NA,2**6,nc)
coverages          <- matrix(NA,2**6,nc)
se_robs            <- matrix(NA,2**6,nc)
se_emps            <- matrix(NA,2**6,nc)

empirical_betas    <- rep(NA, 2**6)
empirical_beta_sds <- rep(NA, 2**6)

setup_count = 0
for(setup in setups){
  setup_count = setup_count + 1
  
  data_filepath = paste("Data",setup,"/",sep="")
  
  
  expit = function(x) 1 / (1+exp(-x))
  logit = function(x) log(x) - log(1-x)
  
  include_baseline_infected_in_analysis = TRUE # if this is false, the other is irrelevant.
  include_baseline_infected_in_outcomes = TRUE
  
  
  trials = 1000     # number of observational studies conducted.
  clusters = 100    # number of clusters for each trial
  M = 5             # number of epidemics performed on each cluster, just in case
  m=1               # which of the epidemic runs you want to use
  n = 200           # number of observations per cluster
  scenario1 = FALSE # If TRUE, this should slice the epidemic data to reflect Scenario 1.
  se.robust = matrix(NA, trials,nc)
  trial_means = rep(NA, trials)
  by_trt = matrix(0,1,5)
  by_outcome = matrix(0,1,4)
  
  cluster_statistics = array(0,c(2, nc-3, trials, clusters))
  trial_statistics = array(0, c(6, nc-3, trials))
  correlations = array(0,c(2, nc-3, trials))
  mean_difference = rep(NA, trials)
  
  coefficients <- c()
  trt_stats    <- c()
  
  pb <- txtProgressBar(0,trials, style = 3)
  for(i in 1:trials){
    coefficients_filename = paste(data_filepath, "CRT_RESULT_",i,".txt",sep="")
    trt_stats_filename    = paste(data_filepath, "CRT_stats_",    i,".txt",sep="")
    if(file.exists(coefficients_filename))
      coefficients = rbind(coefficients,read.table(coefficients_filename,header=T))
    if(file.exists(trt_stats_filename))
      trt_stats <- rbind(trt_stats, unlist(read.table(trt_stats_filename,header=T)))
    setTxtProgressBar(pb, i)
  };close(pb)
  
  coefficients = coefficients[,-1]
  table_data <- coefficients
  empirical_beta = mean(table_data[,1])
  
  colnames(coefficients) = c(paste("GEE(",1:nc,")_beta",sep=""), paste("GEE(",1:nc,")_sd",sep=""))
  write.table(coefficients, paste(data_filepath,"coefficients.txt",sep=""),sep="\t", row.names = FALSE, col.names = TRUE)
  
  coefficients_summary = function(coefficients_data){
    models = dim(coefficients_data)[2]/2
    individual_betas      = coefficients_data[,1:models]
    individual_se.robusts = coefficients_data[,1:models+models]
    betas = apply(coefficients_data[,1:models],2,mean)
    
    bias = empirical_beta - betas
    se.robust    = apply(coefficients_data[,1:models+models],2,mean)
    se.empirical = apply(coefficients_data[,1:models], 2,sd)
    
    rmse.robust    = sqrt(bias^2 + se.robust   ^2)
    rmse.empirical = sqrt(bias^2 + se.empirical^2)
    
    improvement.robust = round(100*(1 - rmse.robust / rmse.robust[1]),4)
    improvement.empirical = round(100*(1 - rmse.empirical / rmse.empirical[1]),4)
    coverage.robust    = 1-apply((individual_betas + qnorm(.975)*individual_se.robusts < empirical_beta) +
                                   (individual_betas - qnorm(.975)*individual_se.robusts > empirical_beta),2,mean)
    power.robust       = apply((individual_betas + qnorm(.975)*individual_se.robusts < 0) +
                                 (individual_betas - qnorm(.975)*individual_se.robusts > 0),2,mean)
    
    summary = round(cbind(bias, se.robust, se.empirical, rmse.robust, rmse.empirical,
                          improvement.robust, improvement.empirical, coverage.robust, power.robust),3)
    
    rownames(summary) = c("Base", paste("mod_",1:(models-1),sep="")) 
    return(summary)
  }
  
  summary = coefficients_summary(table_data)
  write.table(format(summary, digits=3), file = paste(data_filepath, "table.txt",sep=""), sep=" & ", row.names = FALSE, col.names = FALSE, quote = FALSE,append=FALSE)
  write(paste("empirical_beta:",round(empirical_beta,4)), file = paste(data_filepath, "general_information.txt",sep=""), append=FALSE)

  biases[setup_count,]       <- summary[,"bias"]
  improvements[setup_count,] <- summary[,"improvement.empirical"]
  powers[setup_count,]       <- summary[,"power.robust"]
  coverages[setup_count,]    <- summary[,"coverage.robust"]
  se_robs[setup_count,]      <- summary[,"se.robust"]
  se_emps[setup_count,]      <- summary[,"se.empirical"]
  empirical_betas[setup_count]    <- empirical_beta
}

write.table(100*biases,             file = paste("biases.txt",             sep=""), sep=", ", row.names = FALSE, col.names = FALSE, quote = FALSE,append=FALSE)
write.table(improvements,           file = paste("improvements.txt",       sep=""), sep=", ", row.names = FALSE, col.names = FALSE, quote = FALSE,append=FALSE)
write.table(powers,                 file = paste("powers.txt",             sep=""), sep=", ", row.names = FALSE, col.names = FALSE, quote = FALSE,append=FALSE)
write.table(coverages,              file = paste("coverages.txt",          sep=""), sep=", ", row.names = FALSE, col.names = FALSE, quote = FALSE,append=FALSE)
write.table(100*se_robs,            file = paste("se_robs.txt",            sep=""), sep=", ", row.names = FALSE, col.names = FALSE, quote = FALSE,append=FALSE)
write.table(100*se_emps,            file = paste("se_emps.txt",            sep=""), sep=", ", row.names = FALSE, col.names = FALSE, quote = FALSE,append=FALSE)
write.table(100*empirical_betas,    file = paste("empirical_betas.txt",    sep=""), sep=", ", row.names = FALSE, col.names = FALSE, quote = FALSE,append=FALSE)

setups_num = cbind(
  rep(0:1,each=2**5,times=2**0),
  rep(0:1,each=2**4,times=2**1),
  rep(0:1,each=2**3,times=2**2),
  rep(0:1,each=2**2,times=2**3),
  rep(0:1,each=2**1,times=2**4),
  rep(0:1,each=2**0,times=2**5)
)
setups=apply(setups_num,1,paste,collapse="")

blue=rgb(.3,.4,1)
pdf("Exposure_probabilities.pdf",width=8,height=6)
par(cex=1.5)

for(setup in setups){
  print(setup)
  bins   <- c()
  file_name = paste("Data",setup,"/obs_bins.txt",sep="")
  conn <- file(file_name,open="r")
  linn <-readLines(conn)
  for (i in 1:length(linn)){
    line = as.numeric(strsplit(linn[i], " ")[[1]])
    if(length(line)==408){
      bins = rbind(bins, line)
    }
  }
  close(conn)
  
  params = as.data.frame(read.table(paste("Data",setup,"/obs_parameters.txt",sep="")))
  
  clusters = ncol(bins)/4
  OBS = as.matrix(bins[,1:clusters+0*clusters])
  TRT = as.matrix(bins[,1:clusters+1*clusters])
  
  grid = 10
  cuts = seq(min(OBS), max(OBS), length.out=grid)
  cut_xmeans = (c(0,cuts)/2+c(cuts,0)/2)[2:grid]
  cut_means = rep(NA,grid-1)
  cut_sds = rep(NA, grid-1)
  cut_counts = rep(NA,grid-1)
  cut_fit1 = rep(NA,grid-1)
  cut_fit2 = rep(NA,grid-1)
  cut_true = expit(cut_xmeans*1/2*1)
  for(cut_level in 1:(grid-1)){
    cut_counts[cut_level] = sum(((OBS >= cuts[cut_level])*(OBS < cuts[cut_level+1]))==1)
    cut_means[cut_level] = mean(TRT[((OBS >= cuts[cut_level])*(OBS < cuts[cut_level+1]))==1])
    cut_sds[cut_level] = sd(TRT[((OBS >= cuts[cut_level])*(OBS < cuts[cut_level+1]))==1])
  }
  plot(density(OBS),xlim=range(cuts),
       xlab="Exposure Variable",ylab="Density",lwd=2,col=rgb(.5,.5,.5),main="Distribution of Exposure Variable")
  
  plot(cut_xmeans, cut_means, col=NA, ylim=c(0,1),ylab="Probability", xlab="Exposure Variable",main="Exposure Probabilities")
  lines(cut_xmeans, cut_true, lwd=1.75, pch=16,col="gray",lty="dashed")
  for(i in 1:(grid-1)){lines(c(cut_xmeans[i],cut_xmeans[i]),cut_means[i]+c(-1,1)*cut_sds[i]/sqrt(cut_counts[i]),lwd=2)}
  points(cut_xmeans, cut_means, pch=18,col=blue, cex=1.75)
  points(cut_xmeans, logit(cut_fit1), pch=18, cex=1.75, col="green")
  points(cut_xmeans, logit(cut_fit2), pch=18, cex=1.75, col="red")
  legend("bottomright",legend=c("Fraction Exposed","True Probability"),
         col=c(blue,"gray"),pch=15,pt.cex=1.5,bg="white")
}
dev.off()


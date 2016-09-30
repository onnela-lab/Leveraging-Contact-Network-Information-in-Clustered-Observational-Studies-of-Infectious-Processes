#### CLEAR THE CURRENT SESSION
rm(list=ls())
#### CALL FUNCTIONS FOR PARALLELIZATION
library(batch,lib="/home/pcs18/R/x86_64-unknown-linux-gnu-library/3.1")
library(dplyr)
parseCommandArgs()

#### LIBRARIES AND FUNCTIONS
library(Matrix,lib.loc="/opt/R-3.1.2/lib64/R/library")
library(numDeriv,lib.loc="/home/pcs18/R/x86_64-unknown-linux-gnu-library/3.1")
library(geeDoublyRobust,lib.loc ="./lib")
library(OBSgeeDR,lib.loc ="./lib")
library(mgcv)
library(splines)

expit = function(x) 1 / (1+exp(-x))
logit = function(x) log(x) - log(1-x)
working_correlation_structure<-"exchangeable"

setups=apply(cbind(
  rep(0:1,each=2**5,times=2**0),
  rep(0:1,each=2**4,times=2**1),
  rep(0:1,each=2**3,times=2**2),
  rep(0:1,each=2**2,times=2**3),
  rep(0:1,each=2**1,times=2**4),
  rep(0:1,each=2**0,times=2**5)
),1,paste,collapse="")

trials = 1000
if((seed-1) %/% trials == 0){setups = setups[1:32]}
if((seed-1) %/% trials == 1){setups = setups[33:64]}
if(seed > trials) seed = seed - trials

for(setup in setups){
  nameSIM=paste("./Data",setup,"/CRT",sep="")
  trial<-seed
  print(trial)
  
  system(paste("python data_code.py", seed, setup))
  
  file_name = paste(nameSIM,"_",trial,".txt",sep="")
  while(!file.exists(file_name)){Sys.sleep(2)}
  
  data <- read.table(file_name,header=TRUE)
  print(head(data))
  mod1 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs")))
  
  mod.temp <- glm(Outcome~Degree, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(mod.temp,newdata=data,type="response")
  mod.temp <- glm(Outcome~Degree, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(mod.temp,newdata=data,type="response")
  mod2 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  mod.temp <- glm(Outcome~Mean_Neighbor_Degree, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(mod.temp,newdata=data,type="response")
  mod.temp <- glm(Outcome~Mean_Neighbor_Degree, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(mod.temp,newdata=data,type="response")
  mod3 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  mod.temp <- glm(Outcome~Assortativity, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(mod.temp,newdata=data,type="response")
  mod.temp <- glm(Outcome~Assortativity, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(mod.temp,newdata=data,type="response")
  mod4 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  mod.temp <- glm(Outcome~Sex_Worker, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(mod.temp,newdata=data,type="response")
  mod.temp <- glm(Outcome~Sex_Worker, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(mod.temp,newdata=data,type="response")
  mod5 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  mod.temp <- glm(Outcome~LCC_Size, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(mod.temp,newdata=data,type="response")
  mod.temp <- glm(Outcome~LCC_Size, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(mod.temp,newdata=data,type="response")
  mod6 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  mod.temp <- glm(Outcome~Mean_Component_Size, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(mod.temp,newdata=data,type="response")
  mod.temp <- glm(Outcome~Mean_Component_Size, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(mod.temp,newdata=data,type="response")
  mod7 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  mod.temp <- glm(Outcome~Number_Of_Components, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(mod.temp,newdata=data,type="response")
  mod.temp <- glm(Outcome~Number_Of_Components, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(mod.temp,newdata=data,type="response")
  mod8 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  mod.temp <- glm(Outcome~Node_Component_Size, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(mod.temp,newdata=data,type="response")
  mod.temp <- glm(Outcome~Node_Component_Size, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(mod.temp,newdata=data,type="response")
  mod9 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  mod.temp <- glm(Outcome~Total_Neighbor_Seeds, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(mod.temp,newdata=data,type="response")
  mod.temp <- glm(Outcome~Total_Neighbor_Seeds, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(mod.temp,newdata=data,type="response")
  mod10 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  mod.temp <- glm(Outcome~Total_Cluster_Seeds, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(mod.temp,newdata=data,type="response")
  mod.temp <- glm(Outcome~Total_Cluster_Seeds, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(mod.temp,newdata=data,type="response")
  mod11 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  mod.temp <- glm(Outcome~Mins, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(mod.temp,newdata=data,type="response")
  mod.temp <- glm(Outcome~Mins, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(mod.temp,newdata=data,type="response")
  mod12 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  mod.temp <- glm(Outcome~Sums, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(mod.temp,newdata=data,type="response")
  mod.temp <- glm(Outcome~Sums, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(mod.temp,newdata=data,type="response")
  mod13 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  mod.temp <- glm(Outcome~Degree+Mean_Neighbor_Degree+Assortativity+LCC_Size+Mean_Component_Size+Number_Of_Components+Node_Component_Size+Total_Neighbor_Seeds+Total_Cluster_Seeds+Mins+Sums, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(mod.temp,newdata=data,type="response")
  mod.temp <- glm(Outcome~Degree+Mean_Neighbor_Degree+Assortativity+LCC_Size+Mean_Component_Size+Number_Of_Components+Node_Component_Size+Total_Neighbor_Seeds+Total_Cluster_Seeds+Mins+Sums, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(mod.temp,newdata=data,type="response")
  mod14 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  mod.temp <- glm(Outcome~Degree+Mean_Neighbor_Degree+Assortativity+LCC_Size+Mean_Component_Size+Number_Of_Components+Node_Component_Size+Total_Neighbor_Seeds+Total_Cluster_Seeds+Mins+Sums, data=data[which(data[,"Trt"]==1),], family=binomial(link = "logit"));data[,"newvar_trt"]  <- predict(step(mod.temp),newdata=data,type="response")
  mod.temp <- glm(Outcome~Degree+Mean_Neighbor_Degree+Assortativity+LCC_Size+Mean_Component_Size+Number_Of_Components+Node_Component_Size+Total_Neighbor_Seeds+Total_Cluster_Seeds+Mins+Sums, data=data[which(data[,"Trt"]==0),], family=binomial(link = "logit"));data[,"newvar_ctrl"] <- predict(step(mod.temp),newdata=data,type="response")
  mod15 <- try(geeDREstimationObs(formula<-as.formula("Outcome~Trt"),id="Cluster",data=data,family = gaussian,
                         corstr = working_correlation_structure,sandwich.nuisance=FALSE, fay.adjustment=FALSE, ics=FALSE, model.trt.clusterlevel=TRUE,stepwise.trt=FALSE,stepwise.augmentation=TRUE, print.log=FALSE,
                         nameTRT="Trt",nameMISS="MISSING", nameY="Outcome", model.trt=as.formula("Trt~Obs"),  aug=c(ctrl="newvar_ctrl",trt="newvar_trt")))
  
  result<-c(trial,mod1$beta[2],mod2$beta[2],mod3$beta[2],mod4$beta[2],mod5$beta[2],mod6$beta[2],mod7$beta[2],mod8$beta[2],mod9$beta[2],mod10$beta[2],mod11$beta[2],mod12$beta[2],mod13$beta[2],mod14$beta[2],mod15$beta[2],
    sqrt(mod1$var[2,2]),sqrt(mod2$var[2,2]),sqrt(mod3$var[2,2]),sqrt(mod4$var[2,2]),sqrt(mod5$var[2,2]),sqrt(mod6$var[2,2]),sqrt(mod7$var[2,2]),sqrt(mod8$var[2,2]),sqrt(mod9$var[2,2]),sqrt(mod10$var[2,2]),sqrt(mod11$var[2,2]),sqrt(mod12$var[2,2]),sqrt(mod13$var[2,2]),sqrt(mod14$var[2,2]),sqrt(mod15$var[2,2]))#,sqrt(mod2_1$var[2,2]),sqrt(mod2_2$var[2,2])
  print(result)
  result_name = paste(nameSIM,"_RESULT_",trial,".txt",sep="")
  print(result_name)
  write.table(t(result), result_name, row.names = FALSE)  
  
  A1 = mean(unlist((data %>% filter(Trt == 1) %>% group_by(Cluster) %>% summarise(mean(Outcome)))[,2]))
  A0 = mean(unlist((data %>% filter(Trt == 0) %>% group_by(Cluster) %>% summarise(mean(Outcome)))[,2]))
  stats <- unlist(data %>% group_by(Trt) %>% summarise(sums = sum(Outcome), ns = n()))[3:6]
  
  write.table(stats, paste(nameSIM,"_stats_",trial,".txt",sep=""),row.names = FALSE)
  
  TRT = as.data.frame(data %>% group_by(Cluster) %>% summarise(m=max(Trt)))[,2]
  OBS = as.data.frame(data %>% group_by(Cluster) %>% summarise(m=max(Obs)))[,2]
  TOTAL = round(c(OBS, TRT, c(mod1$ps.trt.model$fitted.values, mod2$ps.trt.model$fitted.values)),2)
  write.table(t(TOTAL),paste("./Data",setup,"/obs_bins.txt",sep=""),append=TRUE, row.names=FALSE, col.names=FALSE)
  
  Sys.sleep(1)
}








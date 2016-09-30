#### PRESETS
#install.packages("OBSgeeDR_0.2.tar.gz", repos = NULL, type="source")
#setwd("C:/Users/Patrick/Desktop")
rm(list = ls())

library(xtable)
library(Matrix)
library(mice)
library(OBSgeeDR)
library(dplyr)

data.india<-read.table("complete_data.txt",header=TRUE)
data.india <- data.india[which(data.india$Self_Help != -1),]
data.india$Self_Help<-ifelse(data.india$Self_Help==-1,NA,data.india$Self_Help)

data.india$Sex<-ifelse(data.india$Sex==-1,NA,data.india$Sex)
data.india$Sex<-ifelse(data.india$Sex==2,0,data.india$Sex)
data.india$Age<-ifelse(data.india$Age==-1,NA,data.india$Age)
data.india$Is_Leader<-ifelse(data.india$Is_Leader==-1,0,data.india$Is_Leader)
data.india$leadertake<-data.india$Is_Leader*data.india$Outcome

temp<-aggregate(data.india[,c("leadertake","Self_Help","Sex","Age","Degree","Mean_Neighbor_Degree","Assortativity","LCC_Size","Mean_Component_Size","Node_Component_Size","Total_Neighbor_Seeds","Total_Cluster_Seeds","Mins","Sums","EXP1","EXP2")], list(data.india[,"Village"]), mean,na.rm=TRUE)
names(temp)<-c("Village","leadertake_G","Self_Help_G","Sex_G","Age_G","Degree_G","Mean_Neighbor_Degree_G","Assortativity_G","LCC_Size_G","Mean_Component_Size_G","Node_Component_Size_G","Total_Neighbor_Seeds_G","Total_Cluster_Seeds_G","Mins_G","Sums_G","EXP1_G","EXP2_G")
data.india<-merge(data.india,temp,by="Village")

exp.info<-aggregate(data.india[,c("EXP1","EXP2")], list(data.india[,"Village"]), mean)
pi.a<-data.frame(c(1))
pi.a$EXP1<-sum(exp.info$EXP1)/length(exp.info$EXP1)
pi.a$EXP2<-sum(exp.info$EXP2)/length(exp.info$EXP2)
head(data.india)

########## ESTIMATION

fay_adjustment = TRUE
ICS = TRUE
step_aug = TRUE

dostatistics<-function(data.india,variableEXP,corstr="independence",print.log=TRUE){
  beta0<-matrix(NA,ncol=4,nrow=4)
  beta1<-matrix(NA,ncol=4,nrow=4)
  noadj<-geeDREstimationObs(formula<-as.formula(paste("Outcome~",variableEXP,sep=" ")),
                            id="Village",
                            data=data.india,
                            family = gaussian,
                            pi.a = pi.a[1,variableEXP],
                            corstr = corstr,
                            nameTRT=variableEXP,
                            nameMISS="MISSING",
                            nameY="Outcome",
                            sandwich.nuisance=FALSE,
                            stepwise.trt=TRUE,
                            stepwise.augmentation=step_aug,
                            fay.adjustment=fay_adjustment,
                            ics=ICS,  
                            model.trt=NULL, print.log=print.log)
  
  
  beta0[1,]<-c("unadj",noadj$beta[1],sqrt(noadj$var[1,1]),noadj$beta[1]/sqrt(noadj$var[1,1]))
  beta1[1,]<-c("unadj",noadj$beta[2],sqrt(noadj$var[2,2]),noadj$beta[2]/sqrt(noadj$var[2,2]))
  
  gee<-geeDREstimationObs(formula<-as.formula(paste("Outcome~",variableEXP,sep=" ")),
                          id<-"Village",
                          data<-data.india,
                          family = gaussian,
                          pi.a = pi.a[1,variableEXP],
                          corstr = corstr,
                          nameTRT=variableEXP,
                          nameMISS="MISSING",
                          nameY="Outcome",
                          sandwich.nuisance=FALSE,
                          fay.adjustment=fay_adjustment,
                          ics=ICS,  
                          model.trt=as.formula(paste(variableEXP,"~Self_Help_G+Sex_G+Age_G+Degree_G+ Mean_Neighbor_Degree_G+ Assortativity_G+ LCC_Size_G +Mean_Component_Size_G+ Node_Component_Size_G",sep="")),
                          model.trt.clusterlevel=TRUE,
                          stepwise.trt=TRUE,
                          stepwise.augmentation=step_aug,
                          print.log=print.log)
  
  beta0[2,]<-c("gee",gee$beta[1],sqrt(gee$var[1,1]),gee$beta[1]/sqrt(gee$var[1,1]))
  beta1[2,]<-c("gee",gee$beta[2],sqrt(gee$var[2,2]),gee$beta[2]/sqrt(gee$var[2,2]))
  
  aug<-geeDREstimationObs(formula<-as.formula(paste("Outcome~",variableEXP,sep=" ")),
                          id<-"Village",
                          data<-data.india,
                          family = gaussian,
                          pi.a = pi.a[1,variableEXP],
                          corstr = corstr,
                          model.augmentation.trt =as.formula(Outcome~Degree+Mean_Neighbor_Degree+Assortativity+LCC_Size+Mean_Component_Size+Node_Component_Size+Total_Neighbor_Seeds+Total_Cluster_Seeds+Mins+Sums),
                          model.augmentation.ctrl=as.formula(Outcome~Degree+Mean_Neighbor_Degree+Assortativity+LCC_Size+Mean_Component_Size+Node_Component_Size+Total_Neighbor_Seeds+Total_Cluster_Seeds+Mins+Sums),
                          nameTRT=variableEXP,
                          nameMISS="MISSING",
                          nameY="Outcome",
                          sandwich.nuisance=FALSE,
                          fay.adjustment=fay_adjustment,
                          ics=ICS,  
                          model.trt=as.formula(paste(variableEXP,"~Self_Help_G+Sex_G+Age_G+Degree_G+ Mean_Neighbor_Degree_G+ Assortativity_G+ LCC_Size_G +Mean_Component_Size_G+ Node_Component_Size_G",sep="")),
                          model.trt.clusterlevel=TRUE,
                          stepwise.trt=TRUE,
                          stepwise.augmentation=step_aug,
                          print.log=print.log)
  
  beta0[3,]<-c("aug",aug$beta[1],sqrt(aug$var[1,1]),aug$beta[1]/sqrt(aug$var[1,1]))
  beta1[3,]<-c("aug",aug$beta[2],sqrt(aug$var[2,2]),aug$beta[2]/sqrt(aug$var[2,2]))
  
  augall<-geeDREstimationObs(formula<-as.formula(paste("Outcome~",variableEXP,sep=" ")),
                             id<-"Village",
                             data<-data.india,
                             family = gaussian,
                             pi.a = pi.a[1,variableEXP],
                             corstr = corstr,
                             model.augmentation.trt =as.formula(Outcome~Sex+Age+Degree+Mean_Neighbor_Degree+Assortativity+LCC_Size+Mean_Component_Size+Node_Component_Size+Total_Neighbor_Seeds+Total_Cluster_Seeds+Mins+Sums),
                             model.augmentation.ctrl=as.formula(Outcome~Sex+Age+Degree+Mean_Neighbor_Degree+Assortativity+LCC_Size+Mean_Component_Size+Node_Component_Size+Total_Neighbor_Seeds+Total_Cluster_Seeds+Mins+Sums),
                             nameTRT=variableEXP,
                             nameMISS="MISSING",
                             nameY="Outcome",
                             sandwich.nuisance=FALSE,
                             fay.adjustment=fay_adjustment,
                             ics=ICS,  
                             model.trt=as.formula(paste(variableEXP,"~Sex_G+Age_G+Degree_G+ Mean_Neighbor_Degree_G+ Assortativity_G+ LCC_Size_G +Mean_Component_Size_G+ Node_Component_Size_G",sep="")),
                             model.trt.clusterlevel=TRUE,
                             stepwise.trt=TRUE,
                             stepwise.augmentation=step_aug,
                             print.log=print.log)
  
  beta0[4,]<-c("aug all",augall$beta[1],sqrt(augall$var[1,1]),augall$beta[1]/sqrt(augall$var[1,1]))
  beta1[4,]<-c("aug all",augall$beta[2],sqrt(augall$var[2,2]),augall$beta[2]/sqrt(augall$var[2,2]))
  
  tablefinal<-rbind(beta0,beta1)
  tablefinal<-cbind(c(rep("beta0",4),rep("beta1",4)),tablefinal,rep(NA,4*2))
  for(i in 1:8){
    tablefinal[i,6]<-pt(abs(as.numeric(tablefinal[i,5])),length(unique(data.india$Village))-2,lower.tail=FALSE)
  }
  tablefinal[,3] <- round(as.numeric(tablefinal[,3]),3)
  tablefinal[,4] <- round(as.numeric(tablefinal[,4]),3)
  tablefinal[,5] <- round(as.numeric(tablefinal[,5]),3)
  tablefinal[,6] <- round(as.numeric(tablefinal[,6]),3)
  tablefinal<-as.data.frame(tablefinal)
  names(tablefinal)<-c("param","methods","estimate","se","wald","pval")
  
  return(list(beta0 = beta0, beta1 = beta1, augall=augall, aug=aug,noadj=noadj,gee=gee,tablefinal=tablefinal))
}


dostatistics2<-function(data.india,variableEXP,corstr="independence",print.log=TRUE){
  beta0<-matrix(NA,ncol=4,nrow=4)
  beta1<-matrix(NA,ncol=4,nrow=4)
  noadj<-geeDREstimationObs(formula<-as.formula(paste("Outcome~",variableEXP,sep=" ")),
                            id<-"Village",
                            data<-data.india,
                            family = gaussian,
                            pi.a = pi.a[1,variableEXP],
                            corstr = corstr,
                            nameTRT=variableEXP,
                            nameMISS="MISSING",
                            nameY="Outcome",
                            sandwich.nuisance=FALSE,
                            fay.adjustment=fay_adjustment,
                            ics=ICS,  
                            model.trt=NULL,
                            stepwise.trt=TRUE,
                            stepwise.augmentation=step_aug,
                            print.log=print.log)
  
  beta0[1,]<-c("unadj",noadj$beta[1],sqrt(noadj$var[1,1]),noadj$beta[1]/sqrt(noadj$var[1,1]))
  beta1[1,]<-c("unadj",noadj$beta[2],sqrt(noadj$var[2,2]),noadj$beta[2]/sqrt(noadj$var[2,2]))
  
  gee<-geeDREstimationObs(formula<-as.formula(paste("Outcome~",variableEXP,sep=" ")),
                          id<-"Village",
                          data<-data.india,
                          family = gaussian,
                          pi.a = pi.a[1,variableEXP],
                          corstr = corstr,
                          nameTRT=variableEXP,
                          nameMISS="MISSING",
                          nameY="Outcome",
                          sandwich.nuisance=FALSE,
                          fay.adjustment=fay_adjustment,
                          ics=ICS,  
                          model.trt=as.formula(paste(variableEXP,"~Sex_G+Age_G+Degree_G+ Mean_Neighbor_Degree_G+ Assortativity_G+ LCC_Size_G +Mean_Component_Size_G+ Node_Component_Size_G",sep="")),
                          model.trt.clusterlevel=TRUE,
                          stepwise.trt=TRUE,
                          stepwise.augmentation=TRUE,
                          print.log=print.log)
  
  beta0[2,]<-c("gee",gee$beta[1],sqrt(gee$var[1,1]),gee$beta[1]/sqrt(gee$var[1,1]))
  beta1[2,]<-c("gee",gee$beta[2],sqrt(gee$var[2,2]),gee$beta[2]/sqrt(gee$var[2,2]))
  
  aug<-geeDREstimationObs(formula<-as.formula(paste("Outcome~",variableEXP,sep=" ")),
                          id<-"Village",
                          data<-data.india,
                          family = gaussian,
                          pi.a = pi.a[1,variableEXP],
                          corstr = corstr,
                          model.augmentation.trt =as.formula(Outcome~Degree+Mean_Neighbor_Degree+Assortativity+LCC_Size+Mean_Component_Size+Node_Component_Size+Total_Neighbor_Seeds+Total_Cluster_Seeds+Mins+Sums),
                          model.augmentation.ctrl=as.formula(Outcome~Degree+Mean_Neighbor_Degree+Assortativity+LCC_Size+Mean_Component_Size+Node_Component_Size+Total_Neighbor_Seeds+Total_Cluster_Seeds+Mins+Sums),
                          nameTRT=variableEXP,
                          nameMISS="MISSING",
                          nameY="Outcome",
                          sandwich.nuisance=FALSE,
                          fay.adjustment=fay_adjustment,
                          ics=ICS,  
                          model.trt=as.formula(paste(variableEXP,"~Sex_G+Age_G+Degree_G+ Mean_Neighbor_Degree_G+ Assortativity_G+ LCC_Size_G +Mean_Component_Size_G+ Node_Component_Size_G",sep="")),
                          model.trt.clusterlevel=TRUE,
                          stepwise.trt=TRUE,
                          stepwise.augmentation=step_aug,
                          print.log=print.log)
  beta0[3,]<-c("aug",aug$beta[1],sqrt(aug$var[1,1]),aug$beta[1]/sqrt(aug$var[1,1]))
  beta1[3,]<-c("aug",aug$beta[2],sqrt(aug$var[2,2]),aug$beta[2]/sqrt(aug$var[2,2]))
  
  
  augall<-geeDREstimationObs(formula<-as.formula(paste("Outcome~",variableEXP,sep=" ")),
                             id<-"Village",
                             data<-data.india,
                             family = gaussian,
                             pi.a = pi.a[1,variableEXP],
                             corstr = corstr,
                             model.augmentation.trt =as.formula(Outcome~Sex+Age+Degree+Mean_Neighbor_Degree+Assortativity+LCC_Size+Mean_Component_Size+Node_Component_Size+Total_Neighbor_Seeds+Total_Cluster_Seeds+Mins+Sums),
                             model.augmentation.ctrl=as.formula(Outcome~Sex+Age+Degree+Mean_Neighbor_Degree+Assortativity+LCC_Size+Mean_Component_Size+Node_Component_Size+Total_Neighbor_Seeds+Total_Cluster_Seeds+Mins+Sums),
                             nameTRT=variableEXP,
                             nameMISS="MISSING",
                             nameY="Outcome",
                             sandwich.nuisance=FALSE,
                             fay.adjustment=fay_adjustment,
                             ics=ICS,  
                             model.trt=as.formula(paste(variableEXP,"~Sex_G+Age_G+Degree_G+Mean_Neighbor_Degree_G+Assortativity_G+LCC_Size_G+Mean_Component_Size_G+Node_Component_Size_G",sep="")),
                             model.trt.clusterlevel=TRUE,
                             stepwise.trt=TRUE,
                             stepwise.augmentation=step_aug,
                             print.log=print.log)
  beta0[4,]<-c("aug all",augall$beta[1],sqrt(augall$var[1,1]),augall$beta[1]/sqrt(augall$var[1,1]))
  beta1[4,]<-c("aug all",augall$beta[2],sqrt(augall$var[2,2]),augall$beta[2]/sqrt(augall$var[2,2]))
  
  
  tablefinal<-rbind(beta0,beta1)
  tablefinal<-cbind(c(rep("beta0",4),rep("beta1",4)),tablefinal,rep(NA,4*2))
  for(i in 1:8){
    tablefinal[i,6]<-pt(abs(as.numeric(tablefinal[i,5])),length(unique(data.india$Village))-2,lower.tail=FALSE)
  }
  tablefinal[,3]<-round(as.numeric(tablefinal[,3]),3)
  tablefinal[,4]<-round(as.numeric(tablefinal[,4]),3)
  tablefinal[,5]<-round(as.numeric(tablefinal[,5]),2)
  tablefinal[,6]<-round(as.numeric(tablefinal[,6]),3)
  
  tablefinal<-as.data.frame(tablefinal)
  names(tablefinal)<-c("param","methods","estimate","se","wald","pval")
  
  return(list(beta0=beta0,beta1=beta1,noadj=noadj,gee=gee,aug=aug,augall=augall,aug=aug,tablefinal=tablefinal))
}

EXP1.results.e_ind<-dostatistics(data.india,"EXP1",corstr="exchangeable")
xtable(EXP1.results.e_ind$tablefinal)
print(EXP1.results.e_ind)

EXP2.results.e_ind<-dostatistics2(data.india,"EXP2",corstr="exchangeable")
xtable(EXP2.results.e_ind$tablefinal)
print(EXP2.results.e_ind)

# Exposure 1
#summary(EXP1.results.e_ind$noadj)
xtable(summary(EXP1.results.e_ind$aug$om.model.trt)$coefficients, digits=4)
xtable(summary(EXP1.results.e_ind$aug$om.model.ctrl)$coefficients, digits=4)
xtable(summary(EXP1.results.e_ind$augall$om.model.trt)$coefficients, digits=4)
xtable(summary(EXP1.results.e_ind$augall$om.model.ctrl)$coefficients, digits=4)

# Exposure 2
#summary(EXP2.results.e_ind$noadj)
xtable(summary(EXP2.results.e_ind$aug$om.model.trt)$coefficients, digits=4)
xtable(summary(EXP2.results.e_ind$aug$om.model.ctrl)$coefficients, digits=4)
xtable(summary(EXP2.results.e_ind$augall$om.model.trt)$coefficients, digits=4)
xtable(summary(EXP2.results.e_ind$augall$om.model.ctrl)$coefficients, digits=4)

EXP1.results.e_ind$tablefinal
EXP2.results.e_ind$tablefinal



# Degree Distribution Plot
xx = sort(data.india$Degree, decreasing = TRUE)
y = table(xx)
x = as.numeric(names(y))
y=as.vector(y)
pdf("Degree_Distribution.pdf",width=8,height=6)
par(cex=1.2)
plot(log10(x),log10(1-cumsum(y)/sum(y)),pch=16,xlab="Log Degree",ylab=c("Log Survival"),main="Karnataka Dataset\nDegree Distribution")
dev.off()












library(xtable)
#setwd("C:/Users/Patrick/Desktop")

digits = 0
nc=15

setups_num = cbind(
  rep(0:1,each=2**5,times=2**0),
  rep(0:1,each=2**4,times=2**1),
  rep(0:1,each=2**3,times=2**2),
  rep(0:1,each=2**2,times=2**3),
  rep(0:1,each=2**1,times=2**4),
  rep(0:1,each=2**0,times=2**5)
)
setups=apply(setups_num,1,paste,collapse="")
colnames = c("UnAdj",paste("F_",1:12,sep=""),"All","Step")
modnames = c("(Intercept)", "High_Degree", "Powerlaw","Assortative","Communities","Deg_Infect","High_Baseline")
condnames = setups

betas        =     t(t(unlist(read.csv("empirical_betas.txt", header=F))));rownames(betas) <- condnames
biases       =     read.csv("biases.txt",          header=F);rownames(biases)       <- condnames;colnames(biases)       <- colnames
se_robs      =     read.csv("se_robs.txt",         header=F);rownames(se_robs)      <- condnames;colnames(se_robs)      <- colnames
se_emps      =     read.csv("se_emps.txt",         header=F);rownames(se_emps)      <- condnames;colnames(se_emps)      <- colnames
improvements =     read.csv("improvements.txt",    header=F);rownames(improvements) <- condnames;colnames(improvements) <- colnames
powers       = 100*read.csv("powers.txt",          header=F);rownames(powers)       <- condnames;colnames(powers)       <- colnames
coverages    = 100*read.csv("coverages.txt",       header=F);rownames(coverages)    <- condnames;colnames(coverages)    <- colnames

#Decide what to print.
inc_to_print = c(1,10,14,15)
biases       = biases[,inc_to_print]
se_robs      = se_robs[,inc_to_print]
se_emps      = se_emps[,inc_to_print]
improvements = improvements[,inc_to_print]
powers       = powers[,inc_to_print]
coverages    = coverages[,inc_to_print]
colnames     = colnames[inc_to_print]

betas_mod        <- array(NA,c(7))
biases_mod       <- array(NA,c(7,length(inc_to_print)));colnames(biases_mod)       <- colnames; rownames(biases_mod)       <- modnames
se_robs_mod      <- array(NA,c(7,length(inc_to_print)));colnames(se_robs_mod)      <- colnames; rownames(se_robs_mod)      <- modnames
se_emps_mod      <- array(NA,c(7,length(inc_to_print)));colnames(se_emps_mod)      <- colnames; rownames(se_emps_mod)      <- modnames
improvements_mod <- array(NA,c(7,length(inc_to_print)));colnames(improvements_mod) <- colnames; rownames(improvements_mod) <- modnames
powers_mod       <- array(NA,c(7,length(inc_to_print)));colnames(powers_mod)       <- colnames; rownames(powers_mod)       <- modnames
coverages_mod    <- array(NA,c(7,length(inc_to_print)));colnames(coverages_mod)    <- colnames; rownames(coverages_mod)    <- modnames

modification = function(y){
  data = data.frame(y=y, setups_num)
  return(summary(lm(y~.,data=data))$coef[1:7,1:2])
}

betas_mod            <- modification(betas)[,1];names(betas_mod) <- modnames
for(i in 1:length(inc_to_print)){
  biases_mod[,i]       <- modification(biases      [,i])[,1]
  se_robs_mod[,i]      <- modification(se_robs     [,i])[,1]
  se_emps_mod[,i]      <- modification(se_emps     [,i])[,1]
  improvements_mod[,i] <- modification(improvements[,i])[,1]
  powers_mod[,i]       <- modification(powers      [,i])[,1]
  coverages_mod[,i]    <- modification(coverages   [,i])[,1]
}

betas_marg           <- mean(betas,na.rm=T)
biases_marg          <- apply(biases,2,function(x){mean(x,na.rm=T)});      names(biases_marg)          <- colnames
percent_biases_marg  <- biases_marg/betas_marg
se_robs_marg         <- apply(se_robs,2,mean,na.rm=T);     names(se_robs_marg)         <- colnames
se_emps_marg         <- apply(se_emps,2,mean,na.rm=T);     names(se_emps_marg)         <- colnames
improvements_marg    <- apply(improvements,2,mean,na.rm=T);names(improvements_marg)    <- colnames
improvements_marg_sd <- apply(improvements,2,sd,na.rm=T);  names(improvements_marg_sd) <- colnames
powers_marg          <- apply(powers,2,mean,na.rm=T);      names(powers_marg)          <- colnames
powers_marg_sd       <- apply(powers,2,sd,na.rm=T);        names(powers_marg_sd)       <- colnames
coverages_marg       <- apply(coverages,2,mean,na.rm=T);   names(coverages_marg)       <- colnames
coverages_marg_sd    <- apply(coverages,2,sd,na.rm=T);     names(coverages_marg_sd)    <- colnames

powers_comb       <- paste(round(powers_marg), paste(paste("(",round(powers_marg_sd),sep=""),")",sep=""),sep="")
improvements_comb <- paste(round(improvements_marg), paste(paste("(",round(improvements_marg_sd),sep=""),")",sep=""),sep="")
coverages_comb <- paste(round(coverages_marg), paste(paste("(",round(coverages_marg_sd),sep=""),")",sep=""),sep="")

marginals <- rbind(
  round(percent_biases_marg,2),round(se_robs_marg,2),
  round(se_emps_marg,2),improvements_comb,
  round(powers_marg,1),coverages_comb
        );rownames(marginals) = c("Bias","RobSE","EmpSE","Improvement","Power","Coverage")


specific = 47
specifics = rbind(
  biases[specific,],
  se_robs[specific,],
  se_emps[specific,],
  improvements[specific,],
  powers[specific,],
  coverages[specific,]
);rownames(specifics) = c("Bias","RobSE","EmpSE","Improvement","Power","Coverage")

print(xtable(t(betas_marg)), file = "marginal_table.txt",append=F)
print(xtable(marginals),     file = "marginal_table.txt",append=T)
print(xtable(specifics),     file = "specific_table.txt",append=T)

print(xtable(t(t(betas_mod)),digits=2),            file = "modification_tables.txt",append=F)
print(xtable(biases_mod,digits=2),            file = "modification_tables.txt",append=T)
print(xtable(se_robs_mod,digits=2),           file = "modification_tables.txt",append=T)
print(xtable(se_emps_mod,digits=2),           file = "modification_tables.txt",append=T)
print(xtable(improvements_mod,digits=digits), file = "modification_tables.txt",append=T)
print(xtable(powers_mod,digits=digits),       file = "modification_tables.txt",append=T)
print(xtable(coverages_mod,digits=digits),    file = "modification_tables.txt",append=T)

print(xtable(betas,digits=2),                 file = "conditional_tables.txt",append=F)
print(xtable(biases,digits=2),                file = "conditional_tables.txt",append=T)
print(xtable(100*biases/betas,digits=digits), file = "conditional_tables.txt",append=T)
print(xtable(se_robs,digits=2),               file = "conditional_tables.txt",append=T)
print(xtable(se_emps,digits=2),               file = "conditional_tables.txt",append=T)
print(xtable(improvements,digits=digits),     file = "conditional_tables.txt",append=T)
print(xtable(powers,digits=digits),           file = "conditional_tables.txt",append=T)
print(xtable(coverages,digits=digits),        file = "conditional_tables.txt",append=T)


included <- read.csv("included_table.txt",header=FALSE)

included_mods = array(NA,c(7,12))
for(i in 1:12) included_mods[,i] = modification(included   [,i])[,1][1:7]
rownames(included_mods) <- modnames
print(xtable(included_mods,digits=2),            file = "inclusion_modifications.txt",append=T)

included_mods = array(NA,c(8,12))
included_sigs = array(NA,c(8,12))
included_both = array(NA,c(8,12))

for(col in 1:12){
  rownames(included_both) <- c(modnames[1], "Exposed", modnames[2:7])
  y=as.vector(as.matrix(included[,col+c(0,12)]))
  expanded  = rbind(cbind(1,setups_num),cbind(0,setups_num))
  data = data.frame(y=y, expanded)
  analysis = summary(lm(y~.,data=data))$coef[1:8,c(1,4)]
  included_mods[,col] = analysis[,1]
  included_sigs[,col] = ifelse(analysis[,2]<0.001, "*", " ")
  for(row in 1:8){
    included_both[row,col] = paste(round(included_mods[row,col],2),included_sigs[row,col],sep="")}
}
print(xtable(included_both), file = "inclusion_mods_and_sigs.txt",append=F)




pdf("Improvements.pdf",width=7,height=6)
mgp=c(2.5,1,0)
colors = c(rep(rgb(.8,.2,.2),4),rep(rgb(.2,.7,.2),4),rep(rgb(.2,.2,.7),4),rep(rgb(.4,.4,.4),2),rep(rgb(.8,.2,.2),2))
apply(improvements,2,mean)
plot(0,0,xaxt="n",xlim=c(1,nc-3),ylim=c(-20,90),xlab="Adjustment Covariate",ylab="Percent Reduction in RMSE",main="Range of Estimation Gains\nBy Adjustment Covariate")
axis(1,at=1:nc)
abline(h=0,col="gray",lty=2,lwd=1.5)
for(i in 1:(nc-1)){
  lines(c(i,i),quantile(improvements[,i+1],c(.1,.9)),lwd=2.5,col=colors[i])
  points(c(i,i),quantile(improvements[,i+1],c(0,1)),pch=16,cex=.5,col=colors[i])
  }
points(apply(improvements,2,mean)[2:nc],pch=18,cex=1.5,col=colors)
dev.off()



blue=rgb(.3,.4,1)
pdf("Exposure_probabilities.pdf",width=8,height=6)
par(cex=1.5)

for(setup in setups[1:3]){
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
  
  clusters = ncol(bins)/4#4 when using DRgeeOBS, and 2 when using geeDR.  You can then compare model fit across the two versions.
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


gray = rgb(.5,.5,.5)
blue = rgb(.2,.2,.9)
red = rgb(.9,.2,.2)
colors = c(blue, red)

labels = c(expression(X^(1)), expression(X^(2)), expression(X^(3)), expression(X^(4)), 
           expression(X^(5)), expression(X^(6)), expression(X^(7)), expression(X^(8)), 
           expression(X^(9)), expression(X^(10)), expression(X^(11)), expression(X^(12)))

pdf("percent_included.pdf",width=8,height=6)
par(cex=1.2)
plot(apply(included,2,median)[1:12],ylim=c(0,1),col=NA,xaxt="n",
     ylab="Fraction of Times Included",xlab="Covariate",main="Fraction of Times\nEach Covariate is Included")
axis(1,at=seq(1,11,2),labels=labels[seq(1,11,2)]);axis(1,at=seq(2,12,2),labels=labels[seq(2,12,2)])
abline(h=c(0,1),lty=2,col=rgb(.8,.8,.8))
for(j in 1:2){
  for(i in 1:12){
    lines(i+c(-1,1)*.2+(j-1.5)*.2,rep(quantile(included[,i+(j-1)*12],.25),2),lwd=2,col=colors[j])
    lines(i+c(-1,1)*.2+(j-1.5)*.2,rep(quantile(included[,i+(j-1)*12],.75),2),lwd=2,col=colors[j])  
    lines(rep(i+(j-1.5)*.2,2),quantile(included[,i+(j-1)*12],probs = c(.25,.75)),col=colors[j],lwd=2)
  }
  for(i in 1:12){  
    points(i+(j-1.5)*.2,min(included[,i+(j-1)*12]),col=colors[j],pch=18,cex=1)
    points(i+(j-1.5)*.2,max(included[,i+(j-1)*12]),col=colors[j],pch=18,cex=1)
  }
  points(1:12+(j-1.5)*.2,apply(included,2,median)[1:12+(j-1)*12],pch=18,cex=2,col=colors[j])
  points(1:12+(j-1.5)*.2,apply(included,2,median)[1:12+(j-1)*12],pch=18,cex=1,col="white")
}
par(cex=1.5)
legend("bottomright", legend=c("Exposed","Unexposed"), col=c(red, blue), pch=18, pt.cex=1.35,bg="white")
dev.off()



cors1 <- sqrt(abs(cor(included[,1:12])))
cors2 <- sqrt(abs(cor(included[,1:12+12])))
pdf("inclusion_correlations.pdf",width=7.2,height=8.2)
par(cex=1.5)
red = rgb(1,0,0);blue=rgb(0,0,1);gray=rgb(.9,.9,.9)
plot(0,0,col=NA,xlim=c(0.5,12.5),ylim=c(.5,12.5),xlab=c("Covariate"),ylab="Covariate",main="Correlation Matrix",xaxt="n",yaxt="n")
#axis(1,at=seq(1,11,2),labels=labels[seq(1,11,2)]);axis(1,at=seq(2,12,2),labels=labels[seq(2,12,2)])
axis(1,at=1:6*2,labels=1:6*2);axis(1,at=0:5*2+1,labels=0:5*2+1)
#axis(1,at=4,labels=expression(X^12))
axis(2,at=seq(1,12,2),labels=c(12,10,8,6,4,2));axis(2,at=seq(2,12,2),labels=c(11,9,7,5,3,1))
abline(13,-1,col=gray)
for(i in 1:12){
  for(j in 1:12){
    if(i!=j){
      x=min(cors1[i,j],cors2[i,j]);y=max(cors1[i,j],cors2[i,j])
      points(i,13-j,cex=y*3.7,pch=16,col=ifelse(x==cors1[i,j],blue,red))
      points(i,13-j,cex=max(0,y*3.7-.5),pch=16,col="white")      
      points(i,13-j,cex=x*3.7,pch=16,col=ifelse(y==cors1[i,j],blue,red))
      points(i,13-j,cex=max(0,x*3.7-.5),pch=16,col="white")
      points(i,13-j,cex=3.7,pch=1,col=rgb(.9,.9,.9),lwd=1)
      }
    }
}
box()
dev.off()

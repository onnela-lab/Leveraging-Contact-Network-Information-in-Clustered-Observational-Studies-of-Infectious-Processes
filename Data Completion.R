data = read.table("incomplete_data.txt",header=TRUE)
leaders = c()
N       = c()
SHGs    = c()
for(village in unique(data$Village)){
  leaders=rbind(leaders,apply(data[which(data$Village==village & data$Is_Leader>0),c("Is_Leader", "Degree")], 2,sum))
  N = c(N,length(which(data$Village==village)))
  SHGs = c(SHGs,mean(data[which(data$Village==village & data$Self_Help != -1),"Self_Help"], na.rm=TRUE))
}
EXP1 = ifelse(leaders[,1]/N > as.numeric(quantile(leaders[,1]/N, 0.75, na.rm=TRUE)), 1, 0) # 1 if fraction of leaders high.
EXP2 = ifelse(SHGs>as.numeric(quantile(SHGs, 0.75, na.rm=TRUE)),1,0) # 1 if fraction of participation in concurrent self help groups is high.
counter = 0
for(village in unique(data$Village)){
  counter = counter + 1
  data[which(data$Village==village),"EXP1"] = EXP1[counter]
  data[which(data$Village==village),"EXP2"] = EXP2[counter]
}
write.table(data, "complete_data.txt", col.names=TRUE, row.names = FALSE)

m0 = apply(normalized[,1:3],1,mean)
m4 = apply(normalized[,4:6],1,mean)
m6 = apply(normalized[,7:8],1,mean)
m8 = apply(normalized[,9:11],1,mean)
m10 = apply(normalized[,12:14],1,mean)
m12 = apply(normalized[,15:16],1,mean)
mean_exps <- cbind(m0,m4,m6,m8,m10,m12)
write.csv(mean_exps,"mean_exps.csv")

xy = apply(mean_exps,1,function(x)lm(x[1:6]~c(1:6)))
#summary(xy[[1]])$coefficients[2,4]    ---extract the p_val for the 1st element in xy
lm_p_val = sapply(xy,function(x) summary(x)$coefficients[2,4])
mean_exps_m1 = cbind(mean_exps,lm_p_val)
#extract the p value

lm_q_val = length(lm_p_val)*mean_exps_m1[,7]/rank(mean_exps_m1[,7])
mean_exps_m2 = cbind(mean_exps_m1,lm_q_val)
write.csv(mean_exps_m2,"mean_exps_m2.csv")
#calculate the q value using FDR

library(mouse4302.db)
source("report10.R")
n_lm = report_10genes(lm_q_val)
write.table(n_lm,"lm_sig_gene.txt",row.names=FALSE)
#10 most significant genes by FDR

save.image("mean_coefficient.RData")
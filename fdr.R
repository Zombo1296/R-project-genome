#2. Calculate the FDR using the t-test results from Assignment 3. Append the ttest.csv file by adding the FDR values
source("FDR.R")
q_val4 = fdr(p_val4)
q_val6 = fdr(p_val6)
q_val8 = fdr(p_val8)
q_val10 = fdr(p_val10)
q_val12 = fdr(p_val12)
q_vals <- cbind(q_val4,q_val6,q_val8,q_val10,q_val12)
write.csv(q_vals,"q_vals.csv")
combined_ttest <- cbind(ttest,q_vals)

write.csv(combined_ttest,"combined_ttest.csv")

#3. Report the 10 most significant genes for each comparison in sig_gene.txt file
library(mouse4302.db)
source("report10.R")
n4_BH = report_10genes(q_val4)
n6_BH = report_10genes(q_val6)
n8_BH = report_10genes(q_val8)
n10_BH = report_10genes(q_val10)
n12_BH = report_10genes(q_val12)
combined_gene <- cbind(n4_BH,n6_BH,n8_BH,n10_BH,n12_BH)
write.table(combined_gene,"sig_gene.txt",row.names=FALSE)

save.image("FDR.RData")

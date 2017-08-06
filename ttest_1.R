library(mouse4302.db)

#norm = read.table("normalized.csv")
#can be omitted by restoring saved image

#fc4 <- apply(normalized[,4:6],1,mean) / apply(normalized[,1:3],1,mean)
#probe.fc4 = names(fc4[fc4>1.5])
#gene.fc4 = unlist(mget(probe.fc4,mouse4302SYMBOL))
#gene.sig = unique(gene.fc4[!is.na(gene.fc4)])

#Comparison between Day0 and Day4
fun4 = function(x) {
t.test(x[1:3],x[4:6])
}
#not x[,1:3] but x[1:3]
x4 = apply(normalized,1,fun4)
p_val4 = unlist(lapply(x4,'[[',3))
#generate a list('x4') then combine the 3rd value of each component to a vector 
probe.p_val4 = names(p_val4[p_val4<0.001])
#p-value < alpha(0.001) as cutoff rather than fold_change > 1.5
gene.p_val4 = unlist(mget(probe.p_val4,mouse4302SYMBOL))
gene.sig4 = unique(gene.p_val4[!is.na(gene.p_val4)])
n_of_s4 = length(gene.sig4)

#Comparison between Day0 and Day6
fun6 = function(x) {
t.test(x[1:3],x[7:8])
}
x6 = apply(normalized,1,fun6)
p_val6 = unlist(lapply(x6,'[[',3))
probe.p_val6 = names(p_val6[p_val6<0.001])
gene.p_val6 = unlist(mget(probe.p_val6,mouse4302SYMBOL))
gene.sig6 = unique(gene.p_val6[!is.na(gene.p_val6)])
n_of_s6 = length(gene.sig6)

#Comparison between Day0 and Day8
fun8 = function(x) {
t.test(x[1:3],x[9:11])
}
x8 = apply(normalized,1,fun8)
p_val8 = unlist(lapply(x8,'[[',3))
probe.p_val8 = names(p_val8[p_val8<0.001])
gene.p_val8 = unlist(mget(probe.p_val8,mouse4302SYMBOL))
gene.sig8 = unique(gene.p_val8[!is.na(gene.p_val8)])
n_of_s8 = length(gene.sig8)

#Comparison between Day0 and Day10
fun10 = function(x) {
t.test(x[1:3],x[12:14])
}
x10 = apply(normalized,1,fun10)
p_val10 = unlist(lapply(x10,'[[',3))
probe.p_val10 = names(p_val10[p_val10<0.001])
gene.p_val10 = unlist(mget(probe.p_val10,mouse4302SYMBOL))
gene.sig10 = unique(gene.p_val10[!is.na(gene.p_val10)])
n_of_s10 = length(gene.sig10)

#Comparison between Day0 and Day12
fun12 = function(x) {
t.test(x[1:3],x[15:16])
}
x12 = apply(normalized,1,fun12)
p_val12 = unlist(lapply(x12,'[[',3))
probe.p_val12 = names(p_val12[p_val12<0.001])
gene.p_val12 = unlist(mget(probe.p_val12,mouse4302SYMBOL))
gene.sig12 = unique(gene.p_val12[!is.na(gene.p_val12)])
n_of_s12 = length(gene.sig12)


ttest <- cbind(p_val4, p_val6, p_val8, p_val10, p_val12)
write.csv(ttest, "ttest.csv")

compare_d46 = intersect(gene.sig4,gene.sig6)
n_of_c46 = length(compare_d46)
compare_d48 = intersect(gene.sig4,gene.sig8)
n_of_c48 = length(compare_d48)
compare_d410 = intersect(gene.sig4,gene.sig10)
n_of_c410 = length(compare_d410)
compare_d412 = intersect(gene.sig4,gene.sig12)
n_of_c412 = length(compare_d412)
compare_d68 = intersect(gene.sig6,gene.sig8)
n_of_c68 = length(compare_d68)
compare_d610 = intersect(gene.sig6,gene.sig10)
n_of_c610 = length(compare_d610)
compare_d612 = intersect(gene.sig6,gene.sig12)
n_of_c612 = length(compare_d612)
compare_d810 = intersect(gene.sig8,gene.sig10)
n_of_c810 = length(compare_d810)
compare_d812 = intersect(gene.sig8,gene.sig12)
n_of_c812 = length(compare_d812)
compare_d1012 = intersect(gene.sig10,gene.sig12)
n_of_c1012 = length(compare_d1012)

t_significant <- rbind(n_of_s4, n_of_s6, n_of_s8, n_of_s10, n_of_s12, n_of_c46, n_of_c48, n_of_c410, n_of_c412, n_of_c68, n_of_c610, n_of_c612, n_of_c810, n_of_c812, n_of_c1012)
write.csv(t_significant, "t_significant.csv")

save.image("t-test_1.RData")

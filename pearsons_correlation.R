normalized_copy = normalized
colname(normalized_copy)=c("Day0_1","Day0_2","Day0_3","Day4_1","Day4_2","Day4_3","Day6_1","Day6_2","Day8_1","Day8_2","Day8_3","Day10_1","Day10_2","Day10_3","Day12_1","Day12_2")

### Euclidean distance
d1 = dist(t(normalized_copy))
clust.E1 = hclust(d1, method="single")
clust.E2 = hclust(d1, method="complete")
clust.E3 = hclust(d1, method="average")
clust.E4 = hclust(d1, method="centroid")
clust.E5 = hclust(d1, method="ward.D")
pdf("Euclidean distance.pdf")
plot(clust.E1, main="Euc Simple")
plot(clust.E2, main="Euc Complete")
plot(clust.E3, main="Euc Average")
plot(clust.E4, main="Euc Centroid")
plot(clust.E5, main="Euc Ward")
dev.off()

### Pearsons's correlation
d2 = as.dist(1-cor(normalized_copy, method="pearson", use="na.or.complete"))
clust.P1 = hclust(d2, method="single")
clust.P2 = hclust(d2, method="complete")
clust.P3 = hclust(d2, method="average")
clust.P4 = hclust(d2, method="centroid")
clust.P5 = hclust(d2, method="ward.D")
pdf("Pearson distance.pdf")
plot(clust.P1, main="Euc Simple")
plot(clust.P2, main="Euc Complete")
plot(clust.P3, main="Euc Average")
plot(clust.P4, main="Euc Centroid")
plot(clust.P5, main="Euc Ward")
dev.off()

save.image("Pearsons_correlation.RData")
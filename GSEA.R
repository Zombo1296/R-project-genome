# Heatmap
library(gplots)
all.sd = apply(normalized, 1, sd)
temp = normalized[all.sd>1,]
temp.mean = apply(temp, 1, mean)
temp.sd = apply(temp, 1, sd)
x = (temp - temp.mean)*(2/temp.sd)

distfun = function(x) {as.dist(1-cor(t(x), method="pearson"))}
clustfun = function(x) {hclust(x, method="ward.D")}
pdf("heatmap_1.pdf")
heatmap.2(t(x), hclustfun=clustfun, distfun=distfun,col=greenred(75), trace="none", dendrogram="row",density.info="none", margins=c(4, 7), labCol="", key=T, symkey=F)
dev.off()

# Extract information from 'c5.bp.v5.0.symbols.gmt' into a list
lines = readLines("c5.bp.v5.0.symbols.gmt")
lines.list = strsplit(lines, "\t")
n = length(lines)

# Remove column 1/2 from list
for (i in 1:n) {
    names(lines.list)[i] = lines.list[[i]][1]
    lines.list[[i]] = lines.list[[i]][-c(1,2)]
}

# Get names for all the 45000 genes, including redundant genes
library(mouse4302.db)
gene.normalized.1 = unlist(mget(rownames(normalized),mouse4302SYMBOL))
#gene.normalized.2 = unique(gene.normalized.1[!is.na(gene.normalized.1)])
gene.normalized.3 = toupper(gene.normalized.1)

# Installing the GSA package, refer to the following address: http://stackoverflow.com/questions/1474081/how-do-i-install-an-r-package-from-source
# Compare Day0-Day4 (1,1,1,2,2,2) as an example
library(GSA)
GObp.gsa = GSA(normalized[,1:6], c(1,1,1,2,2,2), genesets=lines.list, genenames=gene.normalized.3, resp.type="Two class unpaired") 

GObp.hi = names(lines.list)[GObp.gsa$pvalues.hi<0.05]
GObp.hi.p = GObp.gsa$pvalues.hi[GObp.gsa$pvalues.hi<0.05]
GObp.hi.p = GObp.hi.p[!is.na(GObp.hi)]
GObp.hi = GObp.hi[!is.na(GObp.hi)]

GObp.lo = names(lines.list)[GObp.gsa$pvalues.lo<0.05]
GObp.lo.p = GObp.gsa$pvalues.lo[GObp.gsa$pvalues.lo<0.05]
GObp.lo.p = GObp.lo.p[!is.na(GObp.lo)]
GObp.lo = GObp.lo[!is.na(GObp.lo)]

n1 = list(up=GObp.hi,down=GObp.lo)
n2 = max(length(GObp.hi),length(GObp.lo))
n3 = sapply(n1,'[',1:n2)
write.csv(n3,"GSA_Day0_vs_Day6.csv")
rm(GObp.hi,GObp.lo,n1,n2,n3)
	
# A function to process several comparisons (needs improvement)
day0=normalized[,1:3]
day4=normalized[,4:6]
day6=normalized[,7:8]
day8=normalized[,9:11]
day10=normalized[,12:14]
day12=normalized[,15:16]
Xs = list(x4=cbind(day0,day4),x6=cbind(day0,day6),x8=cbind(day0,day8),x10=cbind(day0,day10),x12=cbind(day0,day12))

# Biological Process/ KEGG pathway/ Molecular function
data1 = readLines("c5.bp.v5.0.symbols.gmt")				#826 lines
data2 = readLines("c5.mf.v5.0.symbols.gmt")				#397 lines
data3 = readLines("c2.cp.kegg.v5.0.symbols.gmt")		#187 lines
data.sum1 = list(data1,data2,data3)
data.fun=function(x) {
    lines.list = strsplit(x, "\t")
    n = length(x)
    for (i in 1:n) {
        names(lines.list)[i] = lines.list[[i]][1]
        lines.list[[i]] = lines.list[[i]][-c(1,2)]
    }
    c(lines.list)
}
data.sum2 = sapply(data.sum1,data.fun)

normalized = read.csv("normalized.csv",header=TRUE,row.names=1,sep=",")
#source("http://bioconductor.org/biocLite.R")
#biocLite("mouse4302.db")
library(mouse4302.db)
gene.normalized.1 = unlist(mget(rownames(normalized),mouse4302SYMBOL))
gene.normalized.2 = toupper(gene.normalized.1)
day0=normalized[,1:3]
day4=normalized[,4:6]
day6=normalized[,7:8]
day8=normalized[,9:11]
day10=normalized[,12:14]
day12=normalized[,15:16]
Xs = list(x4=cbind(day0,day4),x6=cbind(day0,day6),x8=cbind(day0,day8),x10=cbind(day0,day10),x12=cbind(day0,day12))

#install.packages("http://cran.r-project.org/src/contrib/GSA_1.03.tar.gz", repo=NULL, type="source")
library(GSA)
function.GSA = function(x) {
	if (ncol(x)==5)  y=c(1,1,1,2,2)
	else  y=c(1,1,1,2,2,2)
	aa = 1
	bb = 1
	cc = list()
	for (i in 1:3) {
		gsa.temp1 = GSA(x, y, genesets=data.sum2[[aa]], genenames=gene.normalized.2, resp.type="Two class unpaired", nperms=200)
		gsa.temp2 = GSA.listsets(gsa.temp1, geneset.names=names(data.sum2[[aa]]),FDRcut=0.05)			#class(gsa.2) = list with 5 elements
		gsa.hi = gsa.temp2[[3]][,2]
		gsa.lo = gsa.temp2[[2]][,2]
		cc[[bb]] = list(gsa.hi,gsa.lo)			# a list of vectors
		aa = aa+1
		bb = bb+1
	}
	return (cc)
}

GSA_results = sapply(Xs,function.GSA)			# length(GSA_results)=15, length(GSA_results[[1]])=2
GSA.summary = sapply(GSA_results,function(x) sapply(x,'[',1:50))
write.csv(GSA.summary,"GSA.summary.csv")
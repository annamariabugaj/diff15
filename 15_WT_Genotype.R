library(MASS)
library(class)
library(cluster)
library(impute)
library(Hmisc)
library(splines)
library(WGCNA)

source("NetworkFunctions.txt")

dat1=read.csv("subtable_significant_genes_WT_VSD.csv",header=T)
attach(dat1)
dim(dat1)
head(dat1)

names(dat1);
head(dat1)

dat1 = dat1[, 2:44]
head(dat1)

# This data frame contains the gene expression data.
# By our convention, columns are genes and rows are samples.
datExpr=data.frame(t(dat1[, 2:43]))
dim(datExpr)
class(datExpr)
dimnames(datExpr)[[1]]
names(datExpr) = dat1$X;
rownames(datExpr) = names(dat1[-c(1)]);
dim(datExpr)


indexWT=c(22:42)
indexN3=c(1:21)
indexWT
indexN3

powers1 = c(c(1:10), seq(from=12, to=20, by=2))

RpowerTable=PickSoftThreshold(datExpr[indexWT,], powervector=powers1,
                              RsquaredCut=0.80)[[2]]

collect_garbage()
cex1=0.7
par(mfrow=c(1,2))
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="
     Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],
     labels=powers1,cex=cex1,col="red")
# This line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
plot(RpowerTable[,1], RpowerTable[,5],xlab="Soft Threshold (power)",ylab="Mean
     Connectivity", type="n")
text(RpowerTable[,1], RpowerTable[,5], labels=powers1, cex=cex1,col="red")

powerHuman=6

RpowerTableChimp=PickSoftThreshold(datExpr[indexN3,], powervector=powers1,
                                   RsquaredCut=0.85)[[2]]

powerChimp=powerHuman

## Calculation of the network (adjacency matrix) by raising the absolute value of the correlation #matrix to a power (soft-thresholding with the power adjacency function).
# Human network:
AdjMatHuman =  abs(cor(datExpr[indexWT,] ,use="p"))^powerHuman
diag(AdjMatHuman)=0
# Chimp network:
AdjMatChimp =  abs(cor(datExpr[indexN3,]  ,use="p"))^powerChimp
diag(AdjMatChimp)=0

save(AdjMatHuman, AdjMatChimp, datExpr, file = "WGCNA_Diff_15_adjMat.RData")
## Calculation of the whole network connectivity k:
ConnectivityHuman <- apply(AdjMatHuman,1,sum)
ConnectivityChimp <- apply(AdjMatChimp,1,sum)

## Depiction of scale-free topology:
# The black curve corresponds to scale-free topology and the red curve corresponds to truncated #scale-free topology.
par(mfrow=c(2,1))
ScaleFreePlot1(ConnectivityHuman, AF1=paste("Wildtype ,","power=6"),truncated1=T)
ScaleFreePlot1(ConnectivityChimp, AF1=paste("N3-/- ,","power=6"),truncated1=T)                      

## Scaling k to lie between 0 and 1:
ConnectivityHuman=ConnectivityHuman/max(ConnectivityHuman)
ConnectivityChimp=ConnectivityChimp/max(ConnectivityChimp)
## Comparing gene expression human and chimp:
ExpressionHuman=apply(datExpr[indexWT,],2,mean)
ExpressionHuman=ExpressionHuman/max(ExpressionHuman)
ExpressionChimp=apply(datExpr[indexN3,],2,mean)
ExpressionChimp=ExpressionChimp/max(ExpressionChimp)

cor.test(ExpressionHuman,ExpressionChimp,method="s",use="p", exact = FALSE)

# Now we form a scatter plot between the chimp and human expression values:
par(mfrow=c(1,1))
plot(ExpressionHuman,ExpressionChimp,main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT gene expression", ylab="N3-/- gene expression", sub="rho = 0.99",
     cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ExpressionChimp~ExpressionHuman),col="red",lwd=2)

## Comparing network connectivity between human and chimp:
cor.test(ConnectivityHuman,ConnectivityChimp,method="s",use="p", exact = FALSE)

par(mfrow=c(1,1))
plot(ConnectivityHuman,ConnectivityChimp,main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT network connectivity", ylab="N3-/- network connectivity", sub="rho =
0.82", cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ConnectivityChimp~ConnectivityHuman),col="red",lwd=2)

# As a pre-processing step towards module construction, we restrict the network to genes with #reasonably high connectivity. This does not lead to a big loss of information since module #genes tend to have high connectivity (7). Toward this end, consider the #median connectivity in human and chimp:
median(ConnectivityHuman)
median(ConnectivityChimp)

# This motivates us to restrict the analysis to genes with k > 0.1 in either human or # chimp:
minconnections=.2
rest1=ConnectivityChimp>minconnections | ConnectivityHuman>minconnections
table(rest1)
AdjMatChimprest=AdjMatChimp[rest1,rest1]
AdjMatHumanrest=AdjMatHuman[rest1,rest1]
rm(AdjMatChimp);
rm(AdjMatHuman); 
collect_garbage()

#Module Construction
# The topological overlap of two nodes reflects their similarity in terms of the commonality of #the nodes they connect to, see (6, 10).
## Creating distance matrices based upon the topological overlap of 9630 genes for humans and #chimpanzees:
distTOMChimp <- TOMdist1(AdjMatChimprest)
distTOMHuman <- TOMdist1(AdjMatHumanrest)
collect_garbage()
# here qwe test for the whole network without restricting connectivity
#distTOMChimp <- TOMdist1(AdjMatChimp)
#distTOMHuman <- TOMdist1(AdjMatHuman)
#collect_garbage()
# To group genes with coherent expression profiles into modules, we use average linkage # hierarchical clustering, which uses the topological overlap measure as dissimilarity.
## Performing average linkage hierarchical clustering using these distance matrices:
hierTOMHuman <- hclust(as.dist(distTOMHuman),method="average")
hierTOMChimp <- hclust(as.dist(distTOMChimp),method="average")


# install.packages("flashClust")
# library(flashClust)
# hierTOMHuman <- flashClust(as.dist(distTOMHuman),method="average")
# hierTOMChimp <- flashClust(as.dist(distTOMChimp),method="average")


## Displaying dendrograms:
par(mfrow=c(1,2))
plot(hierTOMHuman,labels=F,main="Wildtype")
plot(hierTOMChimp,labels=F,main="N3-/-")

## Assign colors to modules based upon the height cutoff (h1) and minimum module size #(minsize1):
# Once a dendrogram is obtained from a hierarchical clustering method, we need to choose a
# height cutoff h1 in order to arrive at a clustering. It is a judgement call where to cut the tree
# branches. The height cut-off can be found by inspection: a height cutoff value is chosen in the # dendrogram such that some of the resulting branches correspond to the discrete diagonal blocks # (modules) in the TOM plot.
colorh1=as.character(modulecolor2(hierTOMHuman,h1=.7,minsize1=30))
table(colorh1)

colorh2=as.character(modulecolor2(hierTOMChimp,h1=.7,minsize1=30))
table(colorh2)

colorh12=as.character(modulecolor2(hierTOMHuman,h1=.62,minsize1=30))
table(colorh12)

colorh22=as.character(modulecolor2(hierTOMChimp,h1=.62,minsize1=30))
table(colorh22)

# minsize specifies that a module should contain at least 30 genes.
dev.off()
par(mfrow=c(3,2),mar=c(2,2,2,2))
plot(hierTOMHuman,main="WT hippocampus",labels=F)
abline(h=.62,col="red")
plot(hierTOMChimp,main="N3-/- hippocampus",labels=F)
hclustplot1(hierTOMHuman,colorh12,title1="WT network, WT colors")
hclustplot1(hierTOMChimp,colorh12,title1="N3-/- network, WT colors")
hclustplotn(hierTOMHuman,colorh12)
hclustplotn(hierTOMChimp,colorh12)

#We now create color vectors for all genes in the analysis. 
#The genes that did not pass the minimum connectivity threshold will be assigned the color grey, 
#the color of unassigned genes. 
#At the end we save the module color information; the data is used in Section 2.b.

colorh=as.character(colorh1)
colorhALL=rep("grey", length(ConnectivityHuman))
colorhALL[rest1]=as.character(colorh)
colorh=as.character(colorh2)
colorh2ALL=rep("grey", length(ConnectivityChimp))
colorh2ALL[rest1]=as.character(colorh)
# Give the module colors more descriptive names
colorHuman = colorhALL;
colorChimp = colorh2ALL;
# Also save the indicator of genes that pass the connectivity threshold
inNetwork = rest1
# Save module colors (assignments) and the indicator
save(colorHuman, colorChimp, inNetwork,
     file = "HumanChimp-OldhamAnalysis-colorHuman-colorChimp-inNetwork.RData")


################# ERROR ###############################
TOMplot(distTOMHuman,  hierTOMHuman , colorh)
## Error: C stack usage  7955296 is too close to the limit


## Creating a classical multi-dimensional scaling plots for humans and chimpanzees as another #means of network representation:
cmdChimp=cmdscale(as.dist(distTOMChimp),3)
cmdHuman=cmdscale(as.dist(distTOMHuman),3)
collect_garbage()
pairs(cmdHuman, col=as.character(colorh),  main="WT network, WT colors",
      ylim=c(-.55,.55), xlim=c(-.55, .55))

pairs(cmdChimp, col=as.character(colorh),  main="N3 network, WT colors",
      ylim=c(-.55,.55), xlim=c(-.55, .55))

colorN3 <- as.character(colorh2)

pairs(cmdChimp, col=as.character(colorN3),  main="N3 network, N3 colors",
      ylim=c(-.55,.55), xlim=c(-.55, .55))


# 2 dimensional plot

cmd1 = cmdscale(as.dist(distTOMHuman), 2)
class(cmd1)
head(cmd1)
par(mfrow = c(1, 1))
plot(cmd1, col=as.character(moduleColors), main="MDS plot Wildtype",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")


cmd2 = cmdscale(as.dist(distTOMChimp), 2)
class(cmd2)
head(cmd2)
par(mfrow = c(1, 1))
plot(cmd2, col=as.character(moduleColors), main="MDS plot N3-/-",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")

# Creating heatmaps for each module with samples ordered according to brain region:

indexN3
indexWT

byregionsdat1 = c(23,26,29,32,35,38,41,24,27,30,33,36,39,42,25,28,31,34,37,40,43,2,5,8,11,14,17,20,3,6,9,12,15,18,21,4,7,10,13,16,19,22)
#byregionsdat1 = c(23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 2, 3,4,5,6,7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,22)
orderbyregion=c(1,2,3,4,5,6,7,1,2,3,4,5,6,7,1,2,3,4,5,6,7,1,2,3,4,5,6,7,1,2,3,4,5,6,7,1,2,3,4,5,6,7)
orderbyregion2 = c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)
ordersamples2=order(orderbyregion)
ordersamples3 = order(orderbyregion2)

datExprHumanrest=datExpr[indexWT,rest1]
datExprChimprest=datExpr[indexN3,rest1]

# Now we create a heatmap of the turquoise module genes (rows are genes, columns
#are microarray samples:
library(pheatmap)

par(mar=c(2,3,4,3))
par(oma=c(0,0,2,0))
whichmodule="turquoise"
datcombined=data.frame(rbind(datExprHumanrest , datExprChimprest))[,
                                                                   colorh==whichmodule]
mat1 <- t(scale(datcombined[ordersamples3,]))


data_WT = data.frame(datExprHumanrest)[, colorh12==whichmodule]
mat2 <- t(scale(data_WT))
head(mat2)
# we order the cpolumns in WT table by age and region
order_WT <- as.array(c("WTP81_CA1", "WTP82_CA1", "WTP83_CA1", "WTP84_CA1","WTP81_CA3", "WTP82_CA3", "WTP83_CA3", "WTP84_CA3","WTP81_DG", "WTP82_DG", "WTP83_DG", "WTP84_DG", "SM625_CA1", "SM626_CA1", "SM627_CA1", "SM625_CA3", "SM626_CA3", "SM627_CA3", "SM625_DG", "SM626_DG", "SM627_DG"))
# we order columns in the main dataframe
mat3 <- mat2[ , order_WT]
head(mat3)


pheatmap(mat3, color = colorRampPalette(c("light sky blue", "black", "yellow"))(30), 
         scale = "row", cluster_rows = TRUE,cluster_cols = FALSE,
         show_rownames = T,
         fontsize = 5,
         main = whichmodule)

#repvec=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
#labelHuman=as.character(rep(repvec,1))
#labelHuman
par(mfrow=c(1,1))
plotMat(t(scale(datcombined[ordersamples2,])) , rlabels=  dimnames(datcombined)[[2]]
        ,clabels=dimnames(dat1[orderbyregion])[[2]], ccols=labelHuman,rcols=whichmodule)



# The module eigengenes are used to summarize the expression profiles in each module.
# To compute the module eigengene for each module, we use the first principal component #estimated using the singular value decomposition:
PC1=moduleEigengenes(datExpr[1:42,rest1],colorh1)
#PC2= ModulePrinComps1(datExpr[1:42,rest1],colorh1)
# PC1 is a list with 2 components. The first component is a data frame whose columns #correspond to the module eigengenes.
# Note that the column names list the modules (PC=first principal component=module
# eigengene)
#names(PC2[[1]])
#names(PC2[[3]])
#signif(PC2[[3]],2)

PC=moduleEigengenes(datExpr[1:42,rest1],
                 colorh1,
                 impute = TRUE,
                 excludeGrey = FALSE,
                 nPC = 5)
names(PC[[1]])
PC[3]
head(PC)

#write.table(PC1[[1]],file="First_PC_grey.csv", sep=",",
#            row.names=dimnames(datExpr)[[1]][1:42])
 
dimnames(datExpr)[[1]]


### Finding hub genes 
hubsWT = chooseTopHubInEachModule(
  datExprHumanrest, 
  colorh1, 
  omitColors = "grey", 
  power = 6, 
  type = "unsigned")

hubsWT

hubsN3 = chooseTopHubInEachModule(
  datExprChimprest, 
  colorh1, 
  omitColors = "grey", 
  power = 6, 
  type = "unsigned")

hubsN3



## To calculate within module connectivity (kin) for all genes in each module for human and #chimp:
# 1 Black
AdjMatChimprestblack=AdjMatChimprest[colorh1=="black",colorh1=="black"]
ConnectivityChimpblack=apply(AdjMatChimprestblack,1,sum)
AdjMatHumanrestblack=AdjMatHumanrest[colorh1=="black",colorh1=="black"]
ConnectivityHumanblack=apply(AdjMatHumanrestblack,1,sum)
ConnectivityChimpblack=ConnectivityChimpblack/max(ConnectivityChimpblack)
ConnectivityHumanblack=ConnectivityHumanblack/max(ConnectivityHumanblack)
#2 Blue
AdjMatChimprestblue=AdjMatChimprest[colorh1=="blue",colorh1=="blue"]
ConnectivityChimpblue=apply(AdjMatChimprestblue,1,sum)
AdjMatHumanrestblue=AdjMatHumanrest[colorh1=="blue",colorh1=="blue"]
ConnectivityHumanblue=apply(AdjMatHumanrestblue,1,sum)
ConnectivityChimpblue=ConnectivityChimpblue/max(ConnectivityChimpblue)
ConnectivityHumanblue=ConnectivityHumanblue/max(ConnectivityHumanblue)
#3 Brown
AdjMatChimprestbrown=AdjMatChimprest[colorh1=="brown",colorh1=="brown"]
ConnectivityChimpbrown=apply(AdjMatChimprestbrown,1,sum)
AdjMatHumanrestbrown=AdjMatHumanrest[colorh1=="brown",colorh1=="brown"]

ConnectivityHumanbrown=apply(AdjMatHumanrestbrown,1,sum)
ConnectivityChimpbrown=ConnectivityChimpbrown/max(ConnectivityChimpbrown)
ConnectivityHumanbrown=ConnectivityHumanbrown/max(ConnectivityHumanbrown)
#4 Green
AdjMatChimprestgreen=AdjMatChimprest[colorh1=="green",colorh1=="green"]
ConnectivityChimpgreen=apply(AdjMatChimprestgreen,1,sum)
AdjMatHumanrestgreen=AdjMatHumanrest[colorh1=="green",colorh1=="green"]
ConnectivityHumangreen=apply(AdjMatHumanrestgreen,1,sum)
ConnectivityChimpgreen=ConnectivityChimpgreen/max(ConnectivityChimpgreen)
ConnectivityHumangreen=ConnectivityHumangreen/max(ConnectivityHumangreen)
#5 Purple
AdjMatChimprestpurple=AdjMatChimprest[colorh1=="purple",colorh1=="purple"]
ConnectivityChimppurple=apply(AdjMatChimprestpurple,1,sum)
AdjMatHumanrestpurple=AdjMatHumanrest[colorh1=="purple",colorh1=="purple"]
ConnectivityHumanpurple=apply(AdjMatHumanrestpurple,1,sum)
ConnectivityChimppurple=ConnectivityChimppurple/max(ConnectivityChimppurple)
ConnectivityHumanpurple=ConnectivityHumanpurple/max(ConnectivityHumanpurple)

#6 Magenta
AdjMatChimprestmagenta=AdjMatChimprest[colorh1=="magenta",colorh1=="magenta"]
ConnectivityChimpmagenta=apply(AdjMatChimprestmagenta,1,sum)
AdjMatHumanrestmagenta=AdjMatHumanrest[colorh1=="magenta",colorh1=="magenta"]
ConnectivityHumanmagenta=apply(AdjMatHumanrestmagenta,1,sum)
ConnectivityChimpmagenta=ConnectivityChimpmagenta/max(ConnectivityChimpmagenta)
ConnectivityHumanmagenta=ConnectivityHumanmagenta/max(ConnectivityHumanmagenta)

#7 Pink
AdjMatChimprestpink=AdjMatChimprest[colorh1=="pink",colorh1=="pink"]
ConnectivityChimppink=apply(AdjMatChimprestpink,1,sum)
AdjMatHumanrestpink=AdjMatHumanrest[colorh1=="pink",colorh1=="pink"]
ConnectivityHumanpink=apply(AdjMatHumanrestpink,1,sum)
ConnectivityChimppink=ConnectivityChimppink/max(ConnectivityChimppink)
ConnectivityHumanpink=ConnectivityHumanpink/max(ConnectivityHumanpink)


#8 Red
AdjMatChimprestred=AdjMatChimprest[colorh1=="red",colorh1=="red"]
ConnectivityChimpred=apply(AdjMatChimprestred,1,sum)
AdjMatHumanrestred=AdjMatHumanrest[colorh1=="red",colorh1=="red"]
ConnectivityHumanred=apply(AdjMatHumanrestred,1,sum)
ConnectivityChimpred=ConnectivityChimpred/max(ConnectivityChimpred)
ConnectivityHumanred=ConnectivityHumanred/max(ConnectivityHumanred)
#9 Turquoise
AdjMatChimprestturquoise=AdjMatChimprest[colorh1=="turquoise",colorh1=="turquoise"]
ConnectivityChimpturquoise=apply(AdjMatChimprestturquoise,1,sum)
AdjMatHumanrestturquoise=AdjMatHumanrest[colorh1=="turquoise",colorh1=="turquoise"]
ConnectivityHumanturquoise=apply(AdjMatHumanrestturquoise,1,sum)
ConnectivityChimpturquoise=ConnectivityChimpturquoise/max(ConnectivityChimpturquoise)
ConnectivityHumanturquoise=ConnectivityHumanturquoise/max(ConnectivityHumanturquoise)
#10 Yellow
AdjMatChimprestyellow=AdjMatChimprest[colorh1=="yellow",colorh1=="yellow"]
ConnectivityChimpyellow=apply(AdjMatChimprestyellow,1,sum)
AdjMatHumanrestyellow=AdjMatHumanrest[colorh1=="yellow",colorh1=="yellow"]
ConnectivityHumanyellow=apply(AdjMatHumanrestyellow,1,sum)
ConnectivityChimpyellow=ConnectivityChimpyellow/max(ConnectivityChimpyellow)
ConnectivityHumanyellow=ConnectivityHumanyellow/max(ConnectivityHumanyellow)

# 11 Greenyellow
AdjMatChimprestgreenyellow=AdjMatChimprest[colorh1=="greenyellow",colorh1=="greenyellow"]
ConnectivityChimpgreenyellow=apply(AdjMatChimprestgreenyellow,1,sum)
AdjMatHumanrestgreenyellow=AdjMatHumanrest[colorh1=="greenyellow",colorh1=="greenyellow"]
ConnectivityHumangreenyellow=apply(AdjMatHumanrestgreenyellow,1,sum)
ConnectivityChimpgreenyellow=ConnectivityChimpgreenyellow/max(ConnectivityChimpgreenyellow)
ConnectivityHumangreenyellow=ConnectivityHumangreenyellow/max(ConnectivityHumangreenyellow)


# 12 Grey
AdjMatChimprestgrey=AdjMatChimprest[colorh1=="grey",colorh1=="grey"]
ConnectivityChimpgrey=apply(AdjMatChimprestgrey,1,sum)
AdjMatHumanrestgrey=AdjMatHumanrest[colorh1=="grey",colorh1=="grey"]
ConnectivityHumangrey=apply(AdjMatHumanrestgrey,1,sum)
ConnectivityChimpgrey=ConnectivityChimpgrey/max(ConnectivityChimpgrey)
ConnectivityHumangrey=ConnectivityHumangrey/max(ConnectivityHumangrey)



## To compare the correlation in kin between humans and chimpanzees by module: 
cor.test(ConnectivityHumanblack,ConnectivityChimpblack,method="s", exact = FALSE) 
cor.test(ConnectivityHumanblue,ConnectivityChimpblue,method="s", exact = FALSE) 
cor.test(ConnectivityHumanbrown,ConnectivityChimpbrown,method="s", exact = FALSE) 
cor.test(ConnectivityHumangreen,ConnectivityChimpgreen,method="s", exact = FALSE) 
cor.test(ConnectivityHumanred,ConnectivityChimpred,method="s", exact = FALSE) 
cor.test(ConnectivityHumanturquoise,ConnectivityChimpturquoise,method="s", exact = FALSE) 
cor.test(ConnectivityHumanyellow,ConnectivityChimpyellow,method="s", exact = FALSE)
cor.test(ConnectivityHumanpurple,ConnectivityChimppurple,method="s", exact = FALSE)

cor.test(ConnectivityHumanmagenta,ConnectivityChimpmagenta,method="s", exact = FALSE)
cor.test(ConnectivityHumanpink,ConnectivityChimppink,method="s", exact = FALSE)
cor.test(ConnectivityHumangreenyellow,ConnectivityChimpgreenyellow,method="s", exact = FALSE)
cor.test(ConnectivityHumangrey,ConnectivityChimpgrey,method="s", exact = FALSE)

## Plot magenta connectivity correlation:
par(mfrow=c(1,1))
plot(ConnectivityHumanmagenta,ConnectivityChimpmagenta,main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT magenta module connectivity", ylab="N3-/- magenta module connectivity", sub="rho =
0.72", cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ConnectivityChimpmagenta~ConnectivityHumanmagenta),col="red",lwd=2)

## Plot turquoise connectivity correlation:
par(mfrow=c(1,1))
plot(ConnectivityHumanturquoise,ConnectivityChimpturquoise,main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT turquoise module connectivity", ylab="N3-/- turquoise module connectivity", sub="rho =
0.72", cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ConnectivityChimpturquoise~ConnectivityHumanturquoise),col="red",lwd=2)

## Plot green connectivity correlation:
par(mfrow=c(1,1))
plot(ConnectivityHumangreen,ConnectivityChimpgreen, main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT green module connectivity", ylab="N3-/- green module connectivity", sub="rho =
0.56", cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ConnectivityChimpgreen~ConnectivityHumangreen),col="red",lwd=2)

## Plot black connectivity correlation:
par(mfrow=c(1,1))
plot(ConnectivityHumanblack,ConnectivityChimpblack, main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT black module connectivity", ylab="N3-/- black module connectivity", sub="rho =
0.22", cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ConnectivityChimpblack~ConnectivityHumanblack),col="red",lwd=2)

## Plot blue connectivity correlation:
par(mfrow=c(1,1))
plot(ConnectivityHumanblue,ConnectivityChimpblue, main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT blue module connectivity", ylab="N3-/- blue module connectivity", sub="rho =
0.62", cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ConnectivityChimpblue~ConnectivityHumanblue),col="red",lwd=2)

## Plot brown connectivity correlation:
par(mfrow=c(1,1))
plot(ConnectivityHumanbrown,ConnectivityChimpbrown, main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT brown module connectivity", ylab="N3-/- brown module connectivity", sub="rho =
0.69", cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ConnectivityChimpbrown~ConnectivityHumanbrown),col="red",lwd=2)

## Plot red connectivity correlation:
par(mfrow=c(1,1))
plot(ConnectivityHumanred,ConnectivityChimpred, main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT red module connectivity", ylab="N3-/- red module connectivity", sub="rho =
0.56", cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ConnectivityChimpred~ConnectivityHumanred),col="red",lwd=2)

## Plot yellow connectivity correlation:
par(mfrow=c(1,1))
plot(ConnectivityHumanyellow,ConnectivityChimpyellow, main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT yellow module connectivity", ylab="N3-/- yellow module connectivity", sub="rho =
0.57", cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ConnectivityChimpyellow~ConnectivityHumanyellow),col="red",lwd=2)

## Plot purple connectivity correlation:
par(mfrow=c(1,1))
plot(ConnectivityHumanpurple,ConnectivityChimppurple, main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT purple module connectivity", ylab="N3-/- purple module connectivity", sub="rho =
0.55", cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ConnectivityChimppurple~ConnectivityHumanpurple),col="red",lwd=2)

## Plot pink connectivity correlation:
par(mfrow=c(1,1))
plot(ConnectivityHumanpink,ConnectivityChimppink, main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT pink module connectivity", ylab="N3-/- pink module connectivity", sub="rho =
0.62", cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ConnectivityChimppink~ConnectivityHumanpink),col="red",lwd=2)

## Plot greenyellow connectivity correlation:
par(mfrow=c(1,1))
plot(ConnectivityHumangreenyellow,ConnectivityChimpgreenyellow, main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT greenyellow module connectivity", ylab="N3-/- greenyellow module connectivity", sub="rho =
0.52", cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ConnectivityChimpgreenyellow~ConnectivityHumangreenyellow),col="red",lwd=2)

## Plot grey connectivity correlation:
par(mfrow=c(1,1))
plot(ConnectivityHumangrey,ConnectivityChimpgrey, main="Age pup - adult, region CA1, DG, CA3",
     xlab="WT grey module connectivity", ylab="N3-/- grey module connectivity", sub="rho =
0.08", cex.main=2, cex.lab=2, cex.axis=1.5, cex.sub=1.25)
abline(lm(ConnectivityChimpgrey~ConnectivityHumangrey),col="red",lwd=2)



## To calculate the topological overlap for all genes within each module for human and chimp: # Yellow
distTOMHumanyellow <- TOMdist1(AdjMatHumanrestyellow)
simTOMHumanyellow = 1-distTOMHumanyellow
diag(simTOMHumanyellow)=0
distTOMChimpyellow <- TOMdist1(AdjMatChimprestyellow)
simTOMChimpyellow = 1-distTOMChimpyellow
diag(simTOMChimpyellow)=0
## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# Yellow
simRatioyellow=simTOMHumanyellow/mean(simTOMHumanyellow)/(simTOMHumanyellow/mean(simTOMHumanyellow)+ 
                                                               simTOMChimpyellow/mean(simTOMChimpyellow)+0.00001)
# Comment: we add the negligible constant of 0.00001 to the denominator as a pre-caution to #ensure that the denominator is non-zero.
corrplot(simRatioyellow, method = "circle", col.lim = c(0,0.5))
write.csv(simRatioyellow, file = "simRatioyellow.csv")
bigratio_yellow <- simRatioyellow > 0.65
bigratio_yellow
simRatioyellow
dim(simRatioyellow)
write.csv(bigratio_yellow, file = "bigratio_yellow.csv")



## To calculate the topological overlap for all genes within each module for human and chimp: # Magenta
distTOMHumanmagenta <- TOMdist1(AdjMatHumanrestmagenta)
simTOMHumanmagenta = 1-distTOMHumanmagenta
diag(simTOMHumanmagenta)=0
distTOMChimpmagenta <- TOMdist1(AdjMatChimprestmagenta)
simTOMChimpmagenta = 1-distTOMChimpmagenta
diag(simTOMChimpmagenta)=0
## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# Magenta 
simRatiomagenta=simTOMHumanmagenta/mean(simTOMHumanmagenta)/(simTOMHumanmagenta/mean(simTOMHumanmagenta)+ 
                                                         simTOMChimpmagenta/mean(simTOMChimpmagenta)+0.00001)


## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# Brown
distTOMHumanbrown <- TOMdist1(AdjMatHumanrestbrown)
simTOMHumanbrown = 1-distTOMHumanbrown
diag(simTOMHumanbrown)=0
distTOMChimpbrown <- TOMdist1(AdjMatChimprestbrown)
simTOMChimpbrown = 1-distTOMChimpbrown
diag(simTOMChimpbrown)=0

simRatiobrown=simTOMHumanbrown/mean(simTOMHumanbrown)/(simTOMHumanbrown/mean(simTOMHumanbrown)+ 
                                                               simTOMChimpbrown/mean(simTOMChimpbrown)+0.00001)


# Comment: we add the negligible constant of 0.00001 to the denominator as a pre-caution to #ensure that the denominator is non-zero.
install.packages("corrplot")
library(corrplot)
corrplot(simRatiomagenta, method = "circle", col.lim = c(0,1))
write.csv(simRatiomagenta, file = "simRatiomagenta.csv")
bigratio_magenta <- simRatiomagenta > 0.65
bigratio_magenta
simRatiomagenta
write.csv(simRatiomagenta, file = "simRatiomagenta.csv")
write.csv(bigratio_magenta, file = "bigratio_magenta.csv")



## To calculate the topological overlap for all genes within each module for human and chimp: # Black
distTOMHumanblack <- TOMdist1(AdjMatHumanrestblack)
simTOMHumanblack = 1-distTOMHumanblack
diag(simTOMHumanblack)=0
distTOMChimpblack <- TOMdist1(AdjMatChimprestblack)
simTOMChimpblack = 1-distTOMChimpblack
diag(simTOMChimpblack)=0
## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# Black
simRatioblack=simTOMHumanblack/mean(simTOMHumanblack)/(simTOMHumanblack/mean(simTOMHumanblack)+ 
              simTOMChimpblack/mean(simTOMChimpblack)+0.00001)
# Comment: we add the negligible constant of 0.00001 to the denominator as a pre-caution to #ensure that the denominator is non-zero.
install.packages("corrplot")
library(corrplot)

###### BLACK ######################
corrplot(simRatioblack, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "blue", "yellow"))(200))


bigratio_black <- simRatioblack > 0.65
write.csv(simRatioblack, file = "simRatioblack.csv")
write.csv(bigratio_black, file = "bigRatioblack.csv")
bigratio_black

### we only take the ratios bigger than 0.65 and plot them as "1" and "0"
dat <- bigratio_black
dat <- 1*dat
head(dat)
corrplot(dat, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "dark grey","yellow"))(200))

lowratio_black <- simRatioblack < 0.35
low_black <- lowratio_black
low_black <- 1*low_black
head(low_black)
corrplot(low_black, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "dark grey","blue"))(200))


####### MAGENTA ##################
simRatiomagenta=simTOMHumanmagenta/mean(simTOMHumanmagenta)/(simTOMHumanmagenta/mean(simTOMHumanmagenta)+ 
                                                         simTOMChimpmagenta/mean(simTOMChimpmagenta)+0.00001)
# Comment: we add the negligible constant of 0.00001 to the denominator as a pre-caution to #ensure that the denominator is non-zero.
corrplot(simRatiomagenta, method = "color", col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC", 
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "blue","yellow"))(200))

### we only take the ratios bigger than 0.65 and plot them as "1" and "0"
bigratio_magenta <- simRatiomagenta > 0.65
dat_magenta <- bigratio_magenta
dat_magenta <- 1*dat_magenta
head(dat_magenta)
corrplot(dat_magenta, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "dark grey","yellow"))(200))

lowratio_magenta <- simRatiomagenta < 0.35
low_magenta <- lowratio_magenta
low_magenta <- 1*low_magenta
head(low_magenta)
corrplot(low_magenta, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC",
         tl.cex = 0.6,
         col=colorRampPalette(c("white", "dark grey","blue"))(200))

####### YELLOW ##################
simRatioyellow=simTOMHumanyellow/mean(simTOMHumanyellow)/(simTOMHumanyellow/mean(simTOMHumanyellow)+ 
                                                               simTOMChimpyellow/mean(simTOMChimpyellow)+0.00001)
# Comment: we add the negligible constant of 0.00001 to the denominator as a pre-caution to #ensure that the denominator is non-zero.
corrplot(simRatioyellow, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC", 
         tl.cex = 0.2, 
         col=colorRampPalette(c("white", "blue","yellow"))(200))

### we only take the ratios bigger than 0.65 and plot them as "1" and "0"
bigratio_yellow <- simRatioyellow > 0.65
dat_yellow <- bigratio_yellow
dat_yellow <- 1*dat_yellow
head(dat_yellow)

corrplot(dat_yellow, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC", 
         tl.cex = 0.2, 
         col=colorRampPalette(c("white", "dark grey","yellow"))(200))

### we only take the ratios lower than 0.3 and plot them as "1" and "0"
lowratio_yellow <- simRatioyellow < 0.3
low_yellow <- lowratio_yellow
low_yellow <- 1*low_yellow
head(low_yellow)
write.csv(lowratio_yellow, file = "lowRatioyellow.csv")

corrplot(low_yellow, method = "color", 
         col.lim = c(0,1), 
         diag = FALSE, 
         order = "FPC", 
         tl.cex = 0.2, 
         col=colorRampPalette(c("white", "dark grey","blue"))(200))

## To calculate the topological overlap for all genes within each module for human and chimp: # Turquoise

############ TURQOUISE ################

distTOMHumanturquoise <- TOMdist1(AdjMatHumanrestturquoise)
simTOMHumanturquoise = 1-distTOMHumanturquoise
diag(simTOMHumanturquoise)=0
distTOMChimpturquoise <- TOMdist1(AdjMatChimprestturquoise)
simTOMChimpturquoise = 1-distTOMChimpturquoise
diag(simTOMChimpturquoise)=0
## Creating a matrix for each module to describe the relative topological overlap for all #connections between human and chimp:
# Black
simRatioturquoise=simTOMHumanturquoise/mean(simTOMHumanturquoise)/(simTOMHumanturquoise/mean(simTOMHumanturquoise)+ 
                                                         simTOMChimpturquoise/mean(simTOMChimpturquoise)+0.00001)
# Comment: we add the negligible constant of 0.00001 to the denominator as a pre-caution to #ensure that the denominator is non-zero.








###### EXTRACT GENES FROM MODULES ##########################

black <- names(datExprHumanrest)[colorh1=="black"]
turquoise <- names(datExprHumanrest)[colorh1=="turquoise"]
brown <- names(datExprHumanrest)[colorh1=="brown"]
blue <- names(datExprHumanrest)[colorh1=="blue"]
greenyellow <- names(datExprHumanrest)[colorh1=="greenyellow"]
green <- names(datExprHumanrest)[colorh1=="green"]
yellow <- names(datExprHumanrest)[colorh1=="yellow"]
grey <- names(datExprHumanrest)[colorh1=="grey"]
magenta <- names(datExprHumanrest)[colorh1=="magenta"]
purple <- names(datExprHumanrest)[colorh1=="purple"]
pink <- names(datExprHumanrest)[colorh1=="pink"]
red <- names(datExprHumanrest)[colorh1=="red"]


write.table(black, file = "WT_black.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(turquoise, file = "WT_turquoise.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(brown, file = "WT_brown.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(blue, file = "WT_blue.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(greenyellow, file = "WT_greenyellow.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(green, file = "WT_green.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(yellow, file = "WT_yellow.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(grey, file = "WT_grey.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(magenta, file = "WT_magenta.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(purple, file = "WT_purple.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(pink, file = "WT_pink.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(red, file = "WT_red.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)


blue_cutoff_0.62 <- names(datExprHumanrest)[colorh12=="blue"]
write.table(blue_cutoff_0.62, file = "WT_blue_cutoff_0.62.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

turquoise_cutoff_0.62 <- names(datExprHumanrest)[colorh12=="turquoise"]
write.table(turquoise_cutoff_0.62, file = "WT_turquoise_cutoff_0.62.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

library(tidyverse)
library(ISLR)

# load data object
data(NCI60)

# see structure of object
str(NCI60)

# cancer type
nci.labs <- NCI60$labs

# gene expression data for 6830 
nci.data <- NCI60$data

# check dimension of data set
dim(nci.data)

# check variable names
names(nci.data)

# summarize labels
table(nci.labs)

#######################################
# PCA on Gene Expression Data
#######################################

# look at variance of each gene expression variable
hist(apply(nci.data, 2, var))

# conduct PCA on scaled gene expression data
pr.out <- prcomp(nci.data, scale = TRUE)

# structure of prcomp object
str(pr.out)

# standard deviations for the PCs (square root of lambdas)
pr.out$sdev

# PC loading vectors in columns
pr.out$rotation # shouldn't there only be min(64 -1, 6830) = 63? Yes, but this uses SVD to compute PCs and notice that the eigenvalue for the last PC is 0.

# PC scores for the 64 cells
pr.out$x

# plot first few PCs
Cols <- function(vec) {
 cols <- rainbow(length(unique(vec)))
 return(cols[as.numeric(as.factor(vec))])
}

# plot Z1 vs Z2 and Z1 vs. Z3
#  colors correspond to type of cancer
par(mfrow = c(1, 2))
plot(pr.out$x[, 1:2], col = Cols(nci.labs), pch = 19,
     xlab = "Z1", ylab = "Z2")
plot(pr.out$x[, c(1, 3)], col = Cols(nci.labs), pch = 19,
     xlab = "Z1", ylab = "Z3")

# look at proportion of variance explained by PCs
summary(pr.out)

# plot percent variance explained and cumulative percent variance explained 
pve <- 100 * pr.out$sdev^2 / sum(pr.out$sdev^2)
par(mfrow = c(1, 2))
plot(pve, type = "o", ylab = "PVE",
     xlab = "Principal Component", col = "blue")
plot(cumsum(pve), type = "o", ylab = "Cumulative PVE",
     xlab = "Principal Component", col = "brown3")

#######################################
# Clustering Gene Expression Data
#######################################

# scale the columns of the data matrix to each have variance 1
sd.data <- scale(nci.data)

# do heierarchical clustering on observations
#  using Euclidean distance as the dissimilarity measure

# compute dissimilarity matrix
data.dist <- dist(sd.data) 

# perform clustering
hc.compl <- hclust(data.dist)
hc.avg <- hclust(data.dist, method = "average")
hc.sing <- hclust(data.dist, method = "single")

# look at hc.compl
hc.compl

par(mfrow = c(1, 3))
plot(hc.compl, xlab = "", sub = "", ylab = "",
     labels = nci.labs, main = "Complete Linkage")
plot(hc.avg, labels = nci.labs, main = "Average Linkage",
     xlab = "", sub = "", ylab = "")
plot(hc.sing, labels = nci.labs, main = "Single Linkage", 
     xlab = "", sub = "", ylab = "")

# suppose that we are interested in looking at the 
#   4-cluster solution for complete linkage 
hc.clusters <- cutree(hc.compl, 4)

# see how clusters align with cancer type
table(hc.clusters, nci.labs)

# plot dendrogram with horizontal line that cuts into 4 clusters
par(mfrow = c(1, 1))
plot(hc.compl, labels = nci.labs)
abline(h = 143, col = "red")

# try k-means clustering with 4 clusters
set.seed(2) # need to set seed because of random starts
km.out <- kmeans(sd.data, 4, nstart = 20)
km.clusters <- km.out$cluster

# compare k-means clusters with hierarchical clustering (complete linkage)
table(km.clusters, hc.clusters)

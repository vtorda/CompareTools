#install.packages("fastcluster")
library(fastcluster)
h <- hclust.vector(t(probe1$event_data$gains), method = "ward", metric = "euclidean")
plot(h, labels = FALSE)
rect.hclust(h, k = 25, border = "red")

#lehet hogy manhattan distancet kene szamolni?

# https://bib.dbvis.de/uploadedFiles/155.pdf


h_man <- hclust.vector(t(probe1$event_data$gains), method = "single", metric = "manhattan")
plot(h_man, labels = FALSE)
h_can <- hclust.vector(t(probe1$event_data$gains), method = "single", metric = "canberra")
plot(h_can, labels = FALSE)
#### de amugy teljesen ugy nez ki hogy az euklideszi tavolsag ad empirikusan jo klusterezest


# ki kene probalni pearson correlaciot
#install.packages("devtools")
#library(devtools)
#install_github("anspiess/propagate")
################
###sajnos nem megy, mert hclust kiakasztja a memoriat....szoval ugy nez ki hogy egyenlore csak euklideaszi tavolsaggal fogok szamolni..
###########
###lent ne futtasd!!!!!!!!!!!
library(propagate)
dm <- bigcor(probe1$event_data$gains, fun = "cor", method = "pearson", size = 5000)
BCOR <- dm[1:ncol(probe1$event_data$gains), 1:ncol(dm)]
#h2 <- hclust(as.dist(1-BCOR), method = "ward.D2")
rm(BCOR)
rm(dm)
MAT <- matrix(rnorm(70000), ncol = 700)
colnames(MAT) <- str_c("colname_", 1:700)
rownames(MAT) <- str_c("spec_", 1:100)
COR <- bigcor(MAT, size= 500, fun = "cor")
COR <- COR[1:nrow(COR), 1:ncol(COR)]
COR2 <- cor(MAT)
all.equal(COR, cor(MAT)) # => TRUE

COR[,1]
COR2[,1]
ncol(MAT1)
sapply(1:ncol(COR), function(x) all.equal(COR[,x], unname(COR2[,x])))
plot(COR[,1], COR2[,1])
BCOR <- BCOR[1:5000, 1:ncol(BCOR)]

# meg azt meglehetne csinalni hogy pearson alapjan egy clusterezest pl vagy 10ezer elemre es utana ezekhez clusterezni a maradek pontot.....

# elsokent akkor mintazni kene es abbol szamolni
Sample5K = sample(ncol(probe1$event_data$gains), 5000)
m <- probe1$event_data$gains[,Sample5K]
colnames(m)
dm <- bigcor(m, fun = "cor", method = "pearson", size = 1000)
BCOR <- dm[1:ncol(m), 1:ncol(dm)]
colnames(BCOR) <- colnames(m)
h2 <- fastcluster::hclust(as.dist(1-BCOR), method = "ward.D2")
plot(h2, labels = FALSE)
rect.hclust(h2, k = 50, border = "red")
groups <- cutree(h2, k = 50)
library(class)
#
# x = c(rnorm(250000, 0,0.9), rnorm(350000, 4,1), rnorm(500000, -5,1.1))
# y = c(rnorm(250000, 0,0.9), rnorm(350000, 5.5,1), rnorm(500000,  5,1.1))
# XY = data.frame(x,y)
# Sample5K = sample(length(x), 5000)     ## Downsample
#
# ## Cluster the sample
# DM5K = dist(XY[Sample5K,])
# HC5K = hclust(DM5K, method="single")
# Groups = cutree(HC5K, 8)
# Groups[Groups>4] = 4
# unique(Groups)
# Core = which(Groups<4)
#
# knnClust = knn(XY[Sample5K[Core], ], XY, Groups[Core])

length(Groups[Core])
nrow(XY[Sample5K[Core], ])
knnClust <- knn(t(m), t(probe1$event_data$gains), groups)
length(knnClust)
knnClust[1:10]


phytools::phyl.pairedttest()

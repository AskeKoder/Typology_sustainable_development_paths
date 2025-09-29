#Sensitivity to Imputation method
#Load libraries-----------------------
library(dplyr)
library(NbClust)
library(partitionComparison)

#Functions for clustering --------------------------------------------
#Dissimilarity function
EuclideanDist <- function (X,Y) {
  #X and Y are multivariate time series with a column for each variable and time down the rows
  #X and Y have equal dimensions and no missing data
  sqrt(sum((X-Y)^2 ))
}

#Define clustering function
clusterData <- function(data,m){
  #Set up matrices to store results
  clusteringsPC <- matrix(0,nrow=length(unique(data$Country)),
                          ncol=m)
  rownames(clusteringsPC) <-unique(data$Country)
  nPCs <- rep(0,m)
  nclust <- rep(0,m)
  
  #Cluster each of the m sets
  for (set in 1:m){
    print(paste("Processing set:",set))
    #Extract set
    complete_set <- data[data$imp==set,]%>%
      select(-imp)
    
    #Keep only the indicator columns
    complete_set_numeric <- complete_set[, sapply(complete_set, is.numeric)]
    complete_set_numeric <- complete_set_numeric[,2:ncol(complete_set_numeric)]
    
    #Run PCA
    pca_result <- prcomp(complete_set_numeric, center = TRUE, scale. = TRUE)
    var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
    cumu_var_explained <- cumsum(var_explained)
    
    #Get number of PCs that explain >80% of variance
    nPC <- min(which(cumu_var_explained>0.8))
    nPCs[set] <- nPC
    
    #Select PCs
    print(paste("Principal components used:",nPC))
    reduced_set <- cbind(complete_set[,1:4],pca_result$x[,1:nPC]) 
    
    #Allocate space for distance matrix
    distMat <- matrix(0,nrow=length(unique(reduced_set$Country)),
                      ncol=length(unique(reduced_set$Country)))
    rownames(distMat) <- unique(reduced_set$Country)
    colnames(distMat) <- unique(reduced_set$Country)
    
    #Precompute indices for loops
    countries <- unique(reduced_set$Country)
    rows_by_country <- lapply(countries, function(cty) {
      which(reduced_set$Country == cty)
    })
    nCountries <- length(countries)
    
    #Calculate the distances based on numeric variables (excluding year)
    for (i in 1:nCountries){
      dfi <- reduced_set[rows_by_country[[i]], 5:ncol(reduced_set)]
      for (j in i:nCountries){
        dfj <- reduced_set[rows_by_country[[j]], 5:ncol(reduced_set)]
        distMat[i,j] <- EuclideanDist(dfi,dfj)
      }
    }
    #Put values below the diagonal
    distMat <- distMat + t(distMat)
    #Find optimal number of clusters in set
    nb <- NbClust(data = NULL, diss = as.dist(distMat), distance = NULL,
                  min.nc = 2, max.nc =40, #max number of cluster set to 40 to improve speed
                  index = "cindex",
                  method="ward.D2")
    
    
    nclust[set] <- nb$Best.nc[1] 
    optClustering <- nb$Best.partition
    clusteringsPC[,set] <- optClustering
    print(paste("Set",set, "finished with",nb$Best.nc[1],"clusters." ))
  }
  
  #Median partitioning by BOK (best of k)
  partitionQualityPC <- rep(0,m)
  
  for (partition in 1:m){
    #Compute values for calculating dissimilarity
    pc <- sapply(1:m,function(x){computePairCoefficients(clusteringsPC[,partition],clusteringsPC[,x])})
    N10 <- sapply(1:m,function(x){pc[[x]]}@N10) #number of pairs in 1 but not in 2
    N01 <- sapply(1:m,function(x){pc[[x]]}@N01) #number of pairs in 2 but not in 1
    
    #Compute Mirkin distance fpr "partition"
    partitionQualityPC[partition] <- sum(N10+N01)
  }
  
  #Final Clustering
  finalClusteringPC <- clusteringsPC[,which.min(partitionQualityPC)]
  return(list(
    finalClustering = finalClusteringPC,
    nclust = nclust,
    nPCs = nPCs
  ))
}
#Import different imputed datasets -------------------------------------
data_lag1_lead1 <- read.csv("imputedDataLag1Lead1.csv")%>%
  select(-c(X,.id))%>%
  relocate(.imp, .after = last_col())%>%
  rename(imp = .imp)


data_lag1_lead2 <- read.csv("imputedDataLag1Lead2.csv")%>%
  select(-c(X,.id))%>%
  relocate(.imp, .after = last_col())%>%
  rename(imp = .imp)

data_Lighthouse <- read.csv("imputedDataLightHouse.csv")%>%
  select(-c(X,.id))%>%
  relocate(.imp, .after = last_col())%>%
  rename(imp = .imp)

#Set up tests
testsets <- list(data_lag1_lead1,
                 data_lag1_lead2,
                 data_Lighthouse)

results <- list()
for (i in 1:3){
  results[[i]] <- clusterData(testsets[[i]],m=15)
}




#Analysis
results[[1]]$finalClustering- results[[2]]$finalClustering


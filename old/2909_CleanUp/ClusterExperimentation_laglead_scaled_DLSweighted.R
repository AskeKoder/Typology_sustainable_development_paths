#Load libraries
library(readxl)
library(dplyr)
library(tidyr)
library(partitionComparison)
library(NbClust)
library(ggplot2)
library(reshape2)

#Read imputed data 
data <- read.csv("ImputedDataLag1Lead2_maxit30_scaled.csv")%>%
  select(-c(X,.id))%>%
  relocate(.imp, .after=last_col())

#Read experiment file
Experiments <- read_xlsx("IndicatorsForClusters.xlsx",sheet="Experiments")

#Read raw data for reference
raw_data <- read.csv("extendedDataScaled.csv")%>%
  select(-X)

#Read colonial data for filtering
df_colonial <- read.csv("Colonial history/dfColonial.csv")

#Functions ----------------------
#Dissimilarity function
EuclideanDist <- function (X,Y) {
  #X and Y are multivariate time series with a column for each variable and time down the rows
  #X and Y have equal dimensions and no missing data
  sqrt(sum((X-Y)^2 ))
}

#Define clustering function
clusterData <- function(data,weights=FALSE,m,minInformation,max.nc=40){
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
    
    #Run weighted PCA
    if (length(weights)==1){
    pca_result <- prcomp(complete_set_numeric, center = TRUE, scale. = TRUE)
    var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
    cumu_var_explained <- cumsum(var_explained)
    }else{
      #Center, Scale and apply weights
      complete_set_weighted <- complete_set_numeric%>%
        mutate_all(~scale(.))
      complete_set_weighted <- complete_set_weighted* weights[col(complete_set_weighted)]
      
      #run PCA
      pca_result <- prcomp(complete_set_weighted, center = FALSE, scale. = FALSE)
      var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
      cumu_var_explained <- cumsum(var_explained)
    }
    #Get number of PCs that explain >80% of variance
    nPC <- min(which(cumu_var_explained>minInformation))
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
                  min.nc = 2, max.nc =max.nc, #max number of cluster set to 40 to improve speed
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

#Run experiments selected experiments -----------------------------------------
batch_weighted<- Experiments%>%
  select(c(1,2, "Baseline",
           "BHN_FWB",
           "No_non_DLS",
           "No_non_DLS_BHNFWB"
  ))%>%
  as.matrix()%>%
  na.omit()

#Format names as in data
rownames(batch_weighted) <- colnames(data)[5:(ncol(data)-1)]
rownames(batch_weighted) <- rownames(batch_weighted)[c(6,5,4,3,2,1,
                                                   10,9,8,7,
                                                   11,12,13,
                                                   18,17,16,15,14,
                                                   22,21,20,19,
                                                   23,24,25,
                                                   29,28,27,26,
                                                   33,32,31,30,
                                                   39,38,37,36,35,34,
                                                   44,43,42,41,40,
                                                   48,47,46,45,
                                                   52,51,50,49,
                                                   53,55,57,56,
                                                   54,60,58,59,
                                                   61,62,63)]


#Remove duplicate column
batch_weighted <- batch_weighted%>%
  data.frame()%>%
  select(-c(1))%>%
  mutate(across(2:5,~as.integer(.)))


#Calculate weights 
weights <- batch_weighted%>%
  group_by(Related.DLS.dimension)%>%
  summarize_all(~1/sqrt(sum(.)))%>%
  mutate(across(2:5,~as.numeric(.)))
weights[weights == Inf]  <- 0 #Set weights to zero for categories not included


#Define number of experiments to be conducted
results_batch_weighted <- list()

for (i in 2:ncol(batch_weighted)){
  #Get data for experiment
  id <- which(batch_weighted[,i]>0)
  indicators<- rownames(batch_weighted)[id]
  
  
  #Make dataset for experiment
  expData <- data%>%
    select(c("Country","SPI_countrycode","SPI_year","Region",
             indicators,
             ".imp"))%>%
    rename(imp=.imp)
  
  expWeights <- left_join(batch_weighted[indicators,],weights, by="Related.DLS.dimension")[,4+i]
  #Run test and save results
  results_batch_weighted[[i]] <- clusterData(expData, weights = expWeights, m=15,minInformation = 0.8)
}


results_batch_weighted[[5]]$finalClustering


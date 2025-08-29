#Sensitivity analysis
#Libraries ======================
library(NbClust)
library(dplyr)
library(partitionComparison) #For consensus function
library(ggplot2) 
library(reshape2)

#Load imputed data
data <- complete_long%>% #Load from csv
  select(-c(68:76)) #Remove IEA variables

#Number of indicators
nVar <- ncol(data%>%select(-c(imp,Country_0,SPI_countrycode_1,SPI_year_2,Region_3)))
#Number of years
nYear <-length(unique(data$SPI_year_2))
#Number of Countries
nCountry <- length(unique(data$Country_0))

#Functions ======================

#Dissimilarity function
EuclideanDist <- function (X,Y) {
  #X and Y are multivariate time series with a column for each variable and time down the rows
  #X and Y have equal dimensions and no missing data
  sqrt(sum((X-Y)^2 ))
}

#Define clustering function
clusterData <- function(data,m){
  #Set up matrices to store results
  clusteringsPC <- matrix(0,nrow=length(unique(data$Country_0)),
                          ncol=m)
  rownames(clusteringsPC) <-unique(data$Country_0)
  nPCs <- rep(0,m)
  nclust <- rep(0,m)
  
  #Cluster each of the m sets
  for (set in 1:m){
    print(paste("Processing set:",set))
    #Extract set
    complete_set <- data[data$imp==set,]%>%
      select(-imp)
    ##Do PCA only on SPI indicators
    # complete_set <- complete_set[c(1:56)]
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
    distMat <- matrix(0,nrow=length(unique(reduced_set$Country_0)),
                      ncol=length(unique(reduced_set$Country_0)))
    rownames(distMat) <- unique(reduced_set$Country_0)
    colnames(distMat) <- unique(reduced_set$Country_0)
    
    #Precompute indices for loops
    countries <- unique(reduced_set$Country_0)
    rows_by_country <- lapply(countries, function(cty) {
      which(reduced_set$Country_0 == cty)
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



#Sensitivity to increasingly missing variables ======================
#Sort columns by how much is missing
sortedData <- data[,c(1:56,order(colSums(is.na(data[,57:ncol(data)])))+56)]

#Generate matrix of variables to include in each run
testingSequence <- matrix(0,ncol(sortedData),ncol(sortedData))
rownames(testingSequence) <- colnames(sortedData)
testingSequence[upper.tri(testingSequence, diag = TRUE)] <- 1 

#Include only 1 test with all variables with no missing data
id <- 57 
testingSequence <- testingSequence[,id:ncol(sortedData)]

#Calculate fraction of missing data in each test
fMis = matrix(0,nrow=ncol(testingSequence),ncol=1)
for (i in 1:ncol(testingSequence)){
  totNA <- sum(colSums(is.na(sortedData))[which(testingSequence[,i]>=1)]) #Number of NAs per test
  totVars <- sum(testingSequence[,i]>=1) - 5 #number of indicators - identifiers
  fMis[i] <- totNA/(totVars*nYear*nCountry)
}
missingPercent <- data.frame("run"=1:12,
                             "mis"=fMis)

#Run test
testData <- sortedData[,which(testingSequence[,1]>0)]
test <- clusterData(testData,m=2)

#Set up sensitivity analysis
sensitivityAnalysis <- function(sortedData,testingSequence){
  mtest <- 35
  Nclust <- matrix(0,ncol=ncol(testingSequence), nrow=mtest)
  Clusterings <- matrix(0,ncol=ncol(testingSequence), nrow=nCountry)
  NPCs <- Nclust <- matrix(0,ncol=ncol(testingSequence), nrow=mtest)
  
  for (i in 1:ncol(testingSequence)){
    #Select variables
    testData <- sortedData[,which(testingSequence[,i]>0)]
    testResults <- clusterData(testData,m=mtest)
    Nclust[,i] <- testResults$nclust
    NPCs[,i] <- testResults$nPCs
    Clusterings[,i] <- testResults$finalClustering
    
  }
  return(list(Nclust=Nclust,
              NPCs=NPCs,
              Clusterings=Clusterings))
}

#Calculate an save results
#results <- sensitivityAnalysis(sortedData,testingSequence)
#saveRDS(results, "SensMoreVars.rds") 
results <- readRDS("SensMoreVars.rds")

#Post processing
# Add number of clusters in final clustering
colMax <- function(data) sapply(data.frame(data), max, na.rm=TRUE)
optNum <- colMax(results$Clusterings)
optDf <- data.frame("run"=1:12,
                    "optNum"=optNum)

#Plotting results of sensitivity analysis
df <- melt(results$Nclust)
colnames(df) <- c("Run", "A", "Value")   # better names

ggplot(df, aes(x = factor(A), y = Value)) +
  geom_point(aes(color = "Clusterings"), alpha = 0.1,size=2.5) +
  geom_line(data = optDf, aes(x = run, y = optNum, color = "Final clustering"), size = 1) +
  # Add second dataset, scaled to match first axis
  geom_line(data = missingPercent, aes(x = run, y = mis*100, color = "% Missing values"), size = 1) +
  scale_y_continuous(
    name = "Number of clusters",
    sec.axis = sec_axis(~ ., name = "% Missing values")  # inverse transform
  ) +
  scale_color_manual(values = c("Clusterings" = "black", "Final clustering" = "red", "% Missing values" = "blue")) +
  labs(x = "A", title = "", color = "") +
  theme_minimal()


#Sensitivity to wrong imputations =============================
#Questionable indicators
# Access to electricity
# Satisfied demand for contraception
# Species protection

#Select SPI data
SPIdata <- data %>%
  select(1:56,68)

#Record positions of missing values in the variables to test
NApos <- data.frame(is.na(SPIdata[SPIdata$imp == 0,]))
NApos <- NApos %>%
  mutate(across(-c(Access_to_electricity_16,Species_protection_33,Satisfied_demand_for_contraception_45), ~ FALSE))

#Set up list of experiments
experiments <- data.frame("Indicator"=c(rep("Access_to_electricity_16",2),
                                        rep("Species_protection_33",2),
                                        rep("Satisfied_demand_for_contraception_45",2)),
                          "Scalar"=c(-15,-30,
                                     10,-10,
                                     10,-10))

#Add scalars to generate new sets, where imputations are lower/higher


dataList <- list(rep(SPIdata,nrow(experiments)))
for (i in 1:nrow(experiments)){
  #Get experiment conditions
  ind <- experiments[i,"Indicator"]
  scalar <- experiments[i,"Scalar"]
  
  #Reset data
  dataSens <- SPIdata
  for (j in 1:max(data$imp)){
    #Get imputed set and indicator
    dataToChange <- dataSens[dataSens$imp==j,ind]
    #Alter imputed values
    dataToChange[which(NApos[,ind]),] <- dataToChange[which(NApos[,ind]),] + scalar
    #Put back into place
    dataSens[dataSens$imp==j,ind] <- dataToChange
  }
  #Ensure that values are still limited between 0 and 100
  dataSens[which(dataSens[,ind] < 0),ind] <- 0
  dataSens[which(dataSens[,ind] > 100),ind] <- 100
  
  #Save experiment data in list
  dataList[[i]] <- dataSens
}

#Assess new data
testdat <- data-dataList[[3]]
sum(testdat[,"Species_protection_33"],na.rm=TRUE)

#Run Experiments
mtest <- 35
Nclust <- matrix(0,ncol=nrow(experiments), nrow=mtest)
Clusterings <- matrix(0,ncol=nrow(experiments), nrow=nCountry)
NPCs <- Nclust <- matrix(0,ncol=nrow(experiments), nrow=mtest)

for (i in 1:nrow(experiments)){
  expResults <- clusterData(dataList[[i]], m=mtest)
  Nclust[,i] <- expResults$nclust
  NPCs[,i] <- expResults$nPCs
  Clusterings[,i] <- expResults$finalClustering
}

resultsAnalysis2 <- list(Nclust=Nclust,
                         NPCs=NPCs,
                         Clusterings=Clusterings)


#saveRDS(resultsAnalysis2, "SensImputations.rds") 
resultsAnalysis2<- readRDS("SensImputations.rds")

df <- melt(resultsAnalysis2$Nclust)
colnames(df) <- c("Run", "A", "Value")   # better names

optNum <- colMax(resultsAnalysis2$Clusterings)
optDf <- data.frame("run"=1:nrow(experiments),
                    "optNum"=optNum)

ggplot(df, aes(x = factor(A), y = Value)) +
  geom_point(aes(color = "Clusterings"), alpha = 0.1,size=2.5) +
  geom_line(data = optDf, aes(x = run, y = optNum, color = "Final clustering"), size = 1) +
  # Add second dataset, scaled to match first axis
  scale_y_continuous(
    name = "Number of clusters"
  ) +
  scale_color_manual(values = c("Clusterings" = "black", "Final clustering" = "red", "% Missing values" = "blue")) +
  labs(x = "A", title = "", color = "") +
  theme_minimal()


check <- data-dataList[[1]]

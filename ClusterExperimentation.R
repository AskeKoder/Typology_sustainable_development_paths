#Load libraries
library(readxl)
library(dplyr)
library(partitionComparison)
library(NbClust)
library(ggplot2)
library(reshape2)

#Read imputed data 
data <- read.csv("ImputedDataLag1Lead2_maxit30.csv")%>%
  select(-c(X,.id))%>%
  relocate(.imp, .after=last_col())

#Read experiment file
Experiments <- read_xlsx("IndicatorsForClusters.xlsx",sheet="Experiments")

#Read raw data for reference
raw_data <- read.csv("extendedData.csv")%>%
  select(-X)


#Functions ----------------------
#Dissimilarity function
EuclideanDist <- function (X,Y) {
  #X and Y are multivariate time series with a column for each variable and time down the rows
  #X and Y have equal dimensions and no missing data
  sqrt(sum((X-Y)^2 ))
}

#Define clustering function
clusterData <- function(data,m,minInformation){
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

#Run experiments-----------------------------------------
#Get batch of experiments
batch1 <- Experiments%>%
  select(1,3,4,5,6,7,8)%>%
  as.matrix()%>%
  na.omit() # clear row of unrepresented dimension


#Format names as in data
rownames(batch1) <- colnames(data)[5:(ncol(data)-1)]

#Remove duplicate column
batch1 <- batch1%>%
  data.frame()%>%
  select(-c(1))%>%
  mutate_all(~as.integer(.))

#Define number of experiments to be conducted
results_batch1 <- list()
for (i in 1:ncol(batch1)){
  #Get data for experiment
  id <- which(batch1[,i]>0)
  indicators<- rownames(batch1)[id]
  
  #Make dataset for experiment
  expData <- data%>%
    select(c("Country","SPI_countrycode","SPI_year","Region",
             indicators,
             ".imp"))%>%
    rename(imp=.imp)
  #Run test and save results
  results_batch1[[i]] <- clusterData(expData, m=15,minInformation = 0.8)
}

#saveRDS(results_batch1,"results_batch1.RDS")


#Batch 2 ---------------------------------------
#Get batch of experiments
batch2 <- Experiments%>%
  select(1,9,10)%>%
  as.matrix()%>%
  na.omit() # clear row of unrepresented dimension

#Format names as in data
rownames(batch2) <- colnames(data)[5:(ncol(data)-1)]

#Remove duplicate column
batch2 <- batch2%>%
  data.frame()%>%
  select(-c(1))%>%
  mutate_all(~as.integer(.))

#Define number of experiments to be conducted
results_batch2 <- list()
for (i in 1:ncol(batch2)){
  #Get data for experiment
  id <- which(batch2[,i]>0)
  indicators<- rownames(batch2)[id]
  
  #Make dataset for experiment
  expData <- data%>%
    select(c("Country","SPI_countrycode","SPI_year","Region",
             indicators,
             ".imp"))%>%
    rename(imp=.imp)
  #Run test and save results
  results_batch2[[i]] <- clusterData(expData, m=15,minInformation = 0.95)
}

#saveRDS(results_batch2,"results_batch2.RDS")

#Inspect results --------------------------------------------------
filename <- file.choose()
results <- readRDS(filename)
batch <- batch1

#Calculate
fmis <- data.frame("Exp"=1:length(results),
                   "fmis"=rep(0,length(results)) )
for (i in 1:length(results)){
  #Get data for experiment
  id <- which(batch[,i]>0)
  indicators<- rownames(batch)[id]
  
  #Make dataset for experiment
  raw_indicators <- raw_data%>%
             select(all_of(indicators))
  N <- ncol(raw_indicators)*nrow(raw_indicators) #get total amount of data points
  fmis[i,"fmis"] <- sum(is.na(raw_indicators)) / N * 100 #Get percent that is missing
}

#Extract results
matrices <- lapply(1:3, function(i) {
  do.call(cbind, lapply(results, function(x) x[[i]]))
})

#Visualize how the number of clusters change
#All imputed sets
df1 <- melt(matrices[[2]])
colnames(df1) <- c("Run", "Exp", "Value")

#Consensus clustering
df2 <- melt(lapply(1:length(results), function(i){max(matrices[[1]][,i])}))
colnames(df2) <- c("optNum","Exp")

#Number of PCs
df3 <- melt(matrices[[3]])
colnames(df3) <- c("Run", "Exp", "Value")

ggplot(df1, aes(x = factor(Exp), y = Value)) +
  geom_point(aes(color = "Clusterings"), alpha = 0.1,size=2.5) +
  geom_line(data = df2, aes(x = Exp, y = optNum, color = "Final clustering"), size = 1.5) +
  # Add second dataset, scaled to match first axis
  geom_line(data = df3, aes(x = Exp, y = Value, color = "Number of principal components used"),linetype=2) +
  geom_line(data = fmis, aes(x=Exp, y=fmis, color="% Missing values"),linetype=2)+
  scale_y_continuous(
    name = "Number of clusters",
    sec.axis = sec_axis(~ ., name = "% Missing values")  # inverse transform
  ) +
  scale_color_manual(values = c("Clusterings" = "black", "Final clustering" = "red","% Missing values"="blue" ,"Number of principal components used" = "green")) +
  labs(x = "Last added indicator", title = "First Batch", color = "") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90 )) +
  scale_x_discrete(labels=colnames(batch1))


#Visualize the clusterings---------------------------------------

#Load gis map
library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf",continent = c("south america","oceania","north america", "asia","europe","africa"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'SDZ',"SDN"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'PN1',"PNG"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'PR1',"PRT"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'SSD',"SSD*"))
colnames(world)[57] <- "SPI_countrycode"

#Append clusters to world data
clusterings <- data.frame(SPI_countrycode = unique(data$SPI_countrycode),
           matrices[[1]])
#Name experiments and join
colnames(clusterings)[2:7] <- colnames(batch)
world <- left_join(world, clusterings, 
                   by = "SPI_countrycode")

library(patchwork)
p <- list()
for (i in 1:ncol(batch)){
  p[[i]] <- ggplot() +
    geom_sf(data = world, aes(fill = factor(!!sym(colnames(batch)[i]))), color = "white",size=0.5)+
    labs(title = colnames(batch)[i])+ 
    theme_bw() + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())+
    #guides(fill=guide_legend(title="Cluster"))+
    guides(fill="none")
}
(p[[1]]+p[[2]])/(p[[3]]+p[[4]])/(p[[5]]+p[[6]])


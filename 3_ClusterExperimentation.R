# This script runs clusterings on imputed data using scaled indicators directly from SPI

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

#Read experiment file and set names in the correct order
Experiments <- read_xlsx("IndicatorsForClusters.xlsx",sheet="Experiments")%>%
  na.omit()%>%
  data.frame()
rownames(Experiments) <- colnames(data)[5:(ncol(data)-1)]
rownames(Experiments) <- rownames(Experiments)[c(6,5,4,3,2,1,
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
#Read raw data for reference
raw_data <- read.csv("extendedDataScaled.csv")%>%
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

#Run experiments selected experiments -----------------------------------------
batch_scaled <- Experiments%>%
  select(c(1, "Baseline",
           "BHN_FWB",
           "No_non_DLS",
           "No_nonDLS_limExt",
           "Few.indicators_SPI_preferred",
           "Few.indicators_SPI_preferred_90",
           "Few_indicators_closest_DLS_coverage",
           "Few_indicators_closest_DLS_coverage_90",
           "No_non_DLS_BHNFWB"
           ))%>%
  as.matrix()%>%
  na.omit()

#Remove duplicate column
batch_scaled <- batch_scaled%>%
  data.frame()%>%
  select(-c(1))%>%
  mutate_all(~as.integer(.))

#Define number of experiments to be conducted
results_batch_scaled <- list()
#Define boundary for the pca
minInformation <- c(0.8,0.8,0.8,0.8,0.8,0.9,0.8,0.9,0.8)

for (i in 1:ncol(batch_scaled)){
  #Get data for experiment
  id <- which(batch_scaled[,i]>0)
  indicators<- rownames(batch_scaled)[id]
  
  #Make dataset for experiment
  expData <- data%>%
    select(c("Country","SPI_countrycode","SPI_year","Region",
             indicators,
             ".imp"))%>%
    rename(imp=.imp)
  #Run test and save results
  results_batch_scaled[[i]] <- clusterData(expData, m=15,minInformation = minInformation[i])
}

#saveRDS(results_batch_scaled,"3_ClusterResults.RDS")


#Inspect results --------------------------------------------------
filename <- file.choose()
results <- readRDS(filename)
batch <- batch_scaled #Select which batch to analyse


#Calculate percentage of missing values
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
  geom_point(data = df2, aes(x = Exp, y = optNum, color = "Final clustering"), size = 1.5) +
  # Add second dataset, scaled to match first axis
  geom_line(data = df3, aes(x = Exp, y = Value, color = "Number of principal components used"),linetype=2) +
  geom_line(data = fmis, aes(x=Exp, y=fmis, color="% Missing values"),linetype=2)+
  scale_y_continuous(
    name = "Number of clusters",
    sec.axis = sec_axis(~ ., name = "% Missing values")  # inverse transform
  ) +
  scale_color_manual(values = c("Clusterings" = "black", "Final clustering" = "red","% Missing values"="blue" ,"Number of principal components used" = "green")) +
  labs(x = "", title = "Batch_scaled", color = "") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90 )) +
  scale_x_discrete(labels=colnames(batch))


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
if (batch == batch6){
  clusterings <- data.frame(SPI_countrycode = unique(data$SPI_countrycode)[unique(data$Country)%in% df_colonial$Country],
                            matrices[[1]])
}
#Name experiments and join
colnames(clusterings)[2:(ncol(batch)+1)] <- colnames(batch)
world <- left_join(world, clusterings, 
                   by = "SPI_countrycode")

cluster_names <- 1:8
cluster_colors <- setNames(scales::hue_pal()(9), cluster_names)
library(patchwork)
p <- list()
for (i in 1:ncol(batch)){
  p[[i]] <- ggplot() +
    geom_sf(data = world, aes(fill = factor(!!sym(colnames(batch)[i]))), color = "white",size=0.5)+
    labs(title = gsub("95","90",colnames(batch))[i])+ 
    theme_bw() + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())+
    scale_fill_manual(values = cluster_colors)+
   guides(fill=guide_legend(title="Cluster",ncol=1))
    #guides(fill="none")
}
(p[[1]]+p[[2]]+p[[3]]+p[[4]])/(p[[5]]+p[[6]]+p[[7]]+p[[8]])
p[[1]]+p[[2]]
(p[[1]]+p[[2]])
p[[3]]
p[[1]]

#Join experiments into a dataframe---------------------------------------------------
batch_scaled <- readRDS(file.choose())
batch_scaled <- lapply(1:3, function(i) {
  do.call(cbind, lapply(batch_scaled, function(x) x[[i]]))
})

#Combine into data frame
clusterVariations <- data.frame(batch_scaled[[1]])

colnames(clusterVariations) <- c("Baseline",
                                  "BHN_FWB",
                                  "No_non_DLS",
                                  "No_non_DLS_limExt",
                                  "Few indicators_SPI_preferred",
                                  "Few indicators_SPI_preferred_90",
                                  "Few_indicators_closest_DLS_coverage",
                                  "Few_indicators_closest_DLS_coverage_90",
                                  "No_non_DLS_BHNFWB")

#Add identifier columns
unique(data$Country)==rownames(clusterVariations)
clusterVariations$Country <- unique(data$Country)
clusterVariations$SPI_countrycode <- unique(data$SPI_countrycode)
#Reorder colunms
clusterVariations <- clusterVariations%>%
  relocate(c(Country,SPI_countrycode))

#Save
#saveRDS(clusterVariations, "clusterVariations_laglead_scaled.RDS")


#Grid showing DLS dimension per cluster-----------------------------------------------------

#Comparison of clusters from thesis and from new imputation----------------------------------------------------
filename <- file.choose()
results <- readRDS(filename)
SPI_baseline <- results[[1]]$finalClustering
thesis_clusters <- read.csv("SPIClusters_thesis.csv")
colnames(thesis_clusters) <- c("Country","Original Cluster")

#Join dataset
change <- data.frame("Country"=names(SPI_baseline),"SPI_baseline"=SPI_baseline)
change <- left_join(change,thesis_clusters, by="Country")

#Compare clusterings
cbind(change[which(change$SPI_baseline != change$`Original Cluster`),1],
      change[which(change$SPI_baseline != change$`Original Cluster`),2],
      change[which(change$SPI_baseline != change$`Original Cluster`),3])
      
#Visualize 
scaled_data <- data

#Visualize how this has affected clusters
dimmeans <- scaled_data%>%
  mutate(BHN_mean = rowSums(across(5:22))/18,
         FWB_mean = rowSums(across(23:37))/15,
         OPP_mean = rowSums(across(38:56))/19)
#Join clusters 
dimmeans <- left_join(dimmeans,change, by ="Country")

dimmeans_long <- dimmeans %>%
  pivot_longer(cols = c(BHN_mean, FWB_mean,OPP_mean), names_to = "variable", values_to = "value")

ggplot(dimmeans_long, aes(x = SPI_year, y = value, group=Country,color=factor(SPI_baseline))) +
  geom_line(alpha=1)+
facet_grid(SPI_baseline~variable,labeller = label_value) +
  #theme_minimal() +
  #scale_color_manual(values = cluster_colors)+
  labs(title = "", x = "", y = "")+
  guides(color="none")

ggplot(dimmeans_long, aes(x = SPI_year, y = value, group=Country,color=factor(SPI_baseline))) +
  geom_line(alpha=1)+
  facet_grid(1~variable,labeller = label_value)




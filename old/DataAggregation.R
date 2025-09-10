#Load packages
library(dplyr)
library(mice)

#For clustering test
library(NbClust)
library(partitionComparison)
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

#Data aggregation ========================================

bhn_original <-read.csv("Data/bhn.csv")

#Load imputed data
data <- read.csv("ImputedData.csv")%>%
  select(-X)

#Get raw data
data_raw <- data%>%
  filter(imp==0)%>%
  select(-c("POP_71","GDPPPP_72","FOSSILTES_75","TESJ_77","TFCJ_78",
            "logPOP_71","logGDPPPP_72","logTESJ_77","logTFCJ_78" ))


#Load new sources
LifeExp <- read.csv("Data/LifeExpectancyAtBirth.csv")%>%
  select(Location,Time,Value)
HealthExpend <- read.csv("Data/CurrentHealthExpenditure.csv", skip = 3) %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Year",
    values_to = "Expenditure",
    names_transform = list(Year = ~sub("^X", "", .x))  # remove the X
  )

SuicideRate <- read.csv("Data/SuicideData.csv")

PrimarySchoolTeachers <-read.csv("Data/PrimarySchoolTeachers.csv")

#Load LMIC data
load("Data/Complete DLS data file_country level.RData")
dls_country <- dls_country %>%
  rename(Country_0=country_name,
         SPI_year_2 = year_start_interviews)%>%
  select(-c(wave,dls_index_country,country_code,dls_min10_country,dls_min7_country,dls_min5_country))%>%
  mutate(SPI_year_2 = as.integer(SPI_year_2))
  


#Join new sources with old raw data
new_data <- data_raw%>%
  left_join(dls_country,by=c("Country_0", "SPI_year_2"))


#impute
ini <- mice(new_data, maxit=0)
imp <- mice(new_data,
            m=5,
            maxit=10)
densityplot(imp, as.formula(paste("~",paste0(colnames(new_data)[69:78],collapse="+"))) ,layout = c(3, 4))

#Get imputed data 
imp_data <- complete(imp,"long")%>%
  select(-c(imp,.id))%>%
  rename(imp=.imp)
#Save it
write.csv(imp_data,"imp_data_test.csv")

#Try to cluster 
testCluster <- clusterData(imp_data[,c(1:4,68:78)],m=5)
testCluster$finalClustering




#Plot
library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf")%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'SDZ',"SDN"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'PN1',"PNG"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'PR1',"PRT"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'SSD',"SSD*"))
colnames(world)[57] <- "SPI_countrycode_1"

#Append clusters to world data
world <- left_join(world, data.frame(SPI_countrycode_1 = unique(imp_data$SPI_countrycode_1), ClusterPC = testCluster$finalClustering), by = "SPI_countrycode_1")
world$ClusterPC <- factor(world$ClusterPC)

#Plot
library(ggplot2)
p1 <- ggplot() +
  geom_sf(data = world, aes(fill = ClusterPC), color = "white")+
  labs(title = "Clustering of DLS dimensions")+ 
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  guides(fill=guide_legend(title="Cluster"))
p1

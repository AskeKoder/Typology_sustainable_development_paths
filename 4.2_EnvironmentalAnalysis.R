library(tidyr)
library(dplyr)
library(readxl)
library(countrycode)


#load data
data <- read.csv("ImputedDataLag1Lead2_maxit30_scaled.csv")%>%
  select(-c("X",".id"))%>%
  relocate(.imp, .after=last_col())

#Clusters from experiments
clusterSelection <- "Baseline"
clusters <- readRDS("clusterVariations_laglead_scaled.RDS")%>%
  select(Country,iso3 = SPI_countrycode,cluster = clusterSelection)

#Read experiment file
Experiments <- read_xlsx("IndicatorsForClusters.xlsx",sheet="Experiments")

#link experiments with included indicators
batch_scaled <- Experiments%>%
  select(c("Baseline",
           "BHN_FWB",
           "No_non_DLS",
           "No_nonDLS_limExt",
           "Few_indicators_SPI_preferred",
           "Few_indicators_SPI_preferred_90",
           "Few_indicators_closest_DLS_coverage",
           "Few_indicators_closest_DLS_coverage_90",
           "No_non_DLS_BHNFWB"
  ))%>%
  as.matrix()%>%
  na.omit()

rownames(batch_scaled) <- colnames(data)[5:(ncol(data)-1)]
rownames(batch_scaled) <- rownames(batch_scaled)[c(6,5,4,3,2,1,
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

batch_scaled[,c(1:ncol(batch_scaled))] <- as.numeric(batch_scaled[,c(1:ncol(batch_scaled))])

#Merge data sets
clusteredData <- data %>%
  rename(iso3 = SPI_countrycode)%>%
  left_join(clusters , by= c("Country","iso3"))%>%
  data.frame()%>%
  mutate(cluster = as.factor(cluster))



#Get GDP data
raw_GDP <- read.csv("Data_WellBeing/GDPpercap PPP 2001 international world bank.csv", header = FALSE, stringsAsFactors = FALSE)
GDP <- raw_GDP[-c(1:3), ] #Remove metadata
colnames(GDP) <- raw_GDP[3, ] #Set colnames

#Prepare
GDP <- GDP %>%
  pivot_longer(
    cols = matches("^\\d{4}$"),  # four numbers(\\d{4}) between ^start and $end of string 
    names_to = "Year",
    values_to = "GDP_PPP_current_international_dollars"
  ) %>%
  dplyr::select(
    Country = `Country Name`,
    ISO_Country = `Country Code`,
    Year,
    GDP_PPP_current_international_dollars)%>%
  filter(Year %in% 2000:2020)%>%
  mutate(Country = as.factor(Country),
         ISO_Country = as.factor(ISO_Country),
         Year = as.numeric(Year))%>%
  rename(GDP_PPPcap = GDP_PPP_current_international_dollars)

#Standardize names 
GDP$Country_std  <- countrycode(GDP$Country, origin="country.name",destination="country.name")
GDP$iso3  <- countrycode(GDP$Country, origin="country.name",destination="iso3c")
#Print where the renaming failed
print(unique(GDP[is.na(GDP$Country_std), "Country"]),n=100)  #Only aggregated countries failed. No problem




#Get HDI data
HDI <- read.csv("Data_WellBeing/human-development-index.csv", header = TRUE, stringsAsFactors = FALSE)%>%
  rename(Country=Entity,
         ISO_Country = Code,
         HDI = Human.Development.Index)%>%
  mutate(Country = as.factor(Country),
         ISO_Country = as.factor(ISO_Country))%>%
  filter(Year %in% 2000:2020)

summary(HDI)
str(HDI)
HDI$Country_std  <- countrycode(HDI$Country, origin="country.name",destination="country.name")
HDI$iso3  <- countrycode(HDI$Country, origin="country.name",destination="iso3c")
print(unique(HDI[is.na(HDI$Country_std), "Country"])) #Only aggrregates and micronesia,  which is not included in the clustering




#Environmental data
data <- read_xlsx("Data_Environmental/GCSI_59a_Per_capita_2000-2020_05162025.xlsx") %>%
  pivot_wider(id_cols = c(Region_acronyms, Region_names, Year, Population),
              names_from = Indicator, values_from = "Value_per_Capita")%>%
  rename(Country=Region_names)%>%
  select(-Region_acronyms)


data[data$Country=="Yugoslavia/Serbia (1991/1992)",2] <- "Serbia"
data[data$Country=="Yemen Arab Republic/Yemen (1990/1991)",2] <- "Yemen"

CombinedData <- data %>%
  #Standardize countries
  mutate(Country= countrycode(Country, origin="country.name",destination="country.name"))%>%
  mutate(iso3 = countrycode(Country, origin="country.name",destination="iso3c"))%>%
  
  
  mutate(Year = as.integer(Year))%>%
  rename(GHG=ghg_ar6,
         Biodiversity_Impact = land_threats,
         Scarce_Water_Consumption = water_tot,
         Water_Stress = water_stress_tot,
         Black_Carbon = bc,
         NOx = nox,
         SOx = sox,
         Pm25 = pm25,
         Pm10 = pm10
  )%>%
  mutate(GHG = 1000*GHG)%>% #From Gigagrams/cap to t/cap
  mutate(Scarce_Water_Consumption = 1000000*Scarce_Water_Consumption)%>% #From Mm3 to m3
  mutate(NOx = 1000000*NOx)%>% #From Gigagrams/cap to kg/cap
  
  #Merging data sets
  left_join(clusteredData%>%rename(Year = SPI_year),by=c("Country","iso3","Year"))%>%
  left_join(HDI, by=c("Country","iso3","Year"))%>%
  left_join(GDP,by=c("Country","iso3","Year"))%>%
  dplyr::select(!ends_with(c(".x",".y")))%>%
  #dplyr::select(!ISO_Country)%>%
  # rename(Country = Country_std,
  #        ISO_Country = ISO_Country_std)%>%
  
  #Edit variables
  mutate(GDP_PPP = GDP_PPPcap*Population)%>%
  mutate(Year=as.integer(Year))%>%
  # mutate(Country=as.factor(Country),
  #        ISO_Country=as.factor(ISO_Country))%>%
  
  #Filter data
  filter(Year %in% 2000:2020)%>%
  filter(!is.na(Country))%>%
  filter(!is.na(cluster))%>%
  # dplyr::select(c(Country,ISO_Country,Year,Population,Cluster,ExtCluster,
  #                 GHG,Scarce_Water_Consumption, Biodiversity_Impact,
  #                 NOx,GDP_PPPcap,HDI))%>% #Select indicators
  arrange(Country,Year)%>%
  data.frame()


length(unique(CombinedData$Country))
length(unique(CombinedData$cluster))

#Set up DLS performance ----------------------------------------------
dimIndicators <- data.frame("Housing"="Share_not_in_Slums",
                            "Thermal comfort" = "Access_to_electricity",
                            "Food and Nutrition" = "Nutritional_deficiencies",
                            "Food prep/storage" = "Prevalence_of_cooking_with_coalperbiomass",
                            "Water" = "Safely_Managed_Drinking_Water",
                            "Sanitation" = "Safely_managed_saniation",
                            "Provision of health care" = "Universal_health_coverage",
                            "Education" = "Prim_School_Enroll",
                            "Social connectedness" = "Mobile_and_landline_telephone_subscriptions")


DLSdata <- CombinedData %>%
  mutate(Share_not_in_Slums = 100 - Share_Slums)%>%
  #select(Country,iso3,SPI_year,cluster,.imp, all_of(c(as.matrix(dimIndicators))))%>%
  mutate(min_val = do.call(pmin, c(across(all_of(c(as.matrix(dimIndicators)))), na.rm = TRUE)))

DLSdata$which_min <- ""
for (i in 1:nrow(DLSdata)){
  DLSdata$which_min[i] <- c(as.matrix(dimIndicators))[which(DLSdata[i,c(as.matrix(dimIndicators))]==DLSdata$min_val[i])]
}

#IMPORTANT Getting just one imputation, this needs to be fixed if used 
DLSdata_summarized<- DLSdata%>%
  group_by(Year,Country,cluster)%>%
  mutate(min_est = mean(min_val))%>%
  filter(.imp==1) %>%
  ungroup()%>%
  mutate(Country=factor(Country))


#Fit Mixed MODEL ---------------------------------------------
library(lme4)
fit <- lmer (min_est~ scale(log(GDP_PPPcap),scale=FALSE)*cluster+(scale(log(GDP_PPPcap),scale=FALSE)|Country),
          data=DLSdata_summarized)
summary(fit)

newdata <- expand.grid(
  GDP_PPPcap = seq(min(DLSdata_summarized$GDP_PPPcap, na.rm=TRUE), 
          max(DLSdata_summarized$GDP_PPPcap, na.rm=TRUE), 
          length.out=50),
  cluster = unique(DLSdata_summarized$cluster),
  which_min = unique(DLSdata_summarized$which_min)[1]  # pick one or loop over
)

for (i in 1:11){
  newdata <- newdata%>%
    filter(!(cluster == i & GDP_PPPcap<min(DLSdata_summarized$GDP_PPPcap[DLSdata_summarized$cluster==i], na.rm=TRUE)))%>%
    filter(!(cluster == i & GDP_PPPcap>max(DLSdata_summarized$GDP_PPPcap[DLSdata_summarized$cluster==i], na.rm=TRUE)))
}



preds <- cbind(newdata,"fit"=predict(fit,newdata,re.form = NA))

# preds <- ggeffects::ggpredict(fit, terms = c("GDP_PPPcap", "cluster"))%>%
#   data.frame()
library(ggplot2)
ggplot()+
  geom_line(data=preds, aes(x=log(GDP_PPPcap),y=fit,group=cluster,color=cluster))+
  #geom_ribbon(data=preds, aes(x=log(GDP_PPPcap),ymin=lwr, ymax=upr,fill=cluster),alpha=0.3)+
  geom_point(data=DLSdata_summarized,aes(x=log(GDP_PPPcap),y=min_est,color=cluster),alpha=0.1)
  #coord_cartesian(xlim = range(log(DLSdata_summarized$GDP_PPPcap), na.rm = TRUE))+
  #facet_wrap(~cluster, scales="free_x")






#Fit S curve-------------------------------
#Error function
MSE <- function(obs,preds){
  sum((preds-obs)^2)/length(preds)
}

#S curve
sCurve <- function(mf,k,W,C){
  sCurve <- C/(1+exp(-k*(mf-W)))
}

#sample grid
#k= slope a inflection point
#W= logarithmic value of resource use at inflection point
sampleGrid = expand.grid(
  k = seq(0,100,length.out=400),
  W = seq(min(log(DLSdata_summarized$GHG),na.rm=TRUE),max(log(DLSdata_summarized$GHG),na.rm=TRUE),length.out=300)
)

#Function for fitting model
fitScurve <- function(w, mf, sampleGrid=sampleGrid){
  #List for saving results
  Err  <- rep(NA,nrow(sampleGrid))
  
  #Estimate
  k <- sampleGrid$k
  W <- sampleGrid$W
  for (i in 1:nrow(sampleGrid)){
    est <- sCurve(mf,k=k[i],W=W[i],C=max(w,na.rm=TRUE))
    
    #Calculate Errors
    Err[i] <- MSE(w,est)
  }
  
  #Return parameters with lowest MSE
  sampleGrid[which.min(Err),]
}

#Prediction with S curve
predScurve <- function(w,mf,k,W){
  C <- max(w,na.mr=TRUE)
  w_est <- C/(1+exp(-k*(mf-W)))
}


#Fit model for nutritional deficiencies
params <- fitScurve(DLSdata_summarized$Nutritional_deficiencies,DLSdata_summarized$GHG,sampleGrid = sampleGrid)
#Predict and plot
pred <- predScurve(DLSdata_summarized$Nutritional_deficiencies,DLSdata_summarized$GHG,k=params$k,params$W)
pred

plottingData <- cbind(DLSdata_summarized,pred)
ggplot()+
  geom_path(data=plottingData%>%
              filter(Year %in% c(2000,2020)), aes(x=log(GHG),y=Nutritional_deficiencies,color=cluster,group=Country))+
  geom_point(data=plottingData%>%
              filter(Year %in% c(2000,2020)), aes(x=log(GHG),y=Nutritional_deficiencies,color=cluster,shape=factor(Year)))+
  geom_line(data=plottingData, aes(x=log(GHG),y=pred))




#Calculate Scurve for all DLS dimensions
all_results <- list()
for (xvar in c(as.matrix(dimIndicators))) {
  tmp <- DLSdata_summarized %>%
    select(GHG, x = !!sym(xvar),Year,Country,cluster) %>%
    mutate(variable = xvar)
  
  params <- fitScurve(tmp$x, tmp$GHG, sampleGrid = sampleGrid)
  pred   <- predScurve(tmp$x, tmp$GHG, k = params$k, params$W)
  
  all_results[[xvar]] <- cbind(tmp,pred)
}


#Join results and plot
results <- bind_rows(all_results)%>%
  mutate(cluster=(cluster))%>%
  #filter(cluster==9)
  filter(cluster%in% 1:11)

cluster_names <- 1:11
cluster_colors <- setNames(scales::hue_pal()(length(unique(clusters$cluster))), cluster_names)

ggplot()+
  geom_path(data=results%>%
              filter(Year %in% c(2000,2020)), aes(x=GHG,y=x,color=cluster,group=Country),
            alpha=0.4)+
  geom_point(data=results%>%
               filter(Year %in% c(2000,2020)), aes(x=GHG,y=x,color=cluster,shape=factor(Year)))+
  geom_line(data=results, aes(x=GHG,y=pred))+
  facet_wrap(~variable)+
  #geom_vline(xintercept=3, lty=2)+ #at 3 ton co2/cap
  scale_color_manual(values = cluster_colors)+
  coord_trans(x = "log10")+
  labs(x="t CO_2 cap", y="%")
  #theme_minimal()


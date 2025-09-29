#This script analyses the clusters resulting from using lag lead imputation on scaled indicators
library(dplyr)
library(countrycode) #For standardizing country names

data <- read.csv("ImputedDataLag1Lead2_maxit30_scaled.csv")%>%
  select(-c("X",".id"))%>%
  relocate(.imp, .after=last_col())

#Clusters from experiments
clusterSelection <- "No_non_DLS_BHNFWB"
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
           "Few indicators_SPI_preferred",
           "Few indicators_SPI_preferred_90",
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



#Compare clusters with HDI and GDP--------------------------------------------------
#GDP (PPP) data
raw_GDP <- read.csv("Data_WellBeing/GDPpercap PPP 2001 international world bank.csv", header = FALSE, stringsAsFactors = FALSE)
GDP <- raw_GDP[-c(1:3), ] #Remove metadata
colnames(GDP) <- raw_GDP[3, ] #Set colnames


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
  rename(GDP_PPP = GDP_PPP_current_international_dollars)

summary(GDP)
str(GDP)

#Standardize names 
GDP$Country_std  <- countrycode(GDP$Country, origin="country.name",destination="country.name")
GDP$iso3  <- countrycode(GDP$Country, origin="country.name",destination="iso3c")
#Print where the renaming failed
print(unique(GDP[is.na(GDP$Country_std), "Country"]),n=100)  #Only aggregated countries failed. No problem



#HDI data
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


#Join clusters and indices
HDI_clust <- left_join(clusters,HDI, by=c("iso3","Country"))%>%
  left_join(GDP,by=c("iso3","Country","Year"))%>%
  select(Country,iso3,cluster,Year,HDI,GDP = GDP_PPP)%>%
  mutate(cluster=factor(cluster))%>%
  na.omit()


hulls <- HDI_clust %>%
  group_by(cluster) %>%
  slice(chull(log(GDP), HDI)) %>%
  ungroup()

ggplot(HDI_clust,aes(x=log(GDP) ,y=HDI, color=factor(cluster), group=factor(Country)))+
  geom_point(alpha=0.5)+
  # geom_path(alpha=0.3)+
  geom_polygon(data = hulls, aes(log(GDP), HDI, fill = cluster, group = cluster),
               alpha = 0.2, color = NA)+
  theme_minimal()



#Plot median of indicators over time -----------------------------------------------
selected_vars <- rownames(batch_scaled)[batch_scaled[, clusterSelection] == 1]
set <- clusteredData %>%
  group_by(SPI_year, cluster) %>%
  dplyr::select(c("SPI_year", "cluster",
                  selected_vars))%>%
  summarise(across(all_of(selected_vars),
                   list(
                     median = ~median(., na.rm = TRUE),
                     q1 = ~quantile(., 0.25, na.rm = TRUE),
                     q3 = ~quantile(., 0.75, na.rm = TRUE)
                   ),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")
  

#Get names of variables
variables <- selected_vars


#Convert to long format for plotting
long_median <- set %>%
  select(1,2, ends_with("_median"))%>%
  pivot_longer(cols = ends_with("_median"),
               names_to = "Variable",
               values_to = "Median",
               names_pattern = "(.*)_median")
long_q1 <- set %>%
  select(1,2, ends_with("_q1"))%>%
  pivot_longer(cols = ends_with("_q1"),
               names_to = "Variable",
               values_to = "Q1",
               names_pattern = "(.*)_q1")
long_q3 <- set %>%
  select(1,2, ends_with("_q3"))%>%
  pivot_longer(cols = ends_with("_q3"),
               names_to = "Variable",
               values_to = "Q3",
               names_pattern = "(.*)_q3")
# Combine all
long_set <- long_median %>%
  left_join(long_q1, by = c("SPI_year", "cluster", "Variable")) %>%
  left_join(long_q3, by = c("SPI_year", "cluster", "Variable")) %>%
  mutate(Variable = factor(Variable, levels = variables))


ggplot(long_set, aes(x = SPI_year, y = Median, color = cluster, fill = cluster)) +
  geom_ribbon(aes(ymin = Q1, ymax = Q3), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Variable, scales = "free_y") +
  theme_minimal() +
  # scale_fill_manual(values = cluster_colors)+
  # scale_color_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Year", y = "Cluster Median (with IQR)")+
  theme(plot.title = element_text(size = 5))


long_median_ranked <- long_median %>%
  group_by(SPI_year, Variable) %>%  # Group by year and variable (NOT group)
  mutate(Rank = case_when(
    Variable == "Share_Slums" ~ rank(-Median, ties.method = "min"),   # lower = better
    TRUE                       ~ rank(Median, ties.method = "min")   # higher = better
  ))%>%
  ungroup()

ggplot(long_median_ranked, aes(x = SPI_year, y = Rank, color = cluster, fill = cluster)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Variable, scales = "free_y") +
  ylim(0, 11) +
  theme_minimal() +
  # scale_fill_manual(values = cluster_colors)+
  # scale_color_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Year", y = "Cluster Median (with IQR)")+
  theme(plot.title = element_text(size = 5))

ggplot()+
  geom_density(data=long_median_ranked,aes(x=Rank, group=cluster, color=cluster))

factor(clusters$cluster,levels=c(10,4,1,8,9,7,3,2,5,11,6))



#Plot radarcharts -----------------------------------------------------------
library(fmsb)
plot_list <- list()
for (i in 1:8){
  cluster <- i
  clustering = "No_non_DLS_BHNFWB"
  clusterVariations <- readRDS("clusterVariations_laglead_scaled.RDS")
  
  country_series <- data %>%
    filter(SPI_countrycode %in% clusterVariations[clusterVariations[,clustering]==cluster,"SPI_countrycode"])
  
  #Calculate mean across imputed sets
  country_series_filtered <- country_series%>%
    filter(.imp!=0)%>%
    dplyr::select(!c(Country,SPI_countrycode,Region))%>%
    group_by(SPI_year) %>%
    summarise_all(mean)
  
  #Select indicators
  # id <- which(batch_scaled[,1]>0)
  # indicators<- rownames(batch_scaled)[id]
  
  country_series_filtered <- country_series_filtered%>%
    select(SPI_year,all_of(indicators))
  
  #Add min max levels to data for plotting
  mins <- data %>%
    summarize_all(min,na.rm=TRUE)%>%
    select(SPI_year,all_of(indicators))
  
  maxs <- data %>%
    summarize_all(max,na.rm=TRUE)%>%
    select(SPI_year,all_of(indicators))
  
  plottingData <- rbind(
    maxs[2:35],
    mins[2:35],
    country_series_filtered[2:35]
  )
  
  #Sort columns according to DLS dimension
  order <- order(Experiments[rownames(Experiments)%in% colnames(plottingData),"Related.DLS.dimension"])
  plottingData <- plottingData[,order]%>%
    mutate(Share_Slums = 100-Share_Slums)%>%
    rename(Share_not_in_slums=Share_Slums)
  
  plottingData[1, "Share_not_in_slums"] <- 100
  plottingData[2, "Share_not_in_slums"] <- 0
  #Define color scale
  library(scales)
  cluster_names <- 1:8
  cluster_colors <- setNames(scales::hue_pal()(9), cluster_names)
  
  colors_border <- col_numeric(palette = c("grey",cluster_colors[cluster],cluster_colors[cluster],cluster_colors[cluster],"black"), domain = NULL)(0:20)
  colors_fill <- adjustcolor(colors_border, alpha.f = 0.4)
  
  
  
  par(mai=c(0,0,0.5,0))
  png(filename = paste("Figures/RadarCharts/NoNonDLSBHNFWB/",i,".png"), width=700, height=700)
  radarchart(plottingData, vlabels = colnames(plottingData),
             na.itp = FALSE, col=rep(3,52),
             pcol = colors_border,
             axistype = 0,
             caxislabels = "",
             plwd = 2,
             plty = 1,
             cglty=1,
             cglcol = "grey60",
             pty="",
             title=paste0("Cluster ",cluster),
             cex.main=2)
  dev.off()
  
}
#This script analyses the clusters resulting from using lag lead imputation on scaled indicators
library(countrycode) #For standardizing country names
library(rnaturalearth)
library(ggplot2)
library(tidyr)
library(readxl)
library(dplyr)


data <- read.csv("ImputedDataLag1Lead2_maxit30_scaled.csv")%>%
  select(-c("X",".id"))%>%
  relocate(.imp, .after=last_col())

#Clusters from experiments
clusterSelection <- 'Few indicators_SPI_preferred'
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
  mutate(cluster = as.factor(cluster))%>%
  ungroup()

table(clusters$cluster)

#PLot worldmap--------------------------------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf",continent = c("south america","oceania","north america", "asia","europe","africa"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'SDZ',"SDN"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'PN1',"PNG"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'PR1',"PRT"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'SSD',"SSD*"))
colnames(world)[57] <- "iso3"

#Name experiments and join
world <- left_join(world, clusters, 
                   by = "iso3")

cluster_names <- 1:length(unique(clusters$cluster))
cluster_colors <- setNames(scales::hue_pal()(length(unique(clusters$cluster))), cluster_names)
ggplot() +
    geom_sf(data = world, aes(fill = factor(cluster)), color = "white",size=0.5)+
    labs(title = "SPI Baseline")+ 
    theme_bw() + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())+
    scale_fill_manual(values = cluster_colors)+
    guides(fill=guide_legend(title="Cluster",ncol=2))
  #guides(fill="none")



#Compare clusters with HDI table-------------------------------------------------------------------
#Load HDI data 2022
HDI_raw <- read_xlsx("Data_WellBeing/HDR23-24_Statistical_Annex_HDI_Table (1).xlsx")%>%
  setNames(1:length(.))%>%
  select(!c(4,6,8,10,12,14,15))%>%
  setNames(.[4,])
colnames(HDI_raw)[1] <- "HDI rank"
colnames(HDI_raw)[2] <- "Country"

HDI <- HDI_raw %>%
  slice(-c(1,2,3,4,5,6))%>%
  drop_na()

str(HDI)
HDI$`Human Development Index (HDI)` <- as.numeric(HDI$`Human Development Index (HDI)`)
HDI$`Life expectancy at birth` <- as.numeric(HDI$`Life expectancy at birth`)
HDI$`Expected years of schooling` <- as.numeric(HDI$`Expected years of schooling`)
HDI$`Mean years of schooling` <- as.numeric(HDI$`Mean years of schooling`)
HDI$`Gross national income (GNI) per capita`<- as.numeric(HDI$`Gross national income (GNI) per capita`)
HDI$`GNI per capita rank minus HDI rank`<- as.numeric(HDI$`GNI per capita rank minus HDI rank`)

HDI
str(HDI)

# Original country names (replace with your full list)
original_countries <- c(
  "Switzerland", "Norway", "Iceland", "Hong Kong, China (SAR)", "Denmark", "Sweden",
  "Germany", "Ireland", "Singapore", "Australia", "Netherlands", "Belgium",
  "Finland", "Liechtenstein", "United Kingdom", "New Zealand", "United Arab Emirates",
  "Canada", "Korea (Republic of)", "Luxembourg", "United States", "Austria",
  "Slovenia", "Japan", "Israel", "Malta", "Spain", "France", "Cyprus", "Italy",
  "Estonia", "Czechia", "Greece", "Bahrain", "Andorra", "Poland", "Latvia",
  "Lithuania", "Croatia", "Qatar", "Saudi Arabia", "Portugal", "San Marino",
  "Chile", "Slovakia", "Türkiye", "Hungary", "Argentina", "Kuwait", "Montenegro",
  "Saint Kitts and Nevis", "Uruguay", "Romania", "Antigua and Barbuda",
  "Brunei Darussalam", "Russian Federation", "Bahamas", "Panama", "Oman", "Georgia",
  "Trinidad and Tobago", "Barbados", "Malaysia", "Costa Rica", "Serbia", "Thailand",
  "Kazakhstan", "Seychelles", "Belarus", "Bulgaria", "Palau", "Mauritius", "Grenada",
  "Albania", "China", "Armenia", "Mexico", "Iran (Islamic Republic of)", "Sri Lanka",
  "Bosnia and Herzegovina", "Saint Vincent and the Grenadines", "Dominican Republic",
  "Ecuador", "North Macedonia", "Cuba", "Moldova (Republic of)", "Maldives", "Peru",
  "Azerbaijan", "Brazil", "Colombia", "Libya", "Algeria", "Turkmenistan", "Guyana",
  "Mongolia", "Dominica", "Tonga", "Jordan", "Ukraine", "Tunisia", "Marshall Islands",
  "Paraguay", "Fiji", "Egypt", "Uzbekistan", "Viet Nam", "Saint Lucia", "Lebanon",
  "South Africa", "Palestine, State of", "Indonesia", "Philippines", "Botswana",
  "Jamaica", "Samoa", "Kyrgyzstan", "Belize", "Venezuela (Bolivarian Republic of)",
  "Bolivia (Plurinational State of)", "Morocco", "Nauru", "Gabon", "Suriname", "Bhutan",
  "Tajikistan", "El Salvador", "Iraq", "Bangladesh", "Nicaragua", "Cabo Verde", "Tuvalu",
  "Equatorial Guinea", "India", "Micronesia (Federated States of)", "Guatemala",
  "Kiribati", "Honduras", "Lao People's Democratic Republic", "Vanuatu",
  "Sao Tome and Principe", "Eswatini (Kingdom of)", "Namibia", "Myanmar", "Ghana",
  "Kenya", "Nepal", "Cambodia", "Congo", "Angola", "Cameroon", "Comoros", "Zambia",
  "Papua New Guinea", "Timor-Leste", "Solomon Islands", "Syrian Arab Republic", "Haiti",
  "Uganda", "Zimbabwe", "Nigeria", "Rwanda", "Togo", "Mauritania", "Pakistan",
  "Côte d'Ivoire", "Tanzania (United Republic of)", "Lesotho", "Senegal", "Sudan",
  "Djibouti", "Malawi", "Benin", "Gambia", "Eritrea", "Ethiopia", "Liberia",
  "Madagascar", "Guinea-Bissau", "Congo (Democratic Republic of the)", "Guinea",
  "Afghanistan", "Mozambique", "Sierra Leone", "Burkina Faso", "Yemen", "Burundi",
  "Mali", "Chad", "Niger", "Central African Republic", "South Sudan", "Somalia"
)

# Named vector for replacements
replacements <- c(
  "Hong Kong, China (SAR)" = "Hong Kong",
  "Iran (Islamic Republic of)" = "Iran",
  "Korea (Republic of)" = "Korea, Republic of",
  "Moldova (Republic of)" = "Moldova",
  "Viet Nam" = "Vietnam",
  "Palestine, State of" = "West Bank and Gaza",
  "Eswatini (Kingdom of)" = "Eswatini",
  "Micronesia (Federated States of)" = "Micronesia",
  "Lao People's Democratic Republic" = "Laos",
  "Côte d'Ivoire" = "Côte d'Ivoire",  # already correct
  "Tanzania (United Republic of)" = "Tanzania",
  "Russian Federation" = "Russia",
  "Türkiye" = "Turkey",
  "Syrian Arab Republic" = "Syria",
  "Venezuela (Bolivarian Republic of)" = "Venezuela",
  "Bolivia (Plurinational State of)" = "Bolivia",
  "Gambia" = "Gambia, The",
  "Congo (Democratic Republic of the)" = "Congo, Democratic Republic of",
  "Congo" = "Congo, Republic of",
  "North Macedonia"="Republic of North Macedonia"
)

# Apply replacements
standardized_countries <- ifelse(HDI$Country %in% names(replacements),
                                 replacements[HDI$Country],
                                 HDI$Country)

# Show the result
standardized_countries
HDI$Country <- standardized_countries

#Join HDI data with clusters
HDI_clust <- left_join(HDI, data.frame(Country = clusters$Country,
                                       cluster = clusters$cluster,
                                       by = "Country"))

#Taiwan and North korea are not included in the HDI data
#Calculate means and confidence intervals
CI <- function(X){
  mn <- mean(X)
  se <- sd (X)/sqrt(length(X))
  res <- tibble(
    mean = mn,
    add =  qt(1 - (0.05 / 2), length(X) - 1) * se
  )
  res <- round(res,2)
  paste0(res[1]," ±",res[2])
}

summary <- HDI_clust%>%
  drop_na(cluster)%>%
  group_by(cluster)%>%
  summarise_at(3:7,CI)

for (i in 1:nrow(summary)){
  print(paste(summary[i,],collapse=" & "))
}



#Compare clusters with HDI and GDP over time--------------------------------------------------
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
  geom_path(alpha=0.3)+
  geom_polygon(data = hulls, aes(log(GDP), HDI, group = cluster, fill = cluster),
               alpha = 0.2, color = NA)+
  theme_minimal()+
  guides(fill="none")
  



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

#Filter
long_set_filtered <- long_set%>%
  filter(cluster%in%c(4,8))
  #filter(cluster%in%c(1,7,11))
  #filter(cluster%in%c(9,10))
  #filter(cluster%in%c(2,3,5,6))


cluster_names <- 1:length(unique(clusters$cluster))
cluster_colors <- setNames(scales::hue_pal()(length(unique(clusters$cluster))), cluster_names)
ggplot(long_set_filtered, aes(x = SPI_year, y = Median, color = cluster, fill = cluster)) +
  geom_ribbon(aes(ymin = Q1, ymax = Q3), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Variable) +
  theme_minimal() +
  scale_fill_manual(values = cluster_colors)+
  scale_color_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Year", y = "Cluster Median (with IQR)")+
  theme(plot.title = element_text(size = 5))

#With ranks:
final_ranks <- long_set_filtered %>%
  group_by(Variable) %>%
  filter(SPI_year == max(SPI_year)) %>%          # keep only the last year per Variable
  arrange(desc(Median)) %>%              # highest median = rank 1
  mutate(rank = row_number()) %>%
  ungroup()

label_df <- final_ranks %>%
  mutate(
    label = paste0(cluster),           # e.g., "1. ClusterA"
    x = max(long_set_filtered$SPI_year) + 0.5,         # slightly beyond last year
    y = Median                                    # position based on final median
  )

ggplot(long_set_filtered, aes(x = SPI_year, y = Median, color = cluster, fill = cluster)) +
  geom_ribbon(aes(ymin = Q1, ymax = Q3), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Variable, scales = "free_y") +
  theme_minimal() +
  scale_fill_manual(values = cluster_colors) +
  scale_color_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Year", y = "Cluster Median (with IQR)") +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label, color = cluster),
    hjust = 0, size = 3, show.legend = FALSE
  ) +
  expand_limits(x = max(long_set_filtered$SPI_year) + 1) +  # add space for labels
  theme(plot.title = element_text(size = 5))

#Plot radarcharts -----------------------------------------------------------
library(fmsb)
plot_list <- list()
clustering = "Baseline"
nMax <- length(unique(clusters$cluster))
for (i in 1:nMax){
  cluster <- i
  #Filter out other clusters
  country_series <- clusteredData %>%
    filter(cluster == as.character(i))
  
  #Calculate mean across imputed sets
  country_series_filtered <- country_series%>%
    filter(.imp!=0)%>%
    dplyr::select(!c(Country,iso3,Region,cluster))%>%
    group_by(SPI_year) %>%
    summarise_all(mean)
  
  #Select indicators
  id <- which(batch_scaled[,clustering]>0) #Get indicators from batch
  indicators<- rownames(batch_scaled)[id]
  country_series_filtered <- country_series_filtered%>%
    select(all_of(indicators))
  
  #Add min max levels to data for plotting
  mins <- data %>%
    summarize_all(min,na.rm=TRUE)%>%
    select(all_of(indicators))
  
  maxs <- data %>%
    summarize_all(max,na.rm=TRUE)%>%
    select(all_of(indicators))
  
  plottingData <- rbind(
    maxs[id],
    mins[id],
    country_series_filtered[id]
  )
  
  #Sort columns according to DLS dimension
  # order <- order(Experiments[rownames(Experiments)%in% colnames(plottingData),"Related.DLS.dimension"])
  # plottingData <- plottingData[,order]%>%
  #   mutate(Share_Slums = 100-Share_Slums)%>%
  #   rename(Share_not_in_slums=Share_Slums)
  # 
  # plottingData[1, "Share_not_in_slums"] <- 100
  # plottingData[2, "Share_not_in_slums"] <- 0
  #Define color scale
  library(scales)
  cluster_names <- 1:nMax
  cluster_colors <- setNames(scales::hue_pal()(nMax), cluster_names)
  
  colors_border <- col_numeric(palette = c("grey",cluster_colors[cluster],cluster_colors[cluster],cluster_colors[cluster],"black"), domain = NULL)(0:20)
  colors_fill <- adjustcolor(colors_border, alpha.f = 0.4)
  
  
  
  par(mai=c(0,0,0.5,0))
  par(bg=NA)
  png(filename = paste("Figures/WellBeingAnalysis/Radarchart",clustering,"_Cluster",i,".png"), width=500, height=500)
  radarchart(plottingData, vlabels = 4:55,
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

#SPI Dimension means ---------------------------------------------------------------
dimmeans <- clusteredData%>%
  mutate(BHN_mean = rowSums(across(5:22))/18,
         FWB_mean = rowSums(across(23:37))/15,
         OPP_mean = rowSums(across(38:56))/19)
#Join clusters 
#dimmeans <- left_join(dimmeans,clusters, by =c("Country","iso3"))

dimmeans_long <- dimmeans %>%
  pivot_longer(cols = c(BHN_mean, FWB_mean,OPP_mean), names_to = "variable", values_to = "value")

ggplot(dimmeans_long, aes(x = SPI_year, y = value, group=Country,color=factor(cluster))) +
  geom_line(alpha=1)+
  facet_grid(cluster~variable,labeller = label_value) +
  #theme_minimal() +
  #scale_color_manual(values = cluster_colors)+
  labs(title = "", x = "", y = "")+
  guides(color="none")

ggplot(dimmeans_long, aes(x = SPI_year, y = value, group=Country,color=factor(SPI_baseline))) +
  geom_line(alpha=1)+
  facet_grid(1~variable,labeller = label_value)


#DLS performance measure--------------------------------------------------------
#Idea to measure people brought out of material poverty during the study period


#Connect indicators with DLS dimensions
dimIndicators <- data.frame("Housing"="Share_not_in_Slums",
                            "Thermal comfort" = "Access_to_electricity",
                            "Food and Nutrition" = "Nutritional_deficiencies",
                            "Food prep/storage" = "Prevalence_of_cooking_with_coalperbiomass",
                            "Water" = "Safely_Managed_Drinking_Water",
                            "Sanitation" = "Safely_managed_saniation",
                            "Provision of health care" = "Universal_health_coverage",
                            "Education" = "Prim_School_Enroll",
                            "Social connectedness" = "Mobile_and_landline_telephone_subscriptions")%>%
  select(-Social.connectedness)


DLSdata <- clusteredData %>%
  mutate(Share_not_in_Slums = 100 - Share_Slums)%>%
  select(Country,iso3,SPI_year,cluster,.imp, all_of(c(as.matrix(dimIndicators))))%>%
  mutate(min_val = do.call(pmin, c(across(all_of(c(as.matrix(dimIndicators)))), na.rm = TRUE)))

DLSdata$which_min <- ""
for (i in 1:nrow(DLSdata)){
  DLSdata$which_min[i] <- c(as.matrix(dimIndicators))[which(DLSdata[i,c(as.matrix(dimIndicators))]==DLSdata$min_val[i])]
}

DLSdata_summarized<- DLSdata%>%
  group_by(SPI_year,Country,cluster)%>%
  mutate(min_est = mean(min_val))%>%
  filter(.imp==1)%>%
  left_join(GDP%>%select(c(iso3,GDP_PPP,SPI_year=Year)), by=c("iso3","SPI_year"))


ggplot(DLSdata_summarized%>%
         filter(cluster==2),aes(x = log(GDP_PPP), y=min_est, group = Country, color=cluster,shape=factor(which_min)))+
  geom_path()+
  geom_point()+
  scale_shape_manual(values = c(0:8))#+
  #facet_wrap(~which_min, scales="free_y")



#Plot country medians -----------------------------------------------------------
selected_vars <- rownames(batch_scaled)[batch_scaled[, clusterSelection] == 1]
set <- clusteredData %>%
  group_by(SPI_year, cluster,iso3) %>%
  dplyr::select(c("SPI_year", "cluster","iso3",
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
  select(1,2,3, ends_with("_median"))%>%
  pivot_longer(cols = ends_with("_median"),
               names_to = "Variable",
               values_to = "Median",
               names_pattern = "(.*)_median")
long_q1 <- set %>%
  select(1,2,3, ends_with("_q1"))%>%
  pivot_longer(cols = ends_with("_q1"),
               names_to = "Variable",
               values_to = "Q1",
               names_pattern = "(.*)_q1")
long_q3 <- set %>%
  select(1,2,3, ends_with("_q3"))%>%
  pivot_longer(cols = ends_with("_q3"),
               names_to = "Variable",
               values_to = "Q3",
               names_pattern = "(.*)_q3")
# Combine all
long_set <- long_median %>%
  left_join(long_q1, by = c("SPI_year", "cluster", "Variable","iso3")) %>%
  left_join(long_q3, by = c("SPI_year", "cluster", "Variable","iso3")) %>%
  mutate(Variable = factor(Variable, levels = variables))

#Filter
long_set_filtered <- long_set%>%
  #filter(cluster%in%c(4,8))
filter(cluster%in%c(1,7,11))
#filter(cluster%in%c(9,10))
#filter(cluster%in%c(2,3,6))
  #filter(cluster%in%c(1,2,3,4))
  #filter(cluster%in%c(5,6,7))


cluster_names <- 1:length(unique(clusters$cluster))
cluster_colors <- setNames(scales::hue_pal()(length(unique(clusters$cluster))), cluster_names)
ggplot(long_set_filtered, aes(x = SPI_year, y = Median, color = cluster, fill = cluster,group=iso3)) +
  geom_ribbon(aes(ymin = Q1, ymax = Q3), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8,alpha=0.4) +
  facet_wrap(~ Variable) +
  theme_minimal() +
  #scale_fill_manual(values = cluster_colors)+
  #scale_color_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Year", y = "Country Median (with IQR)")+
  theme(plot.title = element_text(size = 5))





#Compare clusterings in sankey plot-----------------------------------------------
AllClusters <- readRDS("clusterVariations_laglead_scaled.RDS")
library(ggsankey)
df_long <- AllClusters %>%
  make_long(Baseline, 'Few indicators_SPI_preferred')

ggplot(df_long, aes(x = x, 
                    next_x = next_x, 
                    node = node, 
                    next_node = next_node, 
                    fill = factor(node))) +
  geom_sankey(flow.alpha = 0.6, node.color = "gray30") +
  geom_sankey_label(aes(label = node),size = 3, color = "black") +
  theme_sankey(base_size = 16) +
  theme(legend.position = "none")+
  labs(title = "",
       x = "")



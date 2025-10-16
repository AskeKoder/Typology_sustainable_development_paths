#Cluster renaming
clusters <- readRDS("3_clusterVariations.RDS")


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
HDI_clust <- left_join(HDI, clusters,by = "Country")

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

#rename clusters by their mean HDI
for (name in (colnames(clusters)[3:11])){
    summary <- HDI_clust%>%
    drop_na(name)%>%
    group_by(!!sym(name))%>%
    summarise_at(3:7,mean)%>%
    mutate(Rank=rank(-`Human Development Index (HDI)`))
  
  #Rename clusters
  clusters[,name] <- summary$Rank[clusters[,name]]
}


#Join HDI data with clusters
HDI_clust <- left_join(HDI, clusters,by = "Country")
cluster = "DLSFew_coverage"

summary <- HDI_clust%>%
  drop_na(cluster)%>%
  group_by(!!sym(cluster))%>%
  summarise_at(3:7,CI)

for (i in 1:nrow(summary)){
  print(paste(summary[i,],collapse=" & "))
}


#Save ranked clusters
#saveRDS(clusters, "4_RankedClusters.RDS")

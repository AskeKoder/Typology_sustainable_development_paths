library(tidyr)
library(dplyr)
library(readxl)
library(countrycode)
library(ggplot2)


#Load All data ----------------------------------------------------------
WBdata <- read.csv("ImputedDataLag1Lead2_maxit30_scaled.csv")%>%
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

rownames(batch_scaled) <- colnames(WBdata)[5:(ncol(WBdata)-1)]
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

indicators <- rownames(batch_scaled[batch_scaled[,clusterSelection]>0,])

#Merge data sets
clusteredData <- WBdata %>%
  rename(iso3 = SPI_countrycode)%>%
  left_join(clusters , by= c("Country","iso3"))%>%
  data.frame()%>%
  mutate(Cluster = as.factor(cluster))%>%
  select(-cluster) 



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
ENVdata <- read_xlsx("Data_Environmental/GCSI_59a_Per_capita_2000-2020_05162025.xlsx") %>%
  pivot_wider(id_cols = c(Region_acronyms, Region_names, Year, Population),
              names_from = Indicator, values_from = "Value_per_Capita")%>%
  rename(Country=Region_names,
         iso3 = Region_acronyms)


ENVdata[ENVdata$Country=="Yugoslavia/Serbia (1991/1992)",2] <- "Serbia"
ENVdata[ENVdata$Country=="Yemen Arab Republic/Yemen (1990/1991)",2] <- "Yemen"

data_full <- ENVdata %>%
  #Standardize countries
  mutate(Country= countrycode(Country, origin="country.name",destination="country.name"))%>%
  
  
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
  filter(!is.na(Cluster))%>%
  #filter(.imp == 1)%>% #One iteration of imputations are chosen
  # dplyr::select(c(Country,ISO_Country,Year,Population,Cluster,ExtCluster,
  #                 GHG,Scarce_Water_Consumption, Biodiversity_Impact,
  #                 NOx,GDP_PPPcap,HDI))%>% #Select indicators
  arrange(Country,Year)%>%
  data.frame()

data <- data_full%>%
  filter(.imp==1) #Select one imputed set, when these are not needed


length(unique(data$Country))
length(unique(data$Cluster))
data %>%
  group_by(Country) %>%
  summarise(n_obs = n())%>%
  print(n=142)



titles <- c("GHG", "Blue Water Consumption", "Biodiversity Loss", "NOx")
units <- c("t CO2 e / cap", "m3 H20 e / cap", "PDF / cap", "kg NOx / cap")

#Violin plots --------------------------------------------------------------
ilist <- c("GHG",
           "Scarce_Water_Consumption",
           "Biodiversity_Impact",
           "NOx")
#Violin plot
data_norm <- data %>%
  rename_at(vars(ilist), ~ titles)%>%
  mutate(across(all_of(titles), ~ .x / weighted.mean(.x, Population, na.rm = TRUE), .names = "{.col}_norm"))

data_long <- data_norm %>%
  dplyr::select(Cluster, ends_with("_norm")) %>%
  pivot_longer(
    cols = ends_with("_norm"),
    names_to = "Indicator",
    values_to = "Value"
  ) %>%
  mutate(Indicator = gsub("_norm", "", Indicator))  # clean name
ggplot(data_long, aes(x = as.factor(Cluster), y = Value, fill = as.factor(Cluster))) +
  geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +  # Reference line at mean
  facet_wrap(~ Indicator,scale="free_y") +
  labs(
    x = "",
    y = "Relative Value (population wieghed mean = 1)",
    #title = "Normalized Indicator Distributions by Cluster"
  ) +
  #theme_minimal() +
  scale_y_continuous(trans='log10')+
  theme(legend.position = "none")


#Fit environmental models------------------------------------------
#GHG
data_wls <- data %>%
  group_by(Country) %>%
  mutate(sd_within_Country = sd(log(GHG), na.rm = TRUE)) %>%
  ungroup()%>%
  mutate(w = 1 / (sd_within_Country^2))

fitGHG.wls <- lm(log(GHG) ~ Year*Country, data = data_wls, weights = w)
# summary(fitGHG.wls)
# par(mfrow=c(2,2))
# plot(fitGHG.wls)
# anova(fitGHG.wls)
finalGHG <-fitGHG.wls

#Biodiversity
data_wls <- data %>%
  group_by(Country) %>%
  mutate(sd_within_Country = sd(log(Biodiversity_Impact), na.rm = TRUE)) %>%
  ungroup()%>%
  mutate(w = 1 / (sd_within_Country^2))

fitBio.wls <- lm(log(Biodiversity_Impact) ~ Year*Country, data = data_wls, weights = w)
# summary(fitBio.wls)
# par(mfrow=c(2,2))
# plot(fitBio.wls)
# AIC(fitBio.wls)
# AIC(fitBio2)
finalBio <-fitBio.wls

# Water consumption
data_wls <- data%>%
  group_by(Country) %>%
  mutate(sd_within_Country = sd(log(Scarce_Water_Consumption), na.rm = TRUE)) %>%
  ungroup()%>%
  mutate(w = 1 / (sd_within_Country^2))

fitWat.wls <- lm(log(Scarce_Water_Consumption) ~ Year*Country, data = data_wls, weights = w)
# summary(fitWat.wls)
# par(mfrow=c(2,2))
# plot(fitWat.wls)
# plot(fitWat2)
# AIC(fitWat.wls)
# AIC(fitWat2)
finalWat <-fitWat.wls


#Get means from models and plot ---------------------------------------------
#Compute geometric means
means <- data.frame("Year"=rep(2010,length(unique(data$Country))),
                    "Country"=unique(data$Country),
                    "iso3"=unique(data$iso3),
                    "Cluster"=rep("na",length(unique(data$Country))))
#Assign clusters to means 
for (i in 1:length(means$Country)){
  co <- means$Country[i]
  cl <- unique(data[data$Country==co,"Cluster"])
  means$Cluster[i] <- as.character(cl)
}

#Ensure correct format of variables and order of means 
means$Cluster <- factor(means$Cluster, levels=c(1,2,3,4,5,6,7,8,9,10,11))
means$Country <- factor(means$Country, levels = means$Country[order(means$Cluster)])
#means$Country <- factor(means$Country, levels = means$Country[order(means$ExtCluster)])


#Add predicted means and confidence intervals - note we're converting back from the log domain
means <- cbind(means,exp(predict(finalGHG,newdata = means,interval = "confidence")))
colnames(means)[5:7] <- c("fitGHG","lwrGHG","uprGHG")


means <- cbind(means,exp(predict(finalBio,newdata = means,interval = "confidence")))
colnames(means)[8:10] <- c("fitBio","lwrBio","uprBio")

means <- cbind(means,exp(predict(finalWat,newdata = means,interval = "confidence")))
colnames(means)[11:13] <- c("fitWat","lwrWat","uprWat")



#Plot Estimated means of the world
#GHG
cluster_names <- 1:11
cluster_colors <- setNames(scales::hue_pal()(11), cluster_names)
#Order countries
means$iso3 <- factor(means$iso3, levels = means$iso3[order(means$Cluster,means$fitGHG)])
ggplot(means, aes(x =iso3, y = fitGHG, fill = Cluster)) +
  geom_col(width = 0.6) +
  #geom_hline(yintercept = 0.4, color = "black", linetype = "solid")+
  #geom_hline(yintercept = 3.5, color = "black", linetype = "dashed")+
  geom_errorbar(aes(ymin = lwrGHG, ymax = uprGHG), width = 0.2) +
  geom_point(data%>%
               group_by(iso3) %>%
               summarise(mean_value = mean(GHG, na.rm = TRUE)),
             shape="-",
             mapping = aes(x = iso3, y = mean_value, fill = iso3),
             alpha=1,size=4)+
  #geom_hline(yintercept=3)+
  theme_minimal() +
  labs(
    title = "",
    y = "t CO2e/cap",
    x = ""
  ) +
  scale_fill_manual(values = cluster_colors)+
  #scale_color_manual(values = cluster_colors) +
  #guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10))+
  theme(legend.position="bottom")+
  guides(fill = guide_legend(ncol = 11))

#Biodiversity
ggplot(means, aes(x =iso3, y = fitBio, fill = Cluster)) +
  geom_col(width = 0.6) +
  #geom_hline(yintercept = 0.4, color = "black", linetype = "solid")+
  #geom_hline(yintercept = 3.5, color = "black", linetype = "dashed")+
  geom_errorbar(aes(ymin = lwrBio, ymax = uprBio), width = 0.2) +
  geom_point(data%>%
               group_by(iso3) %>%
               summarise(mean_value = mean(Biodiversity_Impact, na.rm = TRUE)),
             shape="-",
             mapping = aes(x = iso3, y = mean_value, fill = iso3),
             alpha=1,size=4)+
  #geom_hline(yintercept=3)+
  theme_minimal() +
  labs(
    title = "",
    y = "PDF / cap",
    x = ""
  ) +
  scale_fill_manual(values = cluster_colors)+
  #scale_color_manual(values = cluster_colors) +
  #guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10))+
  theme(legend.position="bottom")+
  guides(fill = guide_legend(ncol = 11))

#Water
ggplot(means, aes(x =iso3, y = fitWat, fill = Cluster)) +
  geom_col(width = 0.6) +
  #geom_hline(yintercept = 0.4, color = "black", linetype = "solid")+
  #geom_hline(yintercept = 3.5, color = "black", linetype = "dashed")+
  geom_errorbar(aes(ymin = lwrWat, ymax = uprWat), width = 0.2) +
  geom_point(data%>%
               group_by(iso3) %>%
               summarise(mean_value = mean(Scarce_Water_Consumption, na.rm = TRUE)),
             shape="-",
             mapping = aes(x = iso3, y = mean_value, fill = iso3),
             alpha=1,size=4)+
  #geom_hline(yintercept=3)+
  theme_minimal() +
  labs(
    title = "",
    y = "m3 H20 e / cap",
    x = ""
  ) +
  scale_fill_manual(values = cluster_colors)+
  #scale_color_manual(values = cluster_colors) +
  #guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10))+
  theme(legend.position="bottom")+
  guides(fill = guide_legend(ncol = 11))




#Extract growth rates and plot-------------------------------------------
model <- finalWat

#Quick fix due to no GDP data for cuba
nCountries<-length(grep("Year:",names(model$coefficients)))+1

#Extract ML estimates
theta_hat <- rep(0,nCountries)
theta_hat[1] <- model$coefficients[2] #Assign intercept slope to first country
theta_hat[2:nCountries] <- model$coefficients[grep("Year:",names(model$coefficients))] + model$coefficients[2] #intercept + country slope


#(b1+b2) Variance estimate #year:Afghanistan is the intercept
var <- diag(vcov(model))[grep("Year",rownames(vcov(model)))] #Extract diagonal of covariance matrix to get variance of parameters
cov <- vcov(model)[grep("Year:",rownames(vcov(model))),"Year"] #Covariances between slope parameters and intercept slope

#Assign variances to countries:
theta_hat_var <- sapply(1:length(unique(data$Country)),
                        function(x){ if (x==1){
                          var[1]
                        } else{
                          var[1]+var[x]+ 2*cov[x-1]
                        }
                        })

#Name the parameters appropriately
names(theta_hat_var) <- names(var)
names(theta_hat_var)[1] <- "Afghanistan"

#Compute standard errors
theta_hat_SE <- sqrt(theta_hat_var)

#Compute confidence intervals in log domain
df <- model$df.residual #Degrees of freedom
alpha <- 0.05
t_crit <- qt(1 - alpha/2, df)
log_rate <- cbind(theta_hat,theta_hat-t_crit*theta_hat_SE,theta_hat+t_crit*theta_hat_SE)

#Clean up data frame
log_rate <- data.frame(log_rate)
rownames(log_rate) <- gsub("Year:Country","",rownames(log_rate))
colnames(log_rate)<-c("Rate","lwrRate","uprRate")
log_rate$Country <- rownames(log_rate)

#Convert to regular domain
rate <- log_rate%>%
  mutate(Rate=exp(Rate), lwrRate=exp(lwrRate),uprRate=exp(uprRate))

#Convert to percent change
rate[,c(1,2,3)] <- (rate[,c(1,2,3)]-1)*100

#Join with clusters
rates <- left_join(rate,clusters,by="Country")
rates$Country <- factor(rates$Country, levels = rates$Country[order(rates$cluster,rates$Rate)])
rates$iso3 <- factor(rates$iso3, levels = rates$iso3[order(rates$cluster,rates$Rate)])

#Plot
ggplot(rates, aes(x =iso3, y = Rate, fill = factor(cluster))) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lwrRate, ymax = uprRate), width = 0.2) +
  theme_minimal() +
  labs(
    title = "",
    y = paste("yearly %-rate of change in", var),
    x = "Country"
  ) +
  #scale_fill_manual(values = cluster_colors)+
  #scale_color_manual(values = cluster_colors) +
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Test mean/rate scatter plot
test <- left_join(means,rates, by="iso3")
factor <- c(0.2,0.1)
scale_factor <- 0.2*max(test$fitWat, na.rm=TRUE) / max(test$Rate, na.rm=TRUE)
test$iso3 <- factor(test$iso3,levels=test$iso3[order(test$cluster,test$fitWat)])
ggplot(test, aes(x = iso3)) +
  geom_col(aes(y = fitWat, fill = factor(cluster)), width = 0.6) +
  geom_errorbar(aes(ymin = lwrWat, ymax = uprWat), width = 0.2) +
  # Scale y2 up to match y1 axis
  geom_point(aes(y = Rate * scale_factor), color = "black",size=0.9) +
  geom_errorbar(aes(ymin = lwrRate* scale_factor, ymax = uprRate* scale_factor), width = 0.2)+
  scale_y_continuous(
    name = "m3 H2O e /cap",
    sec.axis = sec_axis(~ . / scale_factor, name = "Yearly %-rate of change")
  ) +
  theme_minimal() +
  labs(x = "Country") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none"
  )


#Fitting multivariable cluster model-------------------------------------------------
vars <- c("Cluster","logGDPcap","logGHG","logBio","logWat")
data_mvmodel <- data %>%
  mutate(logGHG = scale(log(GHG)),
         logBio = scale(log(Biodiversity_Impact)),
         logWat = scale(log(Scarce_Water_Consumption)),
         logGDPcap = scale(log(GDP_PPPcap)))
  #filter(Cluster!=11)%>%
  #na.omit()

pairs(data_mvmodel[,vars],col=data_mvmodel$Cluster)
print(psych::corTest(data_mvmodel[,c("logGHG","logBio","logWat")]),short=FALSE)

MvModel <- lm(cbind(logGHG,
                    logBio,
                    logWat)~Cluster*log(GDP_PPPcap),
                           data=data_mvmodel)
summary(MvModel)
car::Anova(MvModel, type="III")
summary.aov(MvModel)

pairwise(MvModel, "Cluster")
#Visualize model

# Basic HE plot
heplots::heplot(MvModel,
                variables=c("logGHG", "logBio"),
                type = "III")

# Add group mean ellipses
heplots::heplot(MvModel, variables=c("logGHG", "logBio"), 
                , terms = "Cluster",
       fill = TRUE, fill.alpha = 0.1)

pairs(MvModel)

can <- candisc::candisc(MvModel, term = "Cluster", type="III",
                        ellipse = TRUE,
                        ellipse.fill = TRUE,          # fill the ellipses
                        ellipse.fill.alpha = 1,     # transparency (0 = transparent, 1 = opaque)
                        ellipse.line.lwd = 2,         # make ellipse lines thicker
                        ellipse.line.lty = 1,         # solid lines
                        #col = c("tomato", "steelblue", "forestgreen"), # group colors
                        pch = 19,                     # filled points
                        level = 0.95,
                        )

# Plot canonical dimensions
plot(can)

scores <- as.data.frame(can$scores)

ggplot(scores, aes(x = Can1, y = Can2, fill = Cluster, color = Cluster)) +
  stat_ellipse(type = "norm", level = 0.95, geom = "polygon", alpha = 0.5) +
  geom_path(shape = 21, size = 3, color = "black",alpha=0.2,group=Country) +
  theme_minimal(base_size = 14) +
  # scale_fill_brewer(palette = "Set1") +
  # scale_color_brewer(palette = "Set1") +
  labs(title = "Canonical Discriminant Analysis (Species separation)",
       x = "Canonical Dimension 1",
       y = "Canonical Dimension 2")


#Emmeans and Emtrends----------------------------------------------
library(emmeans)
library(lme4)

#Fit individual models
mod_logGHG <- lmer(logGHG ~ Cluster * logGDPcap + (logGDPcap|Country), data = data_mvmodel,REML = FALSE)
mod_logBio <- lmer(logBio ~ Cluster * logGDPcap+ (logGDPcap|Country), data = data_mvmodel,REML = FALSE)
mod_logWat <- lmer(logWat ~ Cluster * logGDPcap + (logGDPcap|Country), data = data_mvmodel,REML = FALSE)
model_names <- c("mod_logGHG", "mod_logBio", "mod_logWat")

summary(mod_logGHG)
car::Anova(mod_logGHG,type="III")
plot(ranef(mod_logGHG)$Country[,2])

summary(mod_logBio)
car::Anova(mod_logBio,type="III")

summary(mod_logWat)
car::Anova(mod_logWat,type="III")
mod_logWat <- update(mod_logWat,~.-Cluster:logGDPcap)
car::Anova(mod_logWat,type="III")
mod_logWat <- update(mod_logWat,~.-logGDPcap)
car::Anova(mod_logWat,type="III")

#Refit significant models using REML for unbiased estimates
mod_logGHG <- lmer(logGHG ~ Cluster * logGDPcap + (logGDPcap|Country), data = data_mvmodel,REML = TRUE)
mod_logBio <- lmer(logBio ~ Cluster * logGDPcap+ (logGDPcap|Country), data = data_mvmodel,REML = TRUE)
mod_logWat <- lmer(logWat ~ Cluster + (logGDPcap|Country), data = data_mvmodel,REML = TRUE)
summary(mod_logGHG)

# Initialize lists
emmeans_results <- list()
emtrends_results <- list()

emmeans_df <- list()
emtrends_df <- list()

# Loop through models
for (model_name in model_names) {
  model_obj <- get(model_name)
  
  # Compute emmeans and emtrends
  if (model_name!="mod_logWat"){
    emmeans_res <- emmeans(model_obj, ~ Cluster, by="logGDPcap",pbkrtest.limit = 4000)
  }else{
    emmeans_res <- emmeans(model_obj, ~ Cluster,pbkrtest.limit = 4000)
  }
  # Store results
  emmeans_results[[model_name]] <- emmeans_res
  
  
  # Convert to data frames and add a model column
  emmeans_df[[model_name]] <- as.data.frame(emmeans_res)
  emmeans_df[[model_name]]$model <- model_name
  
  
  if (model_name!="mod_logWat"){
    emtrends_res <- emtrends(model_obj, ~ Cluster, var = "logGDPcap", by="Cluster",pbkrtest.limit = 4000)
    emtrends_results[[model_name]] <- emtrends_res
    emtrends_df[[model_name]] <- as.data.frame(emtrends_res)
    emtrends_df[[model_name]]$model <- model_name
  }
  
}

# Combine all into single data frames
emmeans_all <- do.call(rbind, emmeans_df)
emtrends_all <- do.call(rbind, emtrends_df)

emmeans_all$model= factor(emmeans_all$model, levels=model_names)
emtrends_all$model= factor(emtrends_all$model, levels=model_names)

#plots
p1 <- ggplot(emmeans_all, aes(x = Cluster, y = emmean,color=Cluster)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2,lwd=0.7) +
  theme_minimal() +
  labs(title = "Estimated Marginal Means", y = "EM Mean (log-domain)", x = "") +
  facet_wrap(~ model, scales = "free_y")+
  guides(color="none")

# emtrends plot faceted by model
p2 <- ggplot(emtrends_all, aes(x = Cluster, y = logGDPcap.trend,color=Cluster)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.4,lwd=0.7) +
  theme_minimal() +
  labs(title = "", y = "d log(GHG) / d log(GDPcap)", x = "") +
  facet_wrap(~ model, scales = "free_y")+
  guides(color="none")

# Combine with patchwork
library(patchwork)
(p1 / p2)
p2


#multilevel model with brms-------------------------
library(MCMCglmm)
data_mvmodel <- data_mvmodel%>%
  na.omit()
MvMlmodel <- MCMCglmm::MCMCglmm(cbind(logGHG, logBio,logWat)
                                ~ trait:Cluster:logGDPcap + trait:Cluster+  trait:logGDPcap -1, 
                                random = ~us(1+logGDPcap):Country,
                                rcov = ~us(trait):units,
                                data=data_mvmodel,
                                family=rep("gaussian",3))
summary(MvMlmodel)



#Multilevel long format------------------------------------------------
data_long <- data_mvmodel %>%
  pivot_longer(cols = c(logGHG, logBio, logWat),
               names_to = "variable",
               values_to = "value")

library(lme4)
fit <- lmer(value~variable:(logGDPcap*Cluster) + (variable:logGDPcap|Country),
     data=data_long)
summary(fit)
car::Anova(fit, type="III")


#Compute index of change -----------------------------------------------
index_of_change <- data_full %>%
  group_by(Country,.imp) %>%
  mutate(across(all_of(indicators),
    ~ (last(.) - first(.)),
    .names = "{.col}_change"
  ))

weights <- index_of_change%>%
  group_by(cluster)%>%
  summarise_at(all_of(paste0(indicators,"_change")),mean)

# index_of_change <- index_of_change %>%
#   rowwise() %>%
#   mutate(
#     mean_change = weighted.mean(
#       c_across(ends_with("_change")),
#       w = weights,
#       na.rm = TRUE
#     )
#   ) %>%
#   ungroup()%>%
#   filter(Year==2020) #Select a random year

index_of_change <- index_of_change %>%
  left_join(weights, by = "cluster", suffix = c("", "_w"))

index_of_change <- index_of_change %>%
  filter(Year == 2020)%>% #Select random year as all years are the same
  rowwise() %>%
  mutate(
    mean_change = weighted.mean(
      x = c_across(ends_with("_change") & !ends_with("_w")),  # the indicator changes
      w = 1/(c_across(ends_with("_change_w"))+0.000001),                   # the weights
      na.rm = TRUE
    )
  ) %>%
  ungroup() 

index_of_change$iso3 <- factor(index_of_change$iso3,levels=unique(index_of_change$iso3[order(index_of_change$cluster)]))
ggplot(index_of_change, aes(x=iso3,y=mean_change,color=factor(cluster)))+
  geom_point()+
  geom_jitter(width=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

joined <- left_join(index_of_change,rates,by=c("Country","iso3", "cluster"))

fit <- lm(mean_change~Rate, data=joined)
summary(fit)

joined$iso3 <- factor(joined$iso3, levels=unique(joined$iso3[order(joined$cluster)]))
ggplot(joined,aes(x=Rate,y=mean_change,color=factor(cluster)))+
  geom_point()

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


DLSdata <- data %>%
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


#Fit DLS Mixed MODEL ---------------------------------------------
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

#Fit DLS S curve-------------------------------
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

cluster_names <- 1:7
cluster_colors <- setNames(scales::hue_pal()(length(unique(clusters$cluster))), cluster_names)

ggplot()+
  geom_path(data=results%>%
              filter(Year %in% c(2000,2020)), aes(x=GHG,y=x,color=factor(cluster),group=Country),
            alpha=0.4)+
  geom_point(data=results%>%
               filter(Year %in% c(2000,2020)), aes(x=GHG,y=x,color=factor(cluster),shape=factor(Year)))+
  geom_line(data=results, aes(x=GHG,y=pred))+
  facet_wrap(~variable)+
  #geom_vline(xintercept=3, lty=2)+ #at 3 ton co2/cap
  scale_color_manual(values = cluster_colors)+
  coord_trans(x = "log10")+
  labs(x="t CO_2 cap", y="%")
  #theme_minimal()


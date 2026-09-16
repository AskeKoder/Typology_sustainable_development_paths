#Librarires
library(dplyr)
library(ggplot2)
library(nnet)
library(rnaturalearth)

#load data
folder <- "Data_Colonial/"

#settler mortality
setMort <- read.csv(paste0(folder,"SettlerMortality.csv"), sep=";", dec=",")%>%
  select(cname, ccodealp, ajr_settmort)%>%
  unique()%>%
  rename(Country = cname,
         iso3 = ccodealp,
         settmort = ajr_settmort)

#colonizer
col <- read.csv(paste0(folder,"colonies.csv"))%>%
  filter(!(European.colonial.power..grouped. %in% c("z. Multiple colonizers",
                                                 "zz. Colonizer",
                                                 "zzz. Not colonized",
                                                 "zzzz. No longer colonized")))%>%
  group_by(Entity, Code,European.colonial.power..grouped.) %>%
  summarise(n_timepoints = n()) %>%
  slice_max(n_timepoints, n = 1, with_ties = FALSE)%>%
  rename(Country = Entity,
         iso3 = Code)%>%
  select(iso3,colonizer = European.colonial.power..grouped.)


#Furthermore identify countries that were never colonized
col_never <- read.csv(paste0(folder,"colonies.csv"))%>%
  filter(!(European.colonial.power..grouped. %in% c("z. Multiple colonizers",
                                                    "zz. Colonizer",
                                                    "zzzz. No longer colonized")))%>%
  group_by(Entity, Code,European.colonial.power..grouped.)%>%
  summarise(n_timepoints = n())%>%
  add_count(Entity)%>%
  filter((n==1) & (European.colonial.power..grouped. == "zzz. Not colonized"))


#primary school enrollment 1900
prienr1900 <- read.csv(paste0(folder,"prienr.csv"))%>%
  filter(Year == 1900)%>%
  filter(!(Code %in% c("","OWID_WRL")))%>%
  group_by(Code)%>%
  slice_max(Year)%>%
  ungroup()%>%
  rename(Country = Entity,
         iso3 = Code)%>%
  select(iso3,prienr1900 = Combined.total.net.enrolment.rate..primary..both.sexes )

#Join data 
df <- merge(setMort, col, by=c("iso3"), all=TRUE) %>%
  merge(prienr1900, by=c("iso3"), all=TRUE)%>%
  select(iso3,colonizer, settmort, prienr1900)
df[df$iso3%in%col_never$Code,"colonizer"] <- "Not colonized" #Add the option to not have been colonized

#Load clusters
clusterVariations <- readRDS("clustervariations_laglead_scaled.RDS")%>%
  select(iso3=SPI_countrycode, cluster = 'Few indicators_SPI_preferred')%>%
  mutate(cluster= factor(cluster))




#Group colonizers with less than 5 observations
df <- df %>%
  add_count(colonizer) %>%
  mutate(colonizer = if_else(n < 10, "Other", as.character(colonizer))) %>%
  select(-n)%>%
  mutate(colonizer=factor(colonizer))

#Visualize colonizers
world <- ne_countries(scale = "medium", returnclass = "sf",continent = c("south america","oceania","north america", "asia","europe","africa"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'SDZ',"SDN"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'PN1',"PNG"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'PR1',"PRT"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'SSD',"SSD*"))
colnames(world)[57] <- "iso3"

#Append clusters to world data
world <- left_join(world, df, 
                   by = "iso3")
ggplot() +
  geom_sf(data = world, aes(fill = factor(colonizer)), color = "white",size=0.5)+
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())+
  labs(title="Longest lasting colonizer")+
  guides(fill=guide_legend(title="Colonizer",ncol=1))

#Fit ordinal model---------------------------------------------
library(ordinal)
df_model$cluster <- factor(df_model$cluster, levels=c(10,4,1,8,9,7,3,2,5,11,6))

fitORD <- clm(cluster~colonizer*settmort,
    data=df_model)
summary(fitORD)
car::Anova(fitORD,type="III")


library(ggeffects)

# Generate predictions
preds <- ggpredict(fitORD, terms = c("settmort","colonizer"))

# Plot predicted probabilities
plot(preds)


#Fit model------------------------------------------------------
df_model <- df %>%
  left_join(clusterVariations, by="iso3")%>%
  filter(colonizer!="Not colonized")%>%
  mutate(colonizer=factor(colonizer)) #Filter out countries that have not been colonized

#Overview
table(df_model$cluster,df_model$colonizer)

#Model all data
fit <- multinom(cluster ~ colonizer,
                data=df_model)
summary(fit)           

z <- summary(fit)$coefficients/ summary(fit)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)
car::Anova(fit,type=2)



#Aggregate model------------------------------------------------------
#Aggregate clusters not well represented in colonial data
df_agg <- df_model %>%
  mutate(cluster = if_else(cluster %in%c("2","3","5"), "2,3,5", as.character(cluster)))%>%
  mutate(cluster= factor(cluster),
         colonizer= factor(colonizer))%>%
  filter(!is.na(colonizer))

table(df_agg$cluster, df_agg$colonizer)

fit_agg <- multinom(cluster~colonizer, data=df_agg)
summary(fit_agg) #Parameters do not look good - too many clusters

z <- summary(fit_agg)$coefficients/ summary(fit_agg)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)
car::Anova(fit_agg,type=2)  

#How about settler mortality 
fit_agg_settmort <- multinom(cluster~settmort, data=df_agg)
summary(fit_agg_settmort)
z <- summary(fit_agg_settmort)$coefficients/ summary(fit_agg_settmort)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)
car::Anova(fit_agg_settmort,type=2) #Settler mortality is a highly signifcicant predictor

#What if we include primary school enrollment
car::Anova(update(fit_agg_settmort,~.+settmort*prienr1900),type=2) #Both settler mortality and primary school enrollment seems significant
fit_agg <- update(fit_agg_settmort,~.+prienr1900,maxit=500)
z <-  summary(fit_agg)$coefficients/ summary(fit_agg)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)
car::Anova(fit_agg,type=2) #Both are significant


#Continuous graph
preds <- data.frame(ggeffects::ggemmeans(fit_agg, terms=c("prienr1900")))
ggplot(preds, aes(x=x, y=predicted,color = response.level))+
  geom_line()+
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = response.level),
    alpha = 0.2,
    color = NA
  )+
  labs(x="Primary school enrollment",
       y="Probability",
       fill="Colonizer")+
  theme_minimal()


preds <- data.frame(ggeffects::ggemmeans(fit_agg, terms=c("settmort")))
ggplot(preds, aes(x=x, y=predicted,color = response.level))+
  geom_line()+
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = response.level),
    alpha = 0.2,
    color = NA
  )+
  labs(x="log(Settler mortality)",
       y="Probability",
       fill="Cluster")+
  theme_minimal()+
  guides(color="none")


#Barchart for colonizer
preds <- data.frame(ggeffects::ggemmeans(fit_agg, terms="colonizer"))
ggplot(preds, aes(x=x, y=predicted,fill = response.level))+
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    color = "black",
    width = 0.2,
    position = position_dodge(width = 0.9)
  )+
  labs(x="Colonizer",
       y="Probability",
       fill="Colonizer")+
  theme_minimal()

#Including settler mortality
df_agg2 <- df_agg %>% 
  filter(!is.na(settmort))

fit_agg2 <- multinom(cluster~colonizer+settmort, data=df_agg2)
summary(fit_agg2)

z <- summary(fit_agg2)$coefficients/ summary(fit_agg2)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)
car::Anova(fit_agg2,type=2)  
lmtest::lrtest(update(fit_agg2,~1),fit_agg2)


preds <- data.frame(ggeffects::ggemmeans(fit_agg2, terms=~settmort+colonizer))
ggplot(preds, aes(x=x, y=predicted, color=response.level))+
  geom_line(show.legend = FALSE) +
  geom_point(data=df_agg2%>%
               select(response.level=cluster,
                      group=colonizer,
                      x=settmort)%>%
               na.omit(),
             aes(x=x,y=0, color=response.level),
             shape=4,size=3, stroke = 1.5,alpha=0.7,
             position = position_jitter(width = 0.1, height = 0))+
  facet_wrap(~group)+
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = response.level),
    alpha = 0.2,
    color = NA
  )+
  labs(x="log(Settler mortality)",
       y="Probability",
       fill="Cluster",
       color="Observed clusters")+
  theme_minimal()


#Include primary school enrollment?
fit_agg3 <- update(fit_agg2,~.+prienr1900,maxit=200)
summary(fit_agg3)
z <- summary(fit_agg3)$coefficients/ summary(fit_agg3)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)
car::Anova(fit_agg3,type=2) 

preds <- data.frame(ggeffects::ggemmeans(fit_agg3, terms=~settmort+prienr1900+colonizer))
lmtest::lrtest(fit_agg3,update(fit_agg2,data=df_agg2%>%filter(!is.na(prienr1900))))

#Limt the analysis to clusters represented in africa 
df_filtered <- df_model %>%
  filter(cluster %in% c(1,2,6,7,8))%>%
  mutate(cluster= factor(cluster),
         colonizer= factor(colonizer))%>%
  filter(!is.na(colonizer))
  # filter(!is.na(settmort))

fit_filtered <- multinom(cluster ~ colonizer,
                data=df_filtered,
                maxit=700)

summary(fit_filtered)
z <- summary(fit_filtered)$coefficients/ summary(fit_filtered)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)
car::Anova(fit_filtered,type=2)  

fit_null <- multinom(cluster~1,
                     data=df_filtered)
lmtest::lrtest(fit_null,fit_filtered)

preds <- data.frame(ggeffects::ggemmeans(fit_filtered, terms="colonizer"))

ggplot(preds, aes(x=x, y=predicted,fill = response.level))+
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    color = "black",
    width = 0.2,
    position = position_dodge(width = 0.9)
  )+
  labs(x="Cluster",
       y="Probability",
       fill="Colonizer")+
  theme_minimal()

#Add settler mortality as a variable
df_filtered2 <- df_model %>%
  filter(cluster %in% c(1,2,6,7))%>%
  mutate(cluster= factor(cluster),
         colonizer= factor(colonizer))%>%
  filter(!is.na(colonizer))%>%
  filter(!is.na(settmort))

fit_filtered2 <- multinom(cluster ~ colonizer+settmort,
                         data=df_filtered2,
                         maxit=700)

summary(fit_filtered2)
z <- summary(fit_filtered2)$coefficients/ summary(fit_filtered2)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)
car::Anova(fit_filtered2,type=2)  

fit_null <- multinom(cluster~1,
                     data=df_filtered2)
lmtest::lrtest(fit_null,fit_filtered)



#Add primary school enrollment as a variable
df_filtered2 <- df_model %>%
  filter(cluster %in% c(1,2,6,7))%>%
  mutate(cluster= factor(cluster),
         colonizer= factor(colonizer))%>%
  filter(!is.na(colonizer))%>%
  filter(!is.na(settmort))

fit_filtered2 <- multinom(cluster ~ colonizer+settmort,
                          data=df_filtered2,
                          maxit=700)

summary(fit_filtered2)
z <- summary(fit_filtered2)$coefficients/ summary(fit_filtered2)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)
car::Anova(fit_filtered2,type=2)  

fit_null <- multinom(cluster~1,
                     data=df_filtered2)
lmtest::lrtest(fit_null,fit_filtered)


##Visualize model
# Make prediction grid
newdata <- expand.grid(
  colonizer = levels(df_filtered2$colonizer),
  settmort = seq(min(df_filtered2$settmort,na.rm=TRUE), max(df_filtered$settmort,na.rm=TRUE), length.out = 50)
)

# Get predicted probabilities
preds <- predict(fit_filtered2, newdata = newdata, type = "probs")

# Convert to long format
pred_df <- cbind(newdata, preds) %>%
  tidyr::pivot_longer(cols = -c(colonizer,settmort),
                      names_to = "cluster", values_to = "probability")%>%
  unique()

# Plot
ggplot(pred_df, aes(x = settmort, y = probability, color = cluster)) +
  geom_line() +
  geom_point(data=df_filtered2, aes(x=settmort, y=1, fill=factor(cluster)), color="black",pch=21, alpha=0.5,size=3)+
  facet_wrap(~colonizer) +
  theme_minimal() +
  labs(title = "Predicted Probabilities by Colonizer, Settler mortality and school enrollment (1900)")


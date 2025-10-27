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
clusterVariations <- readRDS("4_RankedClusters.RDS")%>%
  select(iso3=SPI_countrycode, cluster = 'DLSFew_coverage')%>%
  mutate(cluster= factor(cluster))

#Group colonizers with less than 10 observations
df <- df %>%
  add_count(colonizer) %>%
  mutate(colonizer = if_else(n < 10, "Other", as.character(colonizer))) %>%
  select(-n)%>%
  mutate(colonizer=factor(colonizer))%>%
  merge(clusterVariations, by="iso3",all=TRUE)



#Load GDP data
raw_GDP <- read.csv("Data_WellBeing/GDPpercap PPP 2001 international world bank.csv", header = FALSE, stringsAsFactors = FALSE)
GDP <- raw_GDP[-c(1:3), ] #Remove metadata
colnames(GDP) <- raw_GDP[3, ] #Set colnames


GDP <- GDP %>%
  tidyr::pivot_longer(
    cols = matches("^\\d{4}$"),  # four numbers(\\d{4}) between ^start and $end of string 
    names_to = "Year",
    values_to = "GDP_PPP_current_international_dollars"
  ) %>%
  dplyr::select(
    Country = `Country Name`,
    ISO_Country = `Country Code`,
    Year,
    GDP_PPP_current_international_dollars)%>%
  filter(Year == 2020)%>%
  mutate(Country = as.factor(Country),
         ISO_Country = as.factor(ISO_Country),
         Year = as.numeric(Year))%>%
  rename(GDP_PPP = GDP_PPP_current_international_dollars)

summary(GDP)
str(GDP)

#Standardize names 
library(countrycode)
GDP$Country_std  <- countrycode(GDP$Country, origin="country.name",destination="country.name")
GDP$iso3  <- countrycode(GDP$Country, origin="country.name",destination="iso3c")
#Print where the renaming failed
print(unique(GDP[is.na(GDP$Country_std), "Country"]),n=100)  #Only aggregated countries failed. No problem


#Initial asssessment---------------------------------------------------------------------------
df_dual <- df %>%
  left_join(GDP,by="iso3")%>%
  filter(!is.na(colonizer))%>%
  filter(!is.na(cluster))%>%
  filter(!is.na(settmort))%>%
  group_by(cluster)%>%
  mutate(across(c(GDP_PPP,settmort), mean, na.rm = TRUE))%>%
  add_count()

coeff <- 1.7
ggplot(df_dual, aes(x = cluster)) +
  geom_col(aes(y = settmort/n, fill = colonizer)) +
  geom_text(aes(y = (settmort + 0.2), label=n))+
  geom_line(aes(y = log(GDP_PPP)/coeff, group = 1, color = "log(GDP/cap)")) +
  scale_y_continuous(
    name = "log(settmort)",
    sec.axis = sec_axis(~.*coeff, name = "log(GDP/cap)")
  )+
  scale_color_manual(
    name = "",
    values = c("log(GDP/cap)" = "black")
  )+
  theme_minimal()+
  labs(title="", x="Cluster")

#Fit model------------------------------------------------------
df_model <- df %>%
  filter(colonizer!="Not colonized")%>%
  mutate(colonizer=factor(colonizer))%>% #Filter out countries that have not been colonized
  na.omit()

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




fit2 <- multinom(cluster ~ colonizer+ settmort,
                 data=df_model)
summary(fit2)
z <- summary(fit2)$coefficients/ summary(fit2)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)
car::Anova(fit2,type=2)



fit3 <- multinom(cluster ~ colonizer+ settmort+prienr1900,
                 data=df_model,maxit=500)
summary(fit3)
z <- summary(fit3)$coefficients/ summary(fit3)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)
car::Anova(fit3,type=2)

fit4 <- multinom(cluster ~ settmort+prienr1900,
                 data=df_model,maxit=500)
summary(fit4)
z <- summary(fit4)$coefficients/ summary(fit4)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)
car::Anova(fit4,type=2)

#Ensure that datasets are comparable
df_compare <- df_model%>%
  na.omit()
AIC(multinom(cluster ~ settmort+prienr1900,
         data=df_compare,maxit=500),multinom(cluster ~ colonizer+ settmort+prienr1900,
                                           data=df_compare,maxit=500) )
#Without colonizer is the best

# Make prediction grid
newdata <- expand.grid(
  prienr1900 = seq(min(df_model$prienr1900,na.rm=TRUE), max(df_model$prienr1900,na.rm=TRUE), length.out = 50),
  settmort = seq(min(df_model$settmort,na.rm=TRUE), max(df_model$settmort,na.rm=TRUE), length.out = 50)
)

# Get predicted probabilities
preds <- predict(fit4, newdata = newdata, type = "probs")

# Convert to long format
pred_df <- cbind(newdata, preds) %>%
  tidyr::pivot_longer(cols = -c(settmort,prienr1900),
                      names_to = "cluster", values_to = "probability")%>%
  unique()
# 
# preds <- data.frame(ggeffects::ggemmeans(fit4, terms=c("prienr1900","settmort [3:7 by=0.4]")))
# ggplot(preds, aes(x=x, y=predicted,color = response.level))+
#   geom_line()+
#   geom_ribbon(
#     aes(ymin = conf.low, ymax = conf.high, fill = response.level),
#     alpha = 0.2,
#     color = NA
#   )+
#   labs(x="Primary school enrollment",
#        y="Probability")+
#   theme_minimal()+
#   facet_wrap(~group)


library(ggplot2)



ggplot(pred_df, aes(x = prienr1900, y = settmort, z = probability)) +
  geom_contour_filled() +
  geom_point(data=df_model, aes(x=prienr1900, y=settmort,z=0,shape=colonizer),color="salmon")+
  labs(x="Primary school enrollment (1900)",
       y="log(settler mortality)",
       fill="Predicted probability",
       shape = "Observation and \ncolonizer identity",
       color=NULL)+
  theme_minimal() +
  facet_wrap(~ factor(cluster, levels=1:11))
# pred_df_max <- pred_df %>%
#   group_by(colonizer,settmort,prienr1900) %>%
#   slice_max(probability, with_ties = FALSE) %>%
#   ungroup()
# 
# # Plot
# ggplot(pred_df_max, aes(x = settmort, y = prienr1900, fill = cluster)) +
#   geom_tile() +
#   geom_point(data=df_model, aes(x=settmort, y=prienr1900, fill=factor(cluster)), color="black",pch=21, alpha=0.5,size=3)+
#   facet_wrap(~colonizer) +
#   theme_minimal() +
#   labs(title = "Predicted Probabilities by Colonizer, Settler mortality and school enrollment (1900)")
# 
# 
# 
# 
# 
# 
# #Plot as function of settmort (fit 2 due to better parameter estimates)
# preds <- data.frame(ggeffects::ggemmeans(fit2, terms=~settmort+colonizer))
# ggplot(preds, aes(x=x, y=predicted, color=response.level))+
#   geom_line(show.legend = FALSE) +
#   geom_point(data=df_model%>%
#                select(response.level=cluster,
#                       group=colonizer,
#                       x=settmort)%>%
#                na.omit(),
#              aes(x=x,y=0, color=response.level),
#              shape=4,size=3, stroke = 1.5,alpha=0.7,
#              position = position_jitter(width = 0.1, height = 0))+
#   facet_wrap(~group)+
#   geom_ribbon(
#     aes(ymin = conf.low, ymax = conf.high, fill = response.level),
#     alpha = 0.2,
#     color = NA
#   )+
#   labs(x="log(Settler mortality)",
#        y="Probability",
#        fill="Cluster",
#        color="Observed clusters")+
#   theme_minimal()








# #Aggregate model------------------------------------------------------
# #Aggregate clusters not well represented in colonial data
# df_agg <- df_model %>%
#   mutate(cluster = if_else(cluster %in%c("2","3","5"), "2,3,5", as.character(cluster)))%>%
#   mutate(cluster= factor(cluster),
#          colonizer= factor(colonizer))%>%
#   filter(!is.na(colonizer))
# 
# table(df_agg$cluster, df_agg$colonizer)
# 
# fit_agg <- multinom(cluster~colonizer, data=df_agg)
# summary(fit_agg) #Parameters do not look good - too many clusters
# 
# z <- summary(fit_agg)$coefficients/ summary(fit_agg)$standard.errors
# p <- 2*(1-pnorm(abs(z),0,1))
# print(p)
# car::Anova(fit_agg,type=2)  
# 
# #How about settler mortality 
# fit_agg_settmort <- multinom(cluster~settmort, data=df_agg)
# summary(fit_agg_settmort)
# z <- summary(fit_agg_settmort)$coefficients/ summary(fit_agg_settmort)$standard.errors
# p <- 2*(1-pnorm(abs(z),0,1))
# print(p)
# car::Anova(fit_agg_settmort,type=2) #Settler mortality is a highly signifcicant predictor
# 
# #What if we include primary school enrollment
# car::Anova(update(fit_agg_settmort,~.+settmort*prienr1900),type=2) #Both settler mortality and primary school enrollment seems significant
# fit_agg <- update(fit_agg_settmort,~.+prienr1900,maxit=500)
# z <-  summary(fit_agg)$coefficients/ summary(fit_agg)$standard.errors
# p <- 2*(1-pnorm(abs(z),0,1))
# print(p)
# car::Anova(fit_agg,type=2) #Both are significant


# #Continuous graph
# preds <- data.frame(ggeffects::ggemmeans(fit_agg, terms=c("prienr1900")))
# ggplot(preds, aes(x=x, y=predicted,color = response.level))+
#   geom_line()+
#   geom_ribbon(
#     aes(ymin = conf.low, ymax = conf.high, fill = response.level),
#     alpha = 0.2,
#     color = NA
#   )+
#   labs(x="Primary school enrollment",
#        y="Probability",
#        fill="Colonizer")+
#   theme_minimal()
# 
# 
# preds <- data.frame(ggeffects::ggemmeans(fit_agg, terms=c("settmort")))
# ggplot(preds, aes(x=x, y=predicted,color = response.level))+
#   geom_line()+
#   geom_ribbon(
#     aes(ymin = conf.low, ymax = conf.high, fill = response.level),
#     alpha = 0.2,
#     color = NA
#   )+
#   labs(x="log(Settler mortality)",
#        y="Probability",
#        fill="Cluster")+
#   theme_minimal()+
#   guides(color="none")

# 
# #Barchart for colonizer
# preds <- data.frame(ggeffects::ggemmeans(fit_agg, terms="colonizer"))
# ggplot(preds, aes(x=x, y=predicted,fill = response.level))+
#   geom_col(position = position_dodge(width = 0.9)) +
#   geom_errorbar(
#     aes(ymin = conf.low, ymax = conf.high),
#     color = "black",
#     width = 0.2,
#     position = position_dodge(width = 0.9)
#   )+
#   labs(x="Colonizer",
#        y="Probability",
#        fill="Colonizer")+
#   theme_minimal()
# 
# #Including settler mortality
# df_agg2 <- df_agg %>% 
#   filter(!is.na(settmort))
# 
# fit_agg2 <- multinom(cluster~colonizer+settmort, data=df_agg2)
# summary(fit_agg2)
# 
# z <- summary(fit_agg2)$coefficients/ summary(fit_agg2)$standard.errors
# p <- 2*(1-pnorm(abs(z),0,1))
# print(p)
# car::Anova(fit_agg2,type=2)  
# lmtest::lrtest(update(fit_agg2,~1),fit_agg2)
# 
# 
# preds <- data.frame(ggeffects::ggemmeans(fit_agg2, terms=~settmort+colonizer))
# ggplot(preds, aes(x=x, y=predicted, color=response.level))+
#   geom_line(show.legend = FALSE) +
#   geom_point(data=df_agg2%>%
#                select(response.level=cluster,
#                       group=colonizer,
#                       x=settmort)%>%
#                na.omit(),
#              aes(x=x,y=0, color=response.level),
#              shape=4,size=3, stroke = 1.5,alpha=0.7,
#              position = position_jitter(width = 0.1, height = 0))+
#   facet_wrap(~group)+
#   geom_ribbon(
#     aes(ymin = conf.low, ymax = conf.high, fill = response.level),
#     alpha = 0.2,
#     color = NA
#   )+
#   labs(x="log(Settler mortality)",
#        y="Probability",
#        fill="Cluster",
#        color="Observed clusters")+
#   theme_minimal()
# 
# 
# #Include primary school enrollment?
# fit_agg3 <- update(fit_agg2,~.+prienr1900,maxit=200)
# summary(fit_agg3)
# z <- summary(fit_agg3)$coefficients/ summary(fit_agg3)$standard.errors
# p <- 2*(1-pnorm(abs(z),0,1))
# print(p)
# car::Anova(fit_agg3,type=2) 
# 
# preds <- data.frame(ggeffects::ggemmeans(fit_agg3, terms=~settmort+prienr1900+colonizer))
# lmtest::lrtest(fit_agg3,update(fit_agg2,data=df_agg2%>%filter(!is.na(prienr1900))))
# 
# #Limt the analysis to clusters represented in africa 
# df_filtered <- df_model %>%
#   filter(cluster %in% c(1,2,6,7,8))%>%
#   mutate(cluster= factor(cluster),
#          colonizer= factor(colonizer))%>%
#   filter(!is.na(colonizer))
#   # filter(!is.na(settmort))
# 
# fit_filtered <- multinom(cluster ~ colonizer,
#                 data=df_filtered,
#                 maxit=700)
# 
# summary(fit_filtered)
# z <- summary(fit_filtered)$coefficients/ summary(fit_filtered)$standard.errors
# p <- 2*(1-pnorm(abs(z),0,1))
# print(p)
# car::Anova(fit_filtered,type=2)  
# 
# fit_null <- multinom(cluster~1,
#                      data=df_filtered)
# lmtest::lrtest(fit_null,fit_filtered)
# 
# preds <- data.frame(ggeffects::ggemmeans(fit_filtered, terms="colonizer"))
# 
# ggplot(preds, aes(x=x, y=predicted,fill = response.level))+
#   geom_col(position = position_dodge(width = 0.9)) +
#   geom_errorbar(
#     aes(ymin = conf.low, ymax = conf.high),
#     color = "black",
#     width = 0.2,
#     position = position_dodge(width = 0.9)
#   )+
#   labs(x="Cluster",
#        y="Probability",
#        fill="Colonizer")+
#   theme_minimal()
# 
# #Add settler mortality as a variable
# df_filtered2 <- df_model %>%
#   filter(cluster %in% c(1,2,6,7))%>%
#   mutate(cluster= factor(cluster),
#          colonizer= factor(colonizer))%>%
#   filter(!is.na(colonizer))%>%
#   filter(!is.na(settmort))
# 
# fit_filtered2 <- multinom(cluster ~ colonizer+settmort,
#                          data=df_filtered2,
#                          maxit=700)
# 
# summary(fit_filtered2)
# z <- summary(fit_filtered2)$coefficients/ summary(fit_filtered2)$standard.errors
# p <- 2*(1-pnorm(abs(z),0,1))
# print(p)
# car::Anova(fit_filtered2,type=2)  
# 
# fit_null <- multinom(cluster~1,
#                      data=df_filtered2)
# lmtest::lrtest(fit_null,fit_filtered)
# 
# 
# 
# #Add primary school enrollment as a variable
# df_filtered2 <- df_model %>%
#   filter(cluster %in% c(1,2,6,7))%>%
#   mutate(cluster= factor(cluster),
#          colonizer= factor(colonizer))%>%
#   filter(!is.na(colonizer))%>%
#   filter(!is.na(settmort))
# 
# fit_filtered2 <- multinom(cluster ~ colonizer+settmort,
#                           data=df_filtered2,
#                           maxit=700)
# 
# summary(fit_filtered2)
# z <- summary(fit_filtered2)$coefficients/ summary(fit_filtered2)$standard.errors
# p <- 2*(1-pnorm(abs(z),0,1))
# print(p)
# car::Anova(fit_filtered2,type=2)  
# 
# fit_null <- multinom(cluster~1,
#                      data=df_filtered2)
# lmtest::lrtest(fit_null,fit_filtered)
# 
# 
# ##Visualize model
# # Make prediction grid
# newdata <- expand.grid(
#   colonizer = levels(df_filtered2$colonizer),
#   settmort = seq(min(df_filtered2$settmort,na.rm=TRUE), max(df_filtered$settmort,na.rm=TRUE), length.out = 50)
# )
# 
# # Get predicted probabilities
# preds <- predict(fit_filtered2, newdata = newdata, type = "probs")
# 
# # Convert to long format
# pred_df <- cbind(newdata, preds) %>%
#   tidyr::pivot_longer(cols = -c(colonizer,settmort),
#                       names_to = "cluster", values_to = "probability")%>%
#   unique()
# 
# # Plot
# ggplot(pred_df, aes(x = settmort, y = probability, color = cluster)) +
#   geom_line() +
#   geom_point(data=df_filtered2, aes(x=settmort, y=1, fill=factor(cluster)), color="black",pch=21, alpha=0.5,size=3)+
#   facet_wrap(~colonizer) +
#   theme_minimal() +
#   labs(title = "Predicted Probabilities by Colonizer, Settler mortality and school enrollment (1900)")


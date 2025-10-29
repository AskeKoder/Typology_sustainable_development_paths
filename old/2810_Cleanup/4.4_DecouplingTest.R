library(dplyr)

#Read data
data <- read.csv("TFPdata.csv")%>%
  select(-X)

clusters <-readRDS("4_RankedClusters.RDS")%>%
  select(Country,iso3=SPI_countrycode,Cluster = DLSFew_coverage)

#Combine
df <- merge(data,clusters, by=c("iso3","Country"))%>%
  mutate(Cluster = factor(Cluster),
         iso3 = factor(iso3))

table(df$Cluster)/21


#Compute EF-------------------------------------
#World average citizen in 2020
world_citizen <- df%>%
  filter(Year == 2020)%>%
  summarise(across(c(GHG, Biodiversity_Impact, Scarce_Water_Consumption), ~ weighted.mean(., w = Population)))

#Normalize by world average citizen
df<- df %>%
  mutate(
    nGHG = GHG/world_citizen$GHG,
    nBiodiversity_Impact = Biodiversity_Impact / world_citizen$Biodiversity_Impact,
    nScarce_Water_Consumption = Scarce_Water_Consumption / world_citizen$Scarce_Water_Consumption,
    pers.eq.max = pmax(nGHG, nBiodiversity_Impact, nScarce_Water_Consumption),
    pers.eq.min = pmin(nGHG, nBiodiversity_Impact, nScarce_Water_Consumption),
    pers.eq.avg = rowMeans(select(.,c("nGHG", "nBiodiversity_Impact", "nScarce_Water_Consumption")))
  )

df$iso3 <- factor(df$iso3, levels = unique(df$iso3[order(df$Cluster, -df$pers.eq.min)]))

ggplot(df%>%filter(Year==2020), aes(x=iso3, y=pers.eq.min, fill=factor(Cluster)))+
  geom_col()


#Fit models ------------------------------------
fitALL <- lm(cbind(log(GHG),
                   log(Biodiversity_Impact),
                   log(Scarce_Water_Consumption)) ~ Cluster:(Year+TFP+log(GDP_PPPcap)+log(Population)),
             data=df)
summary(fitALL)
anova(fitALL)


fitGHG <- lm(log(GHG) ~ Cluster*(Year+TFP+log(GDP_PPPcap)+log(Population)):Country,
             data=df)
summary(fitGHG)
anova(fitGHG)

library(lme4)
lme.GHG <- lmer(log(GHG) ~ Cluster*(Year+TFP+log(GDP_PPPcap)+log(Population))+(1|iso3),
                data=df,REML=FALSE)
lme.GHG.rs <- lmer(log(GHG) ~ Cluster*(Year+TFP+log(GDP_PPPcap)+log(Population))+(log(GDP_PPPcap)|iso3),
                data=df,REML=FALSE)
car::Anova(lme.GHG, type="III")
car::Anova(lme.GHG.rs, type="III")
anova(lme.GHG,lme.GHG.rs) #We need to includ a random slope

#Model building lmer bottom up GHG ------------------------------------------------------------------------------
df
df_filtered <- df%>%
  filter(!is.na(GDP_PPPcap))%>%
  filter(!is.na(Population))%>%
  filter(!is.na(TFP))

m0 <- lmer(log(GHG) ~ 1 + (1|iso3),
             data=df_filtered,REML=FALSE)
m1 <- lmer(log(GHG) ~ 1 + (log(GDP_PPPcap)|iso3),
             data=df_filtered,REML=FALSE)
summary(m0)
VarCorr(m0)
anova(m0,m1) #Including the random slop is better
VarCorr(m1)
m2 <- lmer(log(GHG) ~ 1 + (scale(log(GDP_PPPcap))+scale(log(Population))|iso3),
           data=df_filtered,REML=FALSE)
anova(m2,m1) #Including more information to the random slope is better
VarCorr(m2)
m3 <- lmer(log(GHG) ~ 1 + (1+scale(Year)+scale(log(GDP_PPPcap))|iso3),
           data=df_filtered,REML=FALSE)
anova(m2,m3) #Including Year rather than population is best, do we need to add both?
VarCorr(m3)
m4 <- lmer(log(GHG) ~ 1 + (1+scale(Year)+scale(log(GDP_PPPcap))+scale(log(Population))|iso3),
           data=df_filtered,REML=FALSE)
anova(m3,m4) #Including both porivdes a moderate improvement to the model, now what about TFP as well
VarCorr(m4)

m5 <- lmer(log(GHG) ~ 1 + (1+scale(Year)+scale(log(GDP_PPPcap))+scale(log(Population))+scale(TFP)|iso3),
           data=df_filtered,REML=FALSE)
anova(m4,m5) #TFP is not needed in the random effects
VarCorr(m5)

#Now we  reset the filter and start including the fixed effects
df_filtered <- df%>%
  filter(!is.na(GDP_PPPcap))%>%
  filter(!is.na(Population))
m4 <- lmer(log(GHG) ~ 1 + (1+scale(Year)+scale(log(GDP_PPPcap))+scale(log(Population))|iso3),
           data=df_filtered,REML=FALSE)
summary(m4)

m4.1 <- lmer(log(GHG) ~ log(GDP_PPPcap) + (1+scale(Year)+scale(log(GDP_PPPcap))+scale(log(Population))|iso3) ,
           data=df_filtered,REML=FALSE)
summary(m4.1)
anova(m4,m4.1)

# m4.2 <- lmer(log(GHG) ~ log(GDP_PPPcap) +TFP + (1+scale(Year)+scale(log(GDP_PPPcap))+scale(log(Population))|iso3) ,
#              data=df_filtered,REML=FALSE) #non-convergence
# summary(m4.2)
# anova(m4.1,m4.2)

m4.3 <- lmer(log(GHG) ~ Cluster*log(GDP_PPPcap) + (1+scale(Year)+scale(log(GDP_PPPcap))+scale(log(Population))|iso3) ,
             data=df_filtered,REML=FALSE)
summary(m4.3)
anova(m4.1,m4.3)


m4.4 <- lmer(log(GHG) ~ Cluster*(log(GDP_PPPcap)) + (1+scale(log(GDP_PPPcap))+scale(log(Population))|iso3) ,
             data=df_filtered,REML=FALSE) #non-convergence
summary(m4.4)
anova(m4.3,m4.4) #Keeping year as a random effect is best



m4.5 <- lmer(log(GHG) ~ Cluster*(log(GDP_PPPcap)+TFP) + (1+scale(Year)+scale(log(GDP_PPPcap))+scale(log(Population))|iso3) ,
             data=df_filtered,REML=FALSE) 
summary(m4.3)
anova(m4.1,m4.3)


#
lmer(log(GHG) ~ 1+scale(log(GDP_PPPcap))+ (1+scale(Year)+log(GDP_PPPcap)+log(Population)|iso3),
     data=df_filtered,REML=FALSE)










lme.GHG.rs.fullmodel <- lmer(log(GHG) ~ Cluster*(Year+TFP+log(GDP_PPPcap)+log(Population))+(Year+TFP+log(GDP_PPPcap)+log(Population)|iso3),
                   data=df,REML=FALSE)
VarCorr(lme.GHG.rs.fullmodel)
summary(lme.GHG.rs2)
car::Anova(lme.GHG.rs2, type="III")


lme.GHG <- lmer(log(GHG) ~ Cluster*(Year+TFP+log(GDP_PPPcap)+log(Population))+(1|iso3),
                data=df,REML=TRUE)

summary(lme.GHG)
library(emmeans)

emmeansList <- list()

# Estimated marginal means per Group
#emmeansList[[1]] <- emmeans(lme.GHG, ~ Cluster)

# Simple slopes for A, B, C, F within each Group
model <- fit.HDI
emmeansList[[1]] <-emtrends(model, ~ Cluster, var = "Year")%>%data.frame()
emmeansList[[2]] <-emtrends(model, ~ Cluster, var = "TFP")%>%data.frame()
emmeansList[[3]] <-emtrends(model, ~ Cluster, var = "log(GDP_PPPcap)")%>%data.frame()
emmeansList[[4]] <-emtrends(model, ~ Cluster, var = "log(Population)")%>%data.frame()

p <- list()
for ( i in 2:3){
  p[[i]] <- ggplot(emmeansList[[i]], aes(x=Cluster, y = !!sym(colnames(emmeansList[[i]])[2])))+
    geom_point(aes(color=factor(Cluster)))+
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL))+
    theme_minimal()
}
library(patchwork)
(p[[1]]+p[[2]])/(p[[3]] + p[[4]]) +  plot_layout(guides = "collect")

p[[2]]+p[[3]]


#Following the approach from the energy article---------------------------------------------------

#Panel Unit root test-----------------------------------------------------------------------------------
library(plm)
df_mod <- df %>%
  mutate(logGHG = log(GHG),
         logGDPcap = log(GDP_PPPcap),
         logPOP = log(Population),
         logTFP = TFP)%>%
  group_by(iso3) %>%                  # ensure differences are within countries
  arrange(Year, .by_group = TRUE) %>%
  mutate(
    dlogGHG = logGHG - dplyr::lag(logGHG,1)                 # first difference
  )%>%
  select(-TFP)

test_frame <- data.frame(Cluster=1:11,basetest_pvalues=rep(NA,11),difftest_pvalues=rep(NA,11))
for (i in 1:11){
        testbase <- purtest(logGHG ~ 1, data = df_mod%>%
                          filter(Cluster==i),
        index = c("iso3","Year"), test = "ips", lags = 2)
        
        testdiff <- purtest(dlogGHG ~ 1, data = df_mod%>%
                          filter(Cluster==i),
                        index = c("iso3","Year"), test = "ips", lags = 2)
        
        test_frame[i,2] <- round(testbase$statistic$p.value,4)
        test_frame[i,3] <- round(testdiff$statistic$p.value,4)
}
#Warnings about time series being short




#Fit
fit.GHG <- lme4::lmer(log(GHG)~Cluster*(log(GDP_PPPcap)+TFP)+(scale(Year)|iso3),
                  data=df)
summary(fit)

fit.HDI <- lme4::lmer(HDI~Cluster*(log(GDP_PPPcap)+TFP)+(scale(Year)|iso3),
                      data=df)
summary(fit)

emmeans::emtrends(fit, )
#FMOLS

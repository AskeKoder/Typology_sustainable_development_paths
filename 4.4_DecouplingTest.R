#Read data
data <- read.csv("TFPdata.csv")%>%
  select(-X)

clusters <-readRDS("4_RankedClusters.RDS")%>%
  select(Country,iso3=SPI_countrycode,Cluster = DLSFew_coverage)

#Combine
df <- merge(data,clusters, by=c("iso3","Country"))%>%
  mutate(Cluster = factor(Cluster))


#Fit models ------------------------------------
fitALL <- lm(cbind(log(GHG),
                   log(Biodiversity_Impact),
                   log(Scarce_Water_Consumption)) ~ Cluster:(Year+TFP+log(GDP_PPPcap)+log(Population)),
             data=df)
summary(fitALL)
anova(fitALL)


fitGHG <- lm(log(GHG) ~ Cluster:(Year+TFP+log(GDP_PPPcap)+log(Population)),
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
anova(lme.GHG,lme.GHG.rs)

lme.GHG <- lmer(log(GHG) ~ Cluster*(Year+TFP+log(GDP_PPPcap)+log(Population))+(1|iso3),
                data=df,REML=TRUE)

summary(lme.GHG)
library(emmeans)

emmeansList <- list()

# Estimated marginal means per Group
#emmeansList[[1]] <- emmeans(lme.GHG, ~ Cluster)

# Simple slopes for A, B, C, F within each Group
emmeansList[[1]] <-emtrends(lme.GHG.rs, ~ Cluster, var = "Year")%>%data.frame()
emmeansList[[2]] <-emtrends(lme.GHG.rs, ~ Cluster, var = "TFP")%>%data.frame()
emmeansList[[3]] <-emtrends(lme.GHG.rs, ~ Cluster, var = "log(GDP_PPPcap)")%>%data.frame()
emmeansList[[4]] <-emtrends(lme.GHG.rs, ~ Cluster, var = "log(Population)")%>%data.frame()


for ( i in 1:4){
  p[[i]] <- ggplot(emmeansList[[i]], aes(x=Cluster, y = !!sym(colnames(emmeansList[[i]])[2])), color=factor(Cluster))+
    geom_point(aes(color=factor(Cluster)))+
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL,color=factor(Cluster)))+
    theme_minimal()
}
library(patchwork)
(p[[1]]+p[[2]])/(p[[3]] + p[[4]]) +  plot_layout(guides = "collect")

library(dplyr)
library(ggplot2)
#Comparison of imputations
Laglead <- read.csv("ImputedDataLag1Lead2_maxit30_scaled.csv")%>%
  select(-c(X,.id))%>%
  relocate(.imp, .after=last_col())%>%
  mutate(Country = factor(Country))

Lighthouse <- read.csv("ImputedDataLightHouse_scaled.csv")%>%
  select(-c(X,.id))%>%
  relocate(.imp, .after=last_col())%>%
  mutate(Country = factor(Country))

raw_data <- read.csv("extendedDataScaled.csv")%>%
  select(-X)%>%
  mutate(Country = factor(Country))

#Read data from subnational survey
load("Data/Complete DLS data file_country level.RData")
dls_country <- dls_country%>%
  rename(Country=country_name)%>%
  mutate(year_start_interviews = as.numeric(year_start_interviews))%>%
  rename(SPI_year = year_start_interviews)%>%
  filter(SPI_year >= 2000)


#Plot comparison for imputed data
var <- colnames(raw_data)[58]
var
ggplot() +
  geom_line(data = Lighthouse[Lighthouse$Country %in% unique(raw_data$Country[which(is.na(raw_data[,var]))]),],
  aes(x = SPI_year, y = !!sym(var), group = .imp),
  alpha = 0.2, color="blue") +
  geom_line(data = Laglead[Laglead$Country %in% unique(raw_data$Country[which(is.na(raw_data[,var]))]),],
            aes(x = SPI_year, y = !!sym(var), group = .imp),
            alpha = 0.2, color="red") +
  geom_point(data = dls_country[dls_country$Country %in% unique(raw_data$Country[which(is.na(raw_data[,var]))]),],
            aes(x = year_start_interviews, y = dim6_sanitation_country*100),
            alpha = 1, color="green")+ 
  geom_line(data=raw_data[raw_data$Country %in% unique(raw_data$Country[which(is.na(raw_data[,var]))]),]%>%
              mutate(.imp=0),
            aes(x = SPI_year, y = !!sym(var), group = .imp),
            alpha = 1, lwd=0.9,color="black") +
  facet_wrap(~ Country) +
  labs(color = "Imputation") +
  theme_minimal() +
  ylim(
    min(Lighthouse[[var]], na.rm = TRUE),
    max(Lighthouse[[var]], na.rm = TRUE)
  )

#Plot comparison for variables in DLS
var <- colnames(raw_data)[58]
var
ggplot() +
  # geom_line(data = Lighthouse[Lighthouse$Country %in% unique(dls_country$Country),],
  #           aes(x = SPI_year, y = !!sym(var), group = .imp),
  #           alpha = 0.2, color="blue") +
  geom_line(data = Laglead[Laglead$Country %in% unique(dls_country$Country),],
            aes(x = SPI_year, y = !!sym(var), group = .imp),
            alpha = 0.2, color="red") +
  geom_point(data = dls_country[dls_country$Country %in% unique(dls_country$Country),],
             aes(x = SPI_year, y = dim6_sanitation_country*100),
             alpha = 1, color="green")+ 
  geom_line(data=raw_data[raw_data$Country %in% unique(dls_country$Country),]%>%
              mutate(.imp=0),
            aes(x = SPI_year, y = !!sym(var), group = .imp),
            alpha = 1, lwd=0.9,color="black") +
  facet_wrap(~ Country) +
  labs(color = "Imputation") +
  theme_minimal() +
  ylim(
    min(Lighthouse[[var]], na.rm = TRUE),
    max(Lighthouse[[var]], na.rm = TRUE)
  )


#-------------------------------------------------------------------
comparisonData_lh <- dls_country%>%
  mutate(across(c(5:ncol(dls_country)), ~ .x * 100))%>%
  left_join(rbind(cbind(raw_data,".imp"=0),Lighthouse), by=c("Country","SPI_year"))
            
comparisonData_ll <- dls_country%>%
  mutate(across(c(5:ncol(dls_country)), ~ .x * 100))%>%
  left_join(rbind(cbind(raw_data,".imp"=0),Laglead), by=c("Country","SPI_year"))


#Variables to compare
variables <- data.frame("Dimension"=c("dim5_water_country",
                                      "dim6_sanitation_country",
                                      "dim8_education_country"),
                        "Indicator"=c("Safely_Managed_Drinking_Water",
                                      "Safely_managed_saniation",
                                      "Prim_School_Enroll"))


#Compute error measures
RMSE<- data.frame("water_lh"=rep(0,16),
                  "sani_lh" = rep(0,16),
                  "edu_lh" = rep(0,16), 
                  "water_ll" = rep(0,16),
                  "sani_ll" = rep(0,16),
                  "edu_ll" = rep(0,16))
for (j in 1:3){
  for (i in 0:15){
    #Calculate error for lh
    RMSE[i+1,j]  <- sqrt(mean(as.matrix((comparisonData_lh[comparisonData_lh$.imp==i,variables[j,"Dimension"]] -
              comparisonData_lh[comparisonData_lh$.imp==i,variables[j,"Indicator"]])^2), na.rm=TRUE)
    )
    #Calculate error for ll
    RMSE[i+1,j+3] <- sqrt(mean(as.matrix((comparisonData_ll[comparisonData_ll$.imp==i,variables[j,"Dimension"]] -
                                    comparisonData_ll[comparisonData_ll$.imp==i,variables[j,"Indicator"]])^2), na.rm=TRUE)
    )
  }
}

par(mfrow=c(3,1))
plot(RMSE$water_lh,type="l", col="blue",ylim=c(0,60))
lines(RMSE$water_ll,type="l", col="red")
plot(RMSE$sani_lh,type="l", col="blue",ylim=c(0,60))
lines(RMSE$sani_ll,type="l", col="red")
plot(RMSE$edu_lh,type="l", col="blue",ylim=c(0,60))
lines(RMSE$edu_ll,type="l", col="red")






#Furthermore we plot the imputations
naIndicators <- names(which(colSums(is.na(raw_data)) > 0))
plot_list <- list()
for (i in 1:length(naIndicators)){
  indicator <- naIndicators[i]
  #Find countries with missing indicators
  countries = unique(raw_data$Country[is.na(raw_data[,indicator])])
  
  p <- ggplot()+
    geom_line(data=Laglead%>%filter(Country%in%countries), aes(x=SPI_year, y=!!sym(indicator), group=.imp),
              color="red",alpha=0.2)+
    geom_line(data=Lighthouse%>%filter(Country%in%countries), aes(x=SPI_year, y=!!sym(indicator), group=.imp),
              color="blue",alpha=0.2)+
    geom_point(data=raw_data%>%filter(Country%in%countries),aes(x=SPI_year, y=!!sym(indicator)),
                                                              color="black")+
    facet_wrap(~Country)+
    ylim(0,100)
  print(p)
  plot_list[[i]] <- p
}

# Save plots to png. Makes a separate file for each plot.
for (i in 1:length(naIndicators)) {
  indicator <- naIndicators[i]
  file_name = paste("Figures/ImputationComparison/", indicator, ".png", sep="")
  png(file_name,width = 1100, height = 800)
  print(plot_list[[i]])
  dev.off()
}

rm(list = ls())

#Import libraries 
library(mice) #For imputations
library(dplyr) #For data management
library(tidyr) #For conversion between wide/long
library(lattice) #For matrixplots/levelplots
library(ggplot2) #for plotting
library(zoo) #For linear interpolation

#Load data
data <- read.csv("extendedDataScaled.csv")%>%
  select(-X)

#Descriptive analysis -----------------------------------
nYears <- length(unique(data$SPI_year))
nCountries <- length(unique(data$Country))
nIndicators <- length(colnames(data)) - 4 #-identifiers

#Overview of missing data per indicator and per country
naIndicators <- data%>%
  is.na()%>%
  colSums()/(nYears*nCountries)
hist(naIndicators,breaks=20, xlab="Fraction of missing values", ylab="Number of indicators")

naOverview <- data%>%
  group_by(Country)%>%
  summarise_all(~sum(is.na(.)))
levelplot(as.matrix(naOverview))
#6 countries are missing the entire set of extension variables, all countries are well observed in SPI

#Correlation matrix 
levelplot(cor(as.matrix(data[,5:67]), use="complete.obs"),xlab="",ylab="",
          scales = list(
            x = list(
              #at = c(-3, -1.5, 0, 1.5, 3),   # tick positions
              labels = 4:66  # tick labels
            )
          ),
          col.regions = colorRampPalette(c("blue","white","red"))(100))



#################################################################################
#Iputation of Share slums by linear interpolation 
##############################################################################
#Imputation of Share slums is linear interpolation if surrounding years are observed
df_interpolated <- data %>%
  group_by(Country) %>%
  arrange(SPI_year) %>%
  mutate(Share_Slums_interp = na.approx(Share_Slums, x = SPI_year, na.rm = FALSE,
                                        #maxgap=1
                                        )) %>%
  ungroup()

#Visually inspect interpolations
C <- unique(data$Country)[140]
ggplot()+
  geom_point(data=df_interpolated,aes(x=SPI_year,y=Share_Slums_interp,color="Interpolated"))+
  geom_point(data=data,aes(x=SPI_year,y=Share_Slums,color="Observed"))+
  labs(title=paste("Interpolated values for"))+
  facet_wrap("Country")

#Incorporate changes
data <- data %>%
  group_by(Country) %>%
  arrange(SPI_year) %>%
  mutate(Share_Slums = na.approx(Share_Slums,
                                 x = SPI_year, 
                                 na.rm = FALSE,
                                 #maxgap=1
                                 )) %>%
  ungroup()

################################################################################
#Setup predictor matrix for mice
###############################################################################
#Convert data to wide format to include time dependency
indicator_list <- colnames(data)[5:ncol(data)]
data_wide <- data %>%
  pivot_wider(names_from = SPI_year, values_from = all_of(indicator_list), names_sep = "yearID")


#1 ) Initialize predictor matrix excluding categorical values
pred <- quickpred(data_wide,
                  mincor=0.6,
                  minpuc =0.9,
                  method = "pearson",
                  exclude = c("Country_0","SPI_countrycode_1")
)
table(rowSums(pred)) #100-1000 parameters for each model is far too much




#2 ) Identify variables that should always act as predictors
colors <-rgb(colSums(is.na(data[,c(2,5:length(data))]))>0,
             colSums(is.na(data[,c(2,5:length(data))]))>0,
             colSums(is.na(data[,c(2,5:length(data))]))>0)

par(mfrow=c(1,1))
corr <- cor(data[,5:ncol(data)],use="complete.obs")
labels <- 4:(ncol(data)-1)
qgraph::qgraph(corr, layout="spring",threshold= 0,
       labels=labels,
       vsize=3.5,repulsion=0.75,
       color = colors)
keep <- colnames(data)[c(5,7,12,22,31,33,40,50)]



#3 ) Setup lag(1) lead (2) moving time window (MTW)
for (var in rownames(pred)) {
  
  #Disregard variables without missing values
  if (sum(is.na(data_wide[,var])) == 0) {next}
  #Disregard non numeric columns
  if (var == "Country_0" || var == "SPI_countrycode_1" || var =="Region_3") {next}
  
  #Get current variable-year
  var_name <- substr(var,1,nchar(var)-4)
  year <- substr(var,nchar(var)-3,nchar(var))
  
  #Exclude all predictors from other years
  othr_yrs <- grep(paste0("*",year),rownames(pred), invert=TRUE)
  pred[var,othr_yrs] <-0
  
  #Set up MTW
  var_othr_yrs <- paste0(var_name, #indicator
                         c(max(c(2000,as.integer(year)-1)), #year before (or 200)
                           min(c(as.integer(year)+2,2020)))) # two years after (or 2020)
  
  #include in predictor matrix
  pred[var,var_othr_yrs] <- 1 #Let lag and lead predict the target
  pred[var,var] <- 0 #The target can't predict itself
  
  #Include "keep" variables for the current year
  pred[var,paste0(keep,"yearID",as.integer(year))] <- 1
  
}
#Avoid using the categorical variables as predictors since we have too many categories
pred[, "Country"] <- 0
pred[, "SPI_countrycode"] <- 0
table(rowSums(pred))


#4 ) Break direct feedback loops

#Go through every coordinate in the predictor matrix
for (i in 1:ncol(pred)){
  for (j in (i+1):ncol(pred)){
    if (j>ncol(pred)){break} #ends loop at i=ncol(pred)
    
    #If variables are used to predict each other
    if (pred[i,j] == 1 && pred[j,i] == 1){
      
      #Get names of both variables
      var_i <- rownames(pred)[i]
      var_j <- rownames(pred)[j]
      
      #If more observations are missing in i than in j
      if (sum(is.na(data_wide[,var_i])) >= sum(is.na(data_wide[,var_j]))){
        
        #Do not use i as a predictor for j
        pred[j,i] <- 0
        print(paste(var_i, "no longer predicts:", var_j))
      }
      else{ #if j has more missing values
        #do not use j as a predictor for i
        pred[i,j] <- 0
        print(paste(var_j, "no longer predicts:", var_i))
      }
    }
    else {next} #Skip if variables are not used to predict each other
  }
}
#Check
table(rowSums(pred)) #Variable years are being predicted with 8-35 predictors

###############################################################################
#Setup mice
###############################################################################
#Initialize
ini <- mice(data_wide,pred=pred,maxit=0)
ini$loggedEvents
unique(ini$loggedEvents[,"meth"]) 
ini$loggedEvents[ini$loggedEvents[,"meth"]=="constant",] #Primary school 2020 has no observations, must be imputed in long format

#Post processing not necessary as we use pmm
post <- ini$post

#Every indicator should be imputed using pmm, unless number number of donors is too small
meth <- ini$method
meth[] <- "pmm"

#Very few observations for PSE 2019 is bad for pmm, therefore this year is also predicted in long format
meth ["Prim_School_EnrollyearID2019"] <- ""
pred [ ,"Prim_School_EnrollyearID2019"] <- 0

#Set seed and settings
seed <- 123
m <- 15
maxit <- 30

#Run
imp<- mice(data_wide,
           defaultMethod = c("pmm", "logreg", "polyreg", "polr"),
           pred = pred,
           seed = seed,
           maxit = maxit,
           m = m,
           method = meth,
           post = post,
           print = TRUE, #running takes some time, so it's nice to follow the process
           remove.collinear = FALSE #FALSE = columns are not automatically deleted due to high correlation, which happens in adjacent years
           )  

#Check logged events
unique(imp$loggedEvents[,"meth"])
imp$loggedEvents[imp$loggedEvents[,"meth"]=="pmm",]  #Some targets have their predictors removed due to collinearity
nrow(imp$loggedEvents)
#Get one iteration of one imputation
check <- imp$loggedEvents[imp$loggedEvents[,"im"]==1,]
nrow(check) #377
overview <- imp$loggedEvents




#################################################################################
#Assess the imputations
#################################################################################

#Convergence plots 
NAlist <-  colnames(data_wide)[colSums(is.na(data_wide))>0]
i<-1

#Manual loop through plots
#Lines should ideally be entangled and without trend towards the end
#However, not meeting those conditions for all variable-years
#is not devastaing as inferential 
#validity does not depend on formal convergence (Buuren 2020)
plot(imp, NAlist[i:(i+17)],layout=c(6,6))
i <- i + 18



#Gather results
complete_wide <- complete(imp,"long")

#Convert back to long format
complete_long <- complete_wide %>%
  pivot_longer(
    cols = -c("Country","SPI_countrycode","Region",".imp",".id"),  # Keep these as identifier columns
    names_to = c(".value", "SPI_year"),  # Split names into variable and year
    names_sep = "yearID"  # Separator used in pivot_wider
  ) %>%
  mutate(SPI_year = as.integer(SPI_year))

complete_long
nrow(complete_long)/m
ncol(complete_long) #Two more columns are added in imputation
colSums(is.na(complete_long)) #NAs are only present in PSE

#Visualize imputations
var <- colnames(data)[58]
var  # variable of interest
ggplot(
  complete_long[complete_long$Country %in% unique(data$Country[which(is.na(data[,var]))]),],
  aes(x = SPI_year, y = !!sym(var), group = .imp)
) +
  #geom_point(alpha = 0.3) +
  geom_line(alpha = 0.15, color="blue") +
  geom_line(data=data[data$Country %in% unique(data$Country[which(is.na(data[,var]))]),]%>%
              mutate(.imp=0), alpha = 1, lwd=0.9,color="black") +
  facet_wrap(~ Country) +
  labs(color = "Imputation") +
  theme_minimal() +
  ylim(
    min(complete_long[[var]], na.rm = TRUE),
    max(complete_long[[var]], na.rm = TRUE)
  )


###################################################################
#Impute remaining values in long format
###################################################################
complete2 <- complete_long
#For each imputed set impute the remaining missing values in PSE
for (i in 1:m){
  #Filter data and create lag column for PSE
  df <- complete2[complete2$.imp==i,]%>%
    group_by(Country) %>%        # replace with your grouping variable(s)
    arrange(SPI_year, .by_group = TRUE) %>%  # optional: ensure proper order
    mutate(PSE_lag = lag(Prim_School_Enroll,2)) %>%
    ungroup()
  
  #Set up predictor matrix
  pred_long <- quickpred(df)
  pred_long[,"Country"] <- 0
  pred_long["Prim_School_Enroll","PSE_lag"] <- 1
  
  #Initialize
  ini_long <- mice(df,maxit=0)
  meth_long <- ini_long$method
  meth_long["PSE_lag"] <- "" #Don't impute lag
  
  #Run
  imp_long <- mice(df,
                   pred = pred_long,
                   method = meth_long,
                   post = ini_long$post,
                   #seed=seed,
                   m=1,
                   maxit=7)
  new_data <- complete(imp_long, "long")
  new_data[,".imp"] <- i
  new_data <- new_data %>%
    select(-PSE_lag)
  complete2[complete2$.imp==i,] <- new_data
}
sum(is.na(complete2))

#Save imputed data
#write.csv(complete2, "2_ImputedData.csv")

#Get data
filename <- file.choose()
complete2 <- read.csv(filename)%>%
  select(-X)


#Check if imputations seem viable
missingIndicators <- colnames(data)[which(colSums(is.na(data))>0)]
var <- missingIndicators[5]
var  # variable of interest


#Visualize imputations, #Lags could be included 
ggplot(
  complete2[complete2$Country %in% unique(data$Country[which(is.na(data[,var]))]),],
  aes(x = SPI_year, y = !!sym(var), group = .imp)
) +
  #geom_point(alpha = 0.3) +
  geom_line(alpha = 0.15, color="blue") +
  geom_line(data=data[data$Country %in% unique(data$Country[which(is.na(data[,var]))]),]%>%
              mutate(.imp=0), alpha = 1, lwd=0.9,color="black") +
  facet_wrap(~ Country) +
  labs(color = "Imputation") +
  theme_minimal() +
  ylim(
    min(complete2[[var]], na.rm = TRUE),
    max(complete2[[var]], na.rm = TRUE)
  )



#Visualizing safely managed drinking water with improved water source
var <- missingIndicators[19]
var 
ggplot(
  complete2[complete2$Country %in% unique(data$Country[which(is.na(data[,var]))]),],
  aes(x = SPI_year, y = !!sym(var), group = .imp)
) +
  #geom_point(alpha = 0.3) +
  geom_line(alpha = 0.15, color="blue") +
  geom_line(data=data[data$Country %in% unique(data$Country[which(is.na(data[,var]))]),]%>%
              mutate(.imp=0), alpha = 1, lwd=0.9,color="black") +
  geom_line(data=data[data$Country %in% unique(data$Country[which(is.na(data[,var]))]),]%>%
              mutate(.imp=0),
            aes(x=SPI_year, y=Improved_water_source_proportion_of_pop*100),
            alpha = 1, lwd=0.9,color="red") +
  facet_wrap(~ Country) +
  labs(color = "Imputation") +
  theme_minimal() +
  ylim(
    min(complete2[[var]], na.rm = TRUE),
    max(complete2[[var]], na.rm = TRUE)
  )


#Visualizing safely managed sanitation with improved sanitatio
var <- missingIndicators[13]
var 
ggplot(
  complete2[complete2$Country %in% unique(data$Country[which(is.na(data[,var]))]),],
  aes(x = SPI_year, y = !!sym(var), group = .imp)
) +
  #geom_point(alpha = 0.3) +
  geom_line(alpha = 0.15, color="blue") +
  geom_line(data=data[data$Country %in% unique(data$Country[which(is.na(data[,var]))]),]%>%
              mutate(.imp=0), alpha = 1, lwd=0.9,color="black") +
  geom_line(data=data[data$Country %in% unique(data$Country[which(is.na(data[,var]))]),]%>%
              mutate(.imp=0),
            aes(x=SPI_year, y=Improved_sanitation_proportion_of_pop*100),
            alpha = 1, lwd=0.9,color="red") +
  facet_wrap(~ Country) +
  labs(color = "Imputation") +
  theme_minimal() +
  ylim(
    min(complete2[[var]], na.rm = TRUE),
    max(complete2[[var]], na.rm = TRUE)
  )

#plot maps of missingness 
sumMis <- data %>%
  group_by(SPI_countrycode) %>%
  summarise(across(everything(), ~ sum(is.na(.))))%>%
  select(SPI_countrycode,Prim_School_Enroll,Safely_managed_saniation,Safely_Managed_Drinking_Water)

library(rnaturalearth)
library(patchwork)
world <- ne_countries(scale = "medium", returnclass = "sf",continent = c("south america","oceania","north america", "asia","europe","africa"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'SDZ',"SDN"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'PN1',"PNG"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'PR1',"PRT"))%>%
  mutate(adm0_iso = replace(adm0_iso,  adm0_iso == 'SSD',"SSD*"))
colnames(world)[57] <- "SPI_countrycode"

world <- left_join(world, sumMis, 
                   by = "SPI_countrycode")
p <- list()
for (i in 2:ncol(sumMis)){
  p[[i-1]] <- ggplot() +
    geom_sf(data = world, aes(fill = !!sym(colnames(sumMis)[i])/21*100), color = "white",size=0.5)+
    labs(title = gsub("95","90",colnames(sumMis))[i])+ 
    theme_bw() + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())
}
p[[1]]/p[[2]]/p[[3]]

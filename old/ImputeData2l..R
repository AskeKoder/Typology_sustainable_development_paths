#Import libraries 
library(mice)
library(miceadds)
library(dplyr)
library(lattice)

#Load data
data <- read.csv("extendedData.csv")%>%
  select(-X)

#Descriptive analysis -----------------------------------
nYears <- length(unique(data$SPI_year))
nCountries <- length(unique(data$Country))
nIndicators <- length(colnames(data)) - 4 #-identifiers

#Overview of missing data per indicator and per countrye
naIndicators <- data%>%
  is.na()%>%
  colSums()/(nYears*nCountries)
hist(naIndicators,breaks=20)

naOverview <- data%>%
  group_by(Country)%>%
  summarise_all(~sum(is.na(.)))

levelplot(as.matrix(naOverview))

#Iputation of Share slums by linear interpolation ----------------------------------------
#Imputation of Share slums is linear interpolation if surrounding years are observed
library(zoo)
df_interpolated <- data %>%
  group_by(Country) %>%
  arrange(SPI_year) %>%
  mutate(Share_Slums_interp = na.approx(Share_Slums, x = SPI_year, na.rm = FALSE)) %>%
  ungroup()

#Visually inspect interpolations
C <- unique(data$Country)[9]
library(ggplot2)
ggplot()+
  geom_point(data=df_interpolated%>%filter(Country==C),aes(x=SPI_year,y=Share_Slums_interp,color="Interpolated"))+
  geom_point(data=data%>%filter(Country==C),aes(x=SPI_year,y=Share_Slums,color="Observed"))+
  labs(title=paste("Interpolated values for",C))

#Incorporate changes
data <- data %>%
group_by(Country) %>%
  arrange(SPI_year) %>%
  mutate(Share_Slums = na.approx(Share_Slums, x = SPI_year, na.rm = FALSE)) %>%
  ungroup()

#Setup mice for imputations---------------------------------------------------
# INPUT: your data frame in long format
# df_long with columns: Group, Year, v1, v2, v3, v4, v5, ...
# Adjust these variables if different
group_col <- "Country"
time_col  <- "SPI_year"

# list of time-varying vars to consider for imputation/predictors:
tv_vars <- colnames(data)[5:ncol(data)]  # extend as needed

# USER TUNABLES
cor_thresh <- 0.25   # min abs(cor) to consider a predictor useful
max_pred   <- 6      # max predictors to keep per target (besides lag/lead and means)

# ---------------------------------------------------------------------
# 1) Preprocess: make sure group is factor, center time, create lag/lead & cluster means and within vars
df <- data %>%
  mutate(
    !!group_col := as.factor(.data[[group_col]])
  ) %>%
  arrange(.data[[group_col]], .data[[time_col]])

# create lag1 and lead1 for each tv var, and cluster means + within variables
for(var in tv_vars){
  lagname  <- paste0(var, "_lag1")
  leadname <- paste0(var, "_lead1")
  meanname <- paste0(var, "_grpmean")
  withinname <- paste0(var, "_within")
  df <- df %>%
    group_by(.data[[group_col]]) %>%
    arrange(.data[[time_col]]) %>%
    mutate(
      !!lagname  := lag(.data[[var]]),
      !!leadname := lead(.data[[var]]),
      !!meanname := mean(.data[[var]], na.rm = TRUE)
    ) %>%
    ungroup()
  # compute within = value - group mean (NA if both NA)
  df[[withinname]] <- df[[var]] - df[[meanname]]
}

# center time to help random slope mixing
df[[paste0(time_col, "_c")]] <- df[[time_col]] - mean(df[[time_col]], na.rm = TRUE)
time_c <- paste0(time_col, "_c")

# ---------------------------------------------------------------------
# 2) Candidate predictors: contemporaneous (same row) tv_vars, lags, leads, group means, baselines
# Build numeric matrix of contemporaneous observed values for correlation screening
cont_vars <- tv_vars
cont_df <- df %>% select(all_of(cont_vars))

# compute pairwise correlations using pairwise complete obs
cors <- cor(cont_df, use = "pairwise.complete.obs")

# helper to choose predictors for a given target variable
choose_predictors <- function(target, cors_mat, tv_vars, cor_thresh, max_pred){
  # exclude self
  cvec <- cors_mat[target, ]
  cvec[target] <- 0
  s <- sort(abs(cvec), decreasing = TRUE)
  # select those above threshold
  picks <- names(s)[s > cor_thresh]
  if(length(picks) > max_pred) picks <- picks[1:max_pred]
  return(picks)
}

# compute chosen predictors per target
chosen_list <- list()
for(v in tv_vars){
  picks <- choose_predictors(v, cors, tv_vars, cor_thresh, max_pred)
  chosen_list[[v]] <- picks
}

# ---------------------------------------------------------------------
# 3) Build mice predictorMatrix with multilevel codes:
# - mark group_col as -2 in the column (group indicator)
# - for each target row (a variable to be imputed), set:
#    * include its lag1 and lead1 (if present) with code 1 (fixed)
#    * include the chosen contemporaneous predictors as fixed (1)
#    * include the group-mean of those predictors as code 3 (contextual)
#    * optionally include time_c as random slope (code 2)
# NOTE: mice expects predictorMatrix rows = variables, cols = variables

m <- make.predictorMatrix(df)
meth <- make.method(df)

# turn everything off to start
m[,] <- 0

# Set group identifier column as special (-2) for all imputation models
m[ , group_col] <- -2
# never impute group or time (original)
m[group_col, ] <- 0
m[time_col, ]  <- 0
m[time_c, ]    <- 0

# Decide which variables we'll impute: typically tv_vars + their lag/lead/within/means
# Ensure lag/lead/mean/within cols exist in df
all_cols <- names(df)

for(target in tv_vars){
  # set method for target to a two-level method
  meth[target] <- "2l.pmm"   # change to "2l.norm" or others if desired
  
  # include centered time as random slope
  m[target, time_c] <- 2
  
  # include lag/lead if present
  lagn  <- paste0(target, "_lag1")
  leadn <- paste0(target, "_lead1")
  if(lagn %in% all_cols) m[target, lagn] <- 1
  if(leadn %in% all_cols) m[target, leadn] <- 1
  
  # include chosen contemporaneous predictors
  preds <- chosen_list[[target]]
  for(p in preds){
    # include the raw predictor as fixed (1)
    if(p %in% all_cols) m[target, p] <- 1
    # include its group mean as contextual fixed effect (3)
    gm <- paste0(p, "_grpmean")
    if(gm %in% all_cols) m[target, gm] <- 3
  }
  
  # include the target's own group mean (context)
  gm_target <- paste0(target, "_grpmean")
  if(gm_target %in% all_cols) m[target, gm_target] <- 3
}

df$Country <- as.integer(df$Country)
imp <- mice(df, method = meth, predictorMatrix = m, m = 5, maxit = 15, printFlag = TRUE)






#Prepare imputation 
ini <- mice(data_wide,pred=pred,maxit=0)
ini$loggedEvents
unique(ini$loggedEvents[,"meth"]) 
ini$loggedEvents[ini$loggedEvents[,"meth"]=="constant",] #Primary school 2020 has no observations, must be imputed in long format

#Post processing not necessary as we use pmm
post <- ini$post

#Every indicator should be imputed using pmm
meth <- ini$method
meth[] <- "pmm"

#Set seed
seed <- 123

m <- 5
maxit <- 5
imp<- mice(data_wide,
           defaultMethod = c("pmm", "logreg", "polyreg", "polr"),
           pred = pred,
           seed = seed,
           maxit = maxit,
           m = m,
           method = meth,
           post = post,
           print = TRUE,
       remove.collinear = FALSE) #FALSE = columns are not automatically deleted due to high correlation 

unique(imp$loggedEvents[,"meth"])
imp$loggedEvents[imp$loggedEvents[,"meth"]=="pmm",]  #Some targets have their predictors removed due to collinearity
nrow(imp$loggedEvents)
#Get one iteration of one imputation
check <- imp$loggedEvents[imp$loggedEvents[,"im"]==1,]
nrow(check) #377

#Inspect quality of the imputations -------------------------------------------------


#Convergence plots
NAlist <-  colnames(data_wide)[colSums(is.na(data_wide))>0]
i<-1
##Loop does not work - plot -save command must be copy-pasted into console repeatedly
# while (i < length(NAlist)){
#   print(i)
plot(imp, NAlist[i:(i+17)],layout=c(6,6),xlim=c(15,25))
#   dev.copy(file=paste0("Ch",i,".png"),device=png, bg="white", width=1764, height = 784) 
#   graphics.off()
i <- i + 18

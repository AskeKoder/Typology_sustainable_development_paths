
library(readxl)
library(readr)
library(dplyr)
library(nnet)
library(ggplot2)

# Load data
agr <- read_excel("xcountry_data_full.xlsx")
names(agr)[names(agr) == 'code'] <- 'iso3'

#New clusters from experiments
clusters <- readRDS("clusterVariations.RDS")%>%
  rename(iso3 = SPI_countrycode)

# Clean ISO3 values
clusters <- clusters %>% mutate(iso3 = toupper(trimws(as.character(iso3))))
agr      <- agr      %>% mutate(iso3 = toupper(trimws(as.character(iso3))))

# Merge
df <- clusters %>% inner_join(agr, by = "iso3")

# Simplified colonizer grouping
df$colonizer <- "Other"
df$colonizer[df$f_brit == 1]   <- "British"
df$colonizer[df$f_french == 1] <- "French"
df$colonizer[df$f_spain == 1]  <- "Spanish"
df$colonizer <- factor(df$colonizer)

# Make clusters a factor
df <- df %>%
  mutate(across(3:15,as.factor))

colnames(df)
# Subset needed variables
df_model <- df %>% select(cluster = "Few indicators_SPI_preferred_90" , colonizer, lcapped, lpd1500s, lat_abst,
                          ruleoflaw, protmiss, prienr1900) %>%
  na.omit()


#Visualize data -----------------------------------------------------------
plot(df_model,col=as.numeric(df_model$cluster))
table(df_model$cluster, df_model$colonizer)


#Count occurences per cluster
df_counts <- df_model %>%
  group_by(cluster, colonizer) %>%
  summarise(n = n(), .groups = "drop")  # count rows per group

ggplot(df_counts, aes(x = cluster, y = n, fill = colonizer)) +
  geom_col(position = "dodge") +
  labs(
    title = "Counts of colonizers within each cluster",
    x = "Cluster",
    y = "Count",
    fill = "Colonizer"
  ) +
  theme_minimal()

#PCA -----------------------------------------------
pca_result <- prcomp(df_model%>%
                       select(is.numeric),scale=TRUE,center = TRUE)

#Scree plot                    
plot(pca_result$sdev^2/sum(pca_result$sdev^2))

#Plotting for the first three dimensions
pc_scores <- cbind(df_model%>%select(!is.numeric),pca_result$x[, 1:3])

library(patchwork)

p1 <- ggplot(data=pc_scores ,
             aes(x=PC1,y=PC2, colour = cluster,shape=colonizer))+
  geom_point()
p2 <- ggplot(data=pc_scores
             , aes(x=PC2,y=PC3, colour = cluster,shape=colonizer))+
  geom_point()
p3 <- ggplot(data=pc_scores
             , aes(x=PC1,y=PC3, colour = cluster,shape=colonizer))+
  geom_point()
p1/p3+p2 +plot_layout(guides = "collect")


#Fit model ------------------------------------------------------------------

#Base model
model_base <- multinom(cluster ~ 1, data = df_model,
                       trace=TRUE)


# Model without one variable at a time:
model_col  <- multinom(cluster ~ colonizer, data = df_model,trace=TRUE)
model_rule  <- multinom(cluster ~ ruleoflaw, data = df_model,trace=TRUE)
model_set  <- multinom(cluster ~ lcapped, data = df_model,trace=TRUE)
model_lat  <- multinom(cluster ~ lat_abst, data = df_model,trace=TRUE)
model_pop <- multinom(cluster ~ lpd1500s, data = df_model,trace=TRUE)
model_prot <- multinom(cluster ~ protmiss, data = df_model, trace=TRUE)
model_school <- multinom(cluster ~ prienr1900, data = df_model, trace=TRUE)


# Compare log-likelihoods
delta_logLik <- function(model_big) {
  2 * (logLik(model_big)-logLik(model_base))
}

var_contrib <- tibble::tibble(
  Variable = c("Colonizer","Rule of law (2005)","Settler mortality", "Latitude", "Pop. density (1500)",
               "Protestant missions", "School enrollment (1900)"),
  DeltaChi2 = c(
    delta_logLik(model_col),
    delta_logLik(model_rule),
    delta_logLik(model_set),
    delta_logLik(model_lat),
    delta_logLik(model_pop),
    delta_logLik(model_prot),
    delta_logLik(model_school)
  )
  
)

# Plot
ggplot(var_contrib, aes(x = reorder(Variable, DeltaChi2), y = DeltaChi2)) +
  geom_col(fill = "#3366CC") +
  coord_flip() +
  labs(title = "Variable Contributions to Cluster Prediction",
       subtitle = "Based on Likelihood Ratio Tests (Δχ²)",
       x = "Variable Removed from Model",
       y = "Δ Log-Likelihood (Chi²)") +
  theme_minimal(base_size = 14)

#Rule of law is the strongest predictors, We include it in the other models
model_base <-multinom(cluster ~ ruleoflaw, data = df_model,trace=TRUE)
model_col  <- multinom(cluster ~ ruleoflaw+colonizer, data = df_model,trace=TRUE)
model_set  <- multinom(cluster ~ ruleoflaw+lcapped, data = df_model,trace=TRUE)
model_lat  <- multinom(cluster ~ ruleoflaw+lat_abst, data = df_model,trace=TRUE, maxit=200)
model_pop <- multinom(cluster ~ ruleoflaw+lpd1500s, data = df_model,trace=TRUE)
model_prot <- multinom(cluster ~ ruleoflaw+protmiss, data = df_model, trace=TRUE)
model_school <- multinom(cluster ~ ruleoflaw+prienr1900, data = df_model, trace=TRUE)

#Model_lat did not converge, remove and rerun evaluations

var_contrib <- tibble::tibble(
  Variable = c("Colonizer","Settler mortality","Latitude" ,"Pop. density (1500)",
               "Protestant missions", "School enrollment (1900)"),
  DeltaChi2 = c(
    delta_logLik(model_col),
    delta_logLik(model_set),
    delta_logLik(model_lat),
    delta_logLik(model_pop),
    delta_logLik(model_prot),
    delta_logLik(model_school)
  )
  
)
# Plot
ggplot(var_contrib, aes(x = reorder(Variable, DeltaChi2), y = DeltaChi2)) +
  geom_col(fill = "#3366CC") +
  coord_flip() +
  labs(title = "Variable Contributions to Cluster Prediction",
       subtitle = "Based on Likelihood Ratio Tests (Δχ²)",
       x = "Variable Removed from Model",
       y = "Δ Log-Likelihood (Chi²)") +
  theme_minimal(base_size = 14)

#Colonizer is the best predictor
model_base <- multinom(cluster ~ ruleoflaw + colonizer, data = df_model,trace=TRUE)
model_set  <- multinom(cluster ~ ruleoflaw + colonizer + lcapped, data = df_model,trace=TRUE,maxit=200)
model_lat  <- multinom(cluster ~ ruleoflaw + colonizer + lat_abst, data = df_model,trace=TRUE,maxit=200)
model_pop <- multinom(cluster ~ ruleoflaw + colonizer + lpd1500s, data = df_model,trace=TRUE,maxit=200)
model_prot <- multinom(cluster ~ ruleoflaw + colonizer + protmiss, data = df_model, trace=TRUE,maxit=200)
model_school <- multinom(cluster ~ ruleoflaw + colonizer + prienr1900, data = df_model, trace=TRUE,maxit=200)


var_contrib <- tibble::tibble(
  Variable = c("Settler mortality","Latitude" ,"Pop. density (1500)",
               "Protestant missions", "School enrollment (1900)"),
  DeltaChi2 = c(
    delta_logLik(model_set),
    delta_logLik(model_lat),
    delta_logLik(model_pop),
    delta_logLik(model_prot),
    delta_logLik(model_school)
  )
  
)
# Plot
ggplot(var_contrib, aes(x = reorder(Variable, DeltaChi2), y = DeltaChi2)) +
  geom_col(fill = "#3366CC") +
  coord_flip() +
  labs(title = "Variable Contributions to Cluster Prediction",
       subtitle = "Based on Likelihood Ratio Tests (Δχ²)",
       x = "Variable Removed from Model",
       y = "Δ Log-Likelihood (Chi²)") +
  theme_minimal(base_size = 14)

#Try adding school enrolment
model_large <- update(model_base,~.+prienr1900,maxit=200)

summary(model_large)
car::Anova(model_large,type=2)

#The final model is
model_final <- multinom(cluster ~ ruleoflaw + colonizer, data = df_model,trace=TRUE)
summary(model_final)
car::Anova(model_final,type=2)


#Visualize model--------------------

# Make prediction grid
newdata <- expand.grid(
  colonizer = levels(df_model$colonizer),
  ruleoflaw = seq(min(df_model$ruleoflaw), max(df_model$ruleoflaw), length.out = 100)
)

# Get predicted probabilities
preds <- predict(model_final, newdata = newdata, type = "probs")

# Convert to long format
pred_df <- cbind(newdata, preds) %>%
  tidyr::pivot_longer(cols = -c(colonizer,ruleoflaw),
                      names_to = "cluster", values_to = "probability")%>%
  unique()

# Plot
ggplot(pred_df, aes(x = ruleoflaw, y = probability, color = cluster)) +
  geom_line(size = 1.2) +
  geom_point(data=df_model, aes(x=ruleoflaw, y=1, color=cluster),alpha=0.5,size=2)+
  facet_wrap(~colonizer) +
  theme_minimal() +
  labs(title = "Predicted Probabilities by Colonizer and Rule of law")

#test accuracy
predictions <- predict(model_final, newdata = df_model)
correctness <- as.numeric(df_model$cluster) - as.numeric(predictions)==0
accuracy <- sum(correctness) / length(correctness)
print(accuracy)





#test using only historical variables---------------------
#Base model
model_base <- multinom(cluster ~ 1, data = df_model,
                       trace=TRUE)


# Model without one variable at a time:
model_col  <- multinom(cluster ~ colonizer, data = df_model,trace=TRUE)
model_rule  <- multinom(cluster ~ ruleoflaw, data = df_model,trace=TRUE)
model_set  <- multinom(cluster ~ lcapped, data = df_model,trace=TRUE)
model_lat  <- multinom(cluster ~ lat_abst, data = df_model,trace=TRUE)
model_pop <- multinom(cluster ~ lpd1500s, data = df_model,trace=TRUE)
model_prot <- multinom(cluster ~ protmiss, data = df_model, trace=TRUE)
model_school <- multinom(cluster ~ prienr1900, data = df_model, trace=TRUE)


# Compare log-likelihoods
delta_logLik <- function(model_big) {
  2 * (logLik(model_big)-logLik(model_base))
}

var_contrib <- tibble::tibble(
  Variable = c("Colonizer","Rule of law (2005)","Settler mortality", "Latitude", "Pop. density (1500)",
               "Protestant missions", "School enrollment (1900)"),
  DeltaChi2 = c(
    delta_logLik(model_col),
    delta_logLik(model_rule),
    delta_logLik(model_set),
    delta_logLik(model_lat),
    delta_logLik(model_pop),
    delta_logLik(model_prot),
    delta_logLik(model_school)
  )
  
)

# Plot
ggplot(var_contrib, aes(x = reorder(Variable, DeltaChi2), y = DeltaChi2)) +
  geom_col(fill = "#3366CC") +
  coord_flip() +
  labs(title = "Variable Contributions to Cluster Prediction",
       subtitle = "Based on Likelihood Ratio Tests (Δχ²)",
       x = "Variable Removed from Model",
       y = "Δ Log-Likelihood (Chi²)") +
  theme_minimal(base_size = 14)

#School enrollment is best historical predictor
model_base <-multinom(cluster ~ prienr1900, data = df_model,trace=TRUE)
model_col  <- multinom(cluster ~ prienr1900+colonizer, data = df_model,trace=TRUE)
model_set  <- multinom(cluster ~ prienr1900+lcapped, data = df_model,trace=TRUE)
model_lat  <- multinom(cluster ~ prienr1900+lat_abst, data = df_model,trace=TRUE, maxit=200)
model_pop <- multinom(cluster ~ prienr1900+lpd1500s, data = df_model,trace=TRUE)
model_prot <- multinom(cluster ~ prienr1900+protmiss, data = df_model, trace=TRUE)

var_contrib <- tibble::tibble(
  Variable = c("Colonizer","Settler mortality","Latitude" ,"Pop. density (1500)",
               "Protestant missions"),
  DeltaChi2 = c(
    delta_logLik(model_col),
    delta_logLik(model_set),
    delta_logLik(model_lat),
    delta_logLik(model_pop),
    delta_logLik(model_prot)
  )
  
)
# Plot
ggplot(var_contrib, aes(x = reorder(Variable, DeltaChi2), y = DeltaChi2)) +
  geom_col(fill = "#3366CC") +
  coord_flip() +
  labs(title = "Variable Contributions to Cluster Prediction",
       subtitle = "Based on Likelihood Ratio Tests (Δχ²)",
       x = "Variable Removed from Model",
       y = "Δ Log-Likelihood (Chi²)") +
  theme_minimal(base_size = 14)


#settler mortality is best historical predictor
model_base <-multinom(cluster ~ prienr1900 + lcapped, data = df_model,trace=TRUE)
model_col  <- multinom(cluster ~ prienr1900 + lcapped+colonizer, data = df_model,trace=TRUE,maxit=200)
model_lat  <- multinom(cluster ~ prienr1900 + lcapped+lat_abst, data = df_model,trace=TRUE, maxit=200)
model_pop <- multinom(cluster ~ prienr1900 + lcapped+lpd1500s, data = df_model,trace=TRUE)
model_prot <- multinom(cluster ~ prienr1900 + lcapped+protmiss, data = df_model, trace=TRUE)

var_contrib <- tibble::tibble(
  Variable = c("Colonizer","Latitude" ,"Pop. density (1500)",
               "Protestant missions"),
  DeltaChi2 = c(
    delta_logLik(model_col),
    delta_logLik(model_lat),
    delta_logLik(model_pop),
    delta_logLik(model_prot)
  )
  
)
# Plot
ggplot(var_contrib, aes(x = reorder(Variable, DeltaChi2), y = DeltaChi2)) +
  geom_col(fill = "#3366CC") +
  coord_flip() +
  labs(title = "Variable Contributions to Cluster Prediction",
       subtitle = "Based on Likelihood Ratio Tests (Δχ²)",
       x = "Variable Removed from Model",
       y = "Δ Log-Likelihood (Chi²)") +
  theme_minimal(base_size = 14)


model_large <- update(model_base,~.+colonizer, maxit=200)
summary(model_large)
car::Anova(model_large,type=2) #Colonizer is not significant

model_final <- model_base
summary(model_final)


predictions <- predict(model_final, newdata = df_model)
correctness <- as.numeric(df_model$cluster) - as.numeric(predictions)==0
accuracy <- sum(correctness) / length(correctness)
print(accuracy)

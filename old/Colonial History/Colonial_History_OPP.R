
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
df_model <- df %>% select(cluster = "Only_OPP" , colonizer, lcapped, lpd1500s, lat_abst,
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
    title = "Counts of B within each A",
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

#colonizer is the strongest predictors, We include it in the other models
model_base <-multinom(cluster ~ prienr1900, data = df_model,trace=TRUE)
model_rule  <- multinom(cluster ~ ruleoflaw+prienr1900, data = df_model,trace=TRUE)
model_set  <- multinom(cluster ~ prienr1900+lcapped, data = df_model,trace=TRUE)
model_lat  <- multinom(cluster ~ prienr1900+lat_abst, data = df_model,trace=TRUE, maxit=200)
model_pop <- multinom(cluster ~ prienr1900+lpd1500s, data = df_model,trace=TRUE)
model_prot <- multinom(cluster ~ prienr1900+protmiss, data = df_model, trace=TRUE)
model_col <- multinom(cluster ~ colonizer+prienr1900, data = df_model, trace=TRUE)

#
var_contrib <- tibble::tibble(
  Variable = c("RuleOflaw","Settler mortality","Latitude" ,"Pop. density (1500)",
               "Protestant missions", "Colonizer"),
  DeltaChi2 = c(
    delta_logLik(model_rule),
    delta_logLik(model_set),
    delta_logLik(model_lat),
    delta_logLik(model_pop),
    delta_logLik(model_prot),
    delta_logLik(model_col)
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

#Rule of law is the best predictor
model_base <- multinom(cluster ~ ruleoflaw + prienr1900, data = df_model,trace=TRUE)
model_set  <- multinom(cluster ~ ruleoflaw + prienr1900 + lcapped, data = df_model,trace=TRUE,maxit=200)
model_lat  <- multinom(cluster ~ ruleoflaw + prienr1900 + lat_abst, data = df_model,trace=TRUE,maxit=200)
model_pop <- multinom(cluster ~ ruleoflaw + prienr1900 + lpd1500s, data = df_model,trace=TRUE,maxit=200)
model_prot <- multinom(cluster ~ ruleoflaw + prienr1900 + protmiss, data = df_model, trace=TRUE,maxit=200)
model_col <- multinom(cluster ~ ruleoflaw + colonizer + prienr1900, data = df_model, trace=TRUE,maxit=200)


var_contrib <- tibble::tibble(
  Variable = c("Settler mortality","Latitude" ,"Pop. density (1500)",
               "Protestant missions", "Colonizer"),
  DeltaChi2 = c(
    delta_logLik(model_set),
    delta_logLik(model_lat),
    delta_logLik(model_pop),
    delta_logLik(model_prot),
    delta_logLik(model_col)
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

#Try adding population density
model_large <- update(model_base,~.+colonizer,maxit=200)
exp((AIC(model_large)-AIC(model_base))/2) #Good chance for improvement with large model

#The final model is the base because the large model yielded NaN estimates
model_final <- model_base 
summary(model_final)

z <- summary(model_final)$coefficients / summary(model_final)$standard.errors
p <- 2*(1-pnorm(abs(z)))
print(p)

#Visualize model------------------------------------------------------------------------

# Make prediction grid
newdata <- expand.grid(
  prienr1900 = seq(min(df_model$prienr1900), max(df_model$prienr1900), length.out = 100),
  ruleoflaw = seq(min(df_model$ruleoflaw), max(df_model$ruleoflaw), length.out = 100)
)

# Get predicted probabilities
preds <- predict(model_final, newdata = newdata, type = "probs")

# Convert to long format
pred_df <- cbind(newdata, preds) %>%
  tidyr::pivot_longer(cols = -c(prienr1900,ruleoflaw),
                      names_to = "cluster", values_to = "probability")%>%
  unique()

# Plot
ggplot() +
  geom_tile(data=pred_df, aes(x = ruleoflaw, y = prienr1900, fill=probability)) +
  geom_point(data=df_model, aes(x=ruleoflaw, y = prienr1900))+
  facet_wrap(~cluster) +
  theme_minimal() +
  labs(title = "Predicted Probabilities by Colonizer and Rule of law")


#test accuracy
predictions <- predict(model_final, newdata = df_model)
correctness <- as.numeric(df_model$cluster) - as.numeric(predictions)==0
accuracy <- sum(correctness) / length(correctness)
print(accuracy)

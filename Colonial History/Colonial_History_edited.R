
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

#Initial model
model1 <- multinom(cluster~colonizer, df_model)

z <- summary(model1)$coefficients / summary(model1)$standard.errors
p <- 2 * (1 - pnorm(abs(z)))
print(p)

#Add ruleoflaw 
model2 <- multinom(cluster~colonizer+ruleoflaw, df_model)

#Check degrees of freedom
model1$edf
model2$edf

# ========== Predict probabilities by colonizer ==========

# Create synthetic observation per colonizer, with average values of other variables
mean_vars <- df_model %>%
  summarise(across(c(lcapped, lpd1500s, lat_abst, ruleoflaw, protmiss, prienr1900), mean))

newdata <- expand.grid(
  colonizer = levels(df_model$colonizer),
  lcapped   = mean_vars$lcapped,
  lpd1500s  = mean_vars$lpd1500s,
  lat_abst  = mean_vars$lat_abst,
  ruleoflaw = mean_vars$ruleoflaw,
  protmiss  = mean_vars$protmiss,
  prienr1900= mean_vars$prienr1900
)

# Predict probabilities
pred_probs <- predict(model1, newdata = newdata, type = "probs")
pred_df <- cbind(newdata["colonizer"], as.data.frame(pred_probs))

# Pivot for plotting
library(tidyr)
plot_df <- pred_df %>%
  pivot_longer(cols = -colonizer, names_to = "cluster", values_to = "probability")

# Plot
ggplot(plot_df, aes(x = cluster, y = probability, fill = colonizer)) +
  geom_col(position = "dodge") +
  labs(title = "Predicted Cluster Probabilities by Colonizer",
       y = "Predicted Probability", x = "Cluster") +
  theme_minimal()


#Base model
model_base <- multinom(cluster ~ colonizer, data = df_model,
                       trace=TRUE)

# Model without one variable at a time:
model_nocol  <- multinom(cluster ~ lcapped + lpd1500s + lat_abst + ruleoflaw + protmiss + prienr1900, data = df_model)
model_noset  <- multinom(cluster ~ colonizer + lpd1500s + lat_abst + ruleoflaw + protmiss + prienr1900, data = df_model)
model_nopop  <- multinom(cluster ~ colonizer + lcapped + lat_abst + ruleoflaw + protmiss + prienr1900, data = df_model)
model_nolat  <- multinom(cluster ~ colonizer + lcapped + lpd1500s + ruleoflaw + protmiss + prienr1900, data = df_model)
model_norule <- multinom(cluster ~ colonizer + lcapped + lpd1500s + lat_abst + protmiss + prienr1900, data = df_model)
model_nomiss <- multinom(cluster ~ colonizer + lcapped + lpd1500s + lat_abst + ruleoflaw + prienr1900, data = df_model)
model_noschool <- multinom(cluster ~ colonizer + lcapped + lpd1500s + lat_abst + ruleoflaw + protmiss, data = df_model)

# Compare log-likelihoods
ll_base <- logLik(model_base)

delta_logLik <- function(model_small) {
  2 * (ll_full - logLik(model_small))
}

delta_logLik <- function(model_big) {
  2 * (logLik(model_big)-ll_base)
}

# Degrees of freedom = number of parameters removed
cat("Likelihood Ratio Tests:\n")
cat("No colonizer:   Δχ² =", round(delta_logLik(model_nocol), 2), "\n")
cat("No lcapped:     Δχ² =", round(delta_logLik(model_noset), 2), "\n")
cat("No lpd1500s:    Δχ² =", round(delta_logLik(model_nopop), 2), "\n")
cat("No latitude:    Δχ² =", round(delta_logLik(model_nolat), 2), "\n")
cat("No ruleoflaw:   Δχ² =", round(delta_logLik(model_norule), 2), "\n")
cat("No protmiss:    Δχ² =", round(delta_logLik(model_nomiss), 2), "\n")
cat("No school1900:  Δχ² =", round(delta_logLik(model_noschool), 2), "\n")





# Compute Δχ² values for each nested model (removing one variable at a time)
delta_logLik <- function(model_small) {
  2 * (logLik(model) - logLik(model_small))
}

var_contrib <- tibble::tibble(
  Variable = c("Colonizer", "Settler mortality", "Pop. density (1500)", "Latitude",
               "Rule of law (2005)", "Protestant missions", "School enrollment (1900)"),
  DeltaChi2 = c(
    delta_logLik(model_nocol),
    delta_logLik(model_noset),
    delta_logLik(model_nopop),
    delta_logLik(model_nolat),
    delta_logLik(model_norule),
    delta_logLik(model_nomiss),
    delta_logLik(model_noschool)
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

#write.csv(df, "df")


model <- multinom(cluster~colonizer+lcapped,
                  df_model, trace=TRUE)
model$edf
modelSum <- summary(model)
z <- modelSum$coefficients / modelSum$standard.errors
p <- 2*(1-pnorm(abs(z)))
print(p)


# Make prediction grid
newdata <- expand.grid(
  colonizer = levels(df_model$colonizer),
  #ruleoflaw = seq(min(df_model$ruleoflaw), max(df_model$ruleoflaw), length.out = 100),
  lcapped = seq(min(df_model$lcapped), max(df_model$lcapped), length.out = 100),
  lat_abst = min(df_model$lat_abst)
  )

# Get predicted probabilities
preds <- predict(model, newdata = newdata, type = "probs")

# Convert to long format
pred_df <- cbind(newdata, preds) %>%
  tidyr::pivot_longer(cols = -c(colonizer,lcapped,lat_abst),
                      names_to = "cluster", values_to = "probability")%>%
  unique()

# Plot
ggplot(pred_df, aes(x = lcapped, y = probability, color = cluster)) +
  geom_line(size = 1.2) +
  facet_wrap(~colonizer) +
  theme_minimal() +
  labs(title = "Predicted Probabilities by Colonizer and Rule of law")






#Iterative model tests -------------------------------------------
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

models <- list(model_base,
               model_col,model_rule,
               model_set,model_lat,
               model_pop,model_prot,
               model_school)

# Compare log-likelihoods
delta_logLik <- function(model_big) {
  2 * (logLik(model_big)-logLik(model_base))
}

var_contrib <- tibble::tibble(
  Variable = c("Colonizer","Rule of law (2005)","Settler mortality", "Latitude", "Pop. density (1500)",
               "Protestant missions", "School enrollment (1900)"),
  
  DeltaChi2 = c(lapply(1:length(models),
                       function(x) delta_logLik(models[[x]]))
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


#Rule of law is the strongest predictor
for (i in 1:length(models)){
  models[[i]] <- update(models[[i]],~.+ruleoflaw)
}


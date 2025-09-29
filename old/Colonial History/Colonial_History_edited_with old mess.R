
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

# Subset needed variables
df_model <- df %>% select(cluster = `Few indicators_SPI_preferred_90` , colonizer, lcapped, lpd1500s, lat_abst,
                          ruleoflaw, protmiss, prienr1900) %>%
  na.omit()

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

ggplot(df_model, aes(x = cluster, y = ruleoflaw, color = colonizer)) +
  geom_boxplot() +
  labs(
    title = "Counts of B within each A",
    x = "Cluster",
    y = "Count",
    fill = "Colonizer"
  ) +
  theme_minimal()

#Chisq test to see if groupings are significantly divided
cluCol <- table(df_model$cluster,df$colonizer)
totalClu<- rowSums(cluCol)
totalCol <- colSums(cluCol)
total <- sum(cluCol)
eMatrix <- matrix(0, nrow=nrow(cluCol),ncol=ncol(cluCol))

for (i in 1:length(totalClu)){
  for (j in 1:length(totalCol)){
    eMatrix[i,j] <- totalClu[i]*totalCol[j] / total
  }
}
#Since the elements of the matrix are less than 5, the chisq is unlikely to be valid
#We can run a simulation, which indicates a significant relations ship  between cluster and colonizer
chisq.test(cluCol, sim=TRUE,  B=20000 ) 
mosaicplot(cluCol, shade=TRUE)

#Check data with PCA
pca_result <- prcomp(df_model%>%
                       select(is.numeric),scale=TRUE,center = TRUE)

#Scree plot                    
plot(pca_result$sdev^2/sum(pca_result$sdev^2))

#Plotting for the first three dimensions
pc_scores <- cbind(df_model%>%select(!is.numeric),pca_result$x[, 1:3])

library(patchwork)
filter <- c("Spanish","British","French","Other")

p1 <- ggplot(data=pc_scores %>%filter(colonizer %in% filter) ,
             aes(x=PC1,y=PC2, colour = cluster,shape=colonizer))+
  geom_point()
p2 <- ggplot(data=pc_scores%>%filter(colonizer %in% filter)
             , aes(x=PC2,y=PC3, colour = cluster,shape=colonizer))+
  geom_point()
p3 <- ggplot(data=pc_scores%>%filter(colonizer %in% filter)
             , aes(x=PC1,y=PC3, colour = cluster,shape=colonizer))+
  geom_point()
p1/p3+p2




#For our case a
model <- multinom(cluster~colonizer, df_model)
table(df_model$cluster, df_model$colonizer)

z <- summary(model)$coefficients / summary(model)$standard.errors
p <- 2 * (1 - pnorm(abs(z)))
print(p)

# Fit multinomial logit model
model <- multinom(cluster ~ colonizer + lcapped + lpd1500s + lat_abst +
                    ruleoflaw + protmiss + prienr1900,
                  data = df_model, trace = FALSE)

#Summary
summary(model)

# P-values
z <- summary(model)$coefficients / summary(model)$standard.errors
p <- 2 * (1 - pnorm(abs(z)))
print(round(p, 4))

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
pred_probs <- predict(model, newdata = newdata, type = "probs")
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


# Full model
model_full <- multinom(cluster ~ colonizer + lcapped + lpd1500s + lat_abst +
                         ruleoflaw + protmiss + prienr1900, data = df_model)

# Model without one variable at a time:
model_nocol  <- multinom(cluster ~ lcapped + lpd1500s + lat_abst + ruleoflaw + protmiss + prienr1900, data = df_model)
model_noset  <- multinom(cluster ~ colonizer + lpd1500s + lat_abst + ruleoflaw + protmiss + prienr1900, data = df_model)
model_nopop  <- multinom(cluster ~ colonizer + lcapped + lat_abst + ruleoflaw + protmiss + prienr1900, data = df_model)
model_nolat  <- multinom(cluster ~ colonizer + lcapped + lpd1500s + ruleoflaw + protmiss + prienr1900, data = df_model)
model_norule <- multinom(cluster ~ colonizer + lcapped + lpd1500s + lat_abst + protmiss + prienr1900, data = df_model)
model_nomiss <- multinom(cluster ~ colonizer + lcapped + lpd1500s + lat_abst + ruleoflaw + prienr1900, data = df_model)
model_noschool <- multinom(cluster ~ colonizer + lcapped + lpd1500s + lat_abst + ruleoflaw + protmiss, data = df_model)

# Compare log-likelihoods
ll_full <- logLik(model_full)

delta_logLik <- function(model_small) {
  2 * (ll_full - logLik(model_small))
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

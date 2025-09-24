
library(readxl)
library(readr)
library(dplyr)
library(nnet)
library(ggplot2)
library(patchwork)

# Load data
agr <- read_excel("xcountry_data_full.xlsx")
names(agr)[names(agr) == 'code'] <- 'iso3'

#New clusters from experiments
clusters <- readRDS("clusterVariations_laglead_scaled.RDS")%>%
  rename(iso3 = SPI_countrycode)

# Clean ISO3 values
clusters <- clusters %>% mutate(iso3 = toupper(trimws(as.character(iso3))))
agr      <- agr      %>% mutate(iso3 = toupper(trimws(as.character(iso3))))

# Merge
df <- clusters %>% inner_join(agr, by = "iso3")
#write.csv(df,"dfColonial.csv")

# Simplified colonizer grouping
df$colonizer <- "Other"
df$colonizer[df$f_brit == 1]   <- "British"
df$colonizer[df$f_french == 1] <- "French"
df$colonizer[df$f_spain == 1]  <- "Spanish"
df$colonizer <- factor(df$colonizer)

# Make clusters a factor
df <- df %>%
  mutate(across(3:10,as.factor))

colnames(df)

clusterSelection = "Few indicators_SPI_preferred"
# Subset needed variables
df_model <- df %>% select(iso3,cluster = clusterSelection , colonizer, lcapped, lpd1500s, lat_abst,
                          ruleoflaw, protmiss, prienr1900) %>%
  na.omit()%>%
  mutate(cluster=factor(cluster,levels=c(10,4,1,8,9,7,3,2,5,11,6)))


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


#Plot dimensions
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
model_set  <- multinom(cluster ~ lcapped, data = df_model,trace=TRUE,maxit=500)
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

#settler mortality is the strongest predictors, We include it in the other models
model_base <-multinom(cluster ~ lcapped, data = df_model,trace=TRUE, maxit=500)
model_school  <- multinom(cluster ~ lcapped+prienr1900, data = df_model,trace=TRUE, maxit=500)
model_col  <- multinom(cluster ~ lcapped+colonizer, data = df_model,trace=TRUE,maxit=500)
model_lat  <- multinom(cluster ~ lcapped+lat_abst, data = df_model,trace=TRUE, maxit=500)
model_pop <- multinom(cluster ~ lcapped+lpd1500s, data = df_model,trace=TRUE, maxit=500)
model_prot <- multinom(cluster ~ lcapped+protmiss, data = df_model, trace=TRUE, maxit=500)
model_rule <- multinom(cluster ~ lcapped+ruleoflaw, data = df_model, trace=TRUE, maxit=500)

#
var_contrib <- tibble::tibble(
  Variable = c("School enrollment 1900","Colonizer","Latitude" ,"Pop. density (1500)",
               "Protestant missions", "Rule of law"),
  DeltaChi2 = c(
    delta_logLik(model_school),
    delta_logLik(model_col),
    delta_logLik(model_lat),
    delta_logLik(model_pop),
    delta_logLik(model_prot),
    delta_logLik(model_rule)
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
model_base <- multinom(cluster ~ lcapped + colonizer, data = df_model,trace=TRUE,maxit=500)
model_rule  <- multinom(cluster ~ lcapped + colonizer + ruleoflaw, data = df_model,trace=TRUE,maxit=500)
model_lat  <- multinom(cluster ~ lcapped + colonizer + lat_abst, data = df_model,trace=TRUE,maxit=10000)
model_pop <- multinom(cluster ~ lcapped + colonizer + lpd1500s, data = df_model,trace=TRUE,maxit=500)
model_prot <- multinom(cluster ~ lcapped + colonizer + protmiss, data = df_model, trace=TRUE,maxit=500)
model_school <- multinom(cluster ~ lcapped + colonizer + prienr1900, data = df_model, trace=TRUE,maxit=200)


var_contrib <- tibble::tibble(
  Variable = c("Ruleoflaw","Latitude" ,"Pop. density (1500)",
               "Protestant missions", "School enrollment 1900"),
  DeltaChi2 = c(
    delta_logLik(model_rule),
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

#Try adding school enrollment
model_large <- update(model_base,~.+lpd1500s,maxit=500)
summary(model_large)

#Type two anova because we don't have an interaction
car::Anova(model_large, type="2")

#The final model is the large model as all variables are significant
model_final <- update(model_large,~.-colonizer)
summary(model_final)
car::Anova(model_final, type="2")

#What if we replace ruleoflaw with lcapped? 
car::Anova(update(model_large,~.-colonizer), type="2") #also significant


z <- summary(model_final)$coefficients / summary(model_final)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)

#Visualize model------------------------------------------------------------------------
# Make prediction grid
newdata <- expand.grid(
  colonizer = levels(df_model$colonizer),
  lpd1500s = seq(min(df_model$lpd1500s), max(df_model$lpd1500s), length.out = 100),
  lcapped = seq(min(df_model$lcapped), max(df_model$lcapped), length.out = 100)
)

# Get predicted probabilities
preds <- predict(model_final, newdata = newdata, type = "probs")

# Convert to long format
pred_df <- cbind(newdata, preds) %>%
  tidyr::pivot_longer(cols = -c(colonizer,lcapped,protmiss),
                      names_to = "cluster", values_to = "probability")%>%
  unique()

# Plot
ggplot(pred_df%>%filter(protmiss==min(protmiss)), aes(x = lcapped, y = probability, color = cluster)) +
  geom_line(size = 1.2) +
  geom_point(data=df_model, aes(x=lcapped, y=1, color=cluster),alpha=0.5,size=2)+
  facet_wrap(~colonizer) +
  theme_minimal() +
  labs(title = "Predicted Probabilities by Colonizer and Rule of law")

#test accuracy
summary(model_final)
predictions <- predict(model_final, newdata = df_model)
correctness <- as.numeric(df_model$cluster) - as.numeric(predictions)==0
accuracy <- sum(correctness) / length(correctness)
print(accuracy)



#Dual axis plot------------------------------------------------------------
df_dual <- df %>%
  select(cluster = SPI_baseline,
         logpgdp05,
         lcapped,
         colonizer)%>%
  group_by(cluster)%>%
  mutate(across(c(logpgdp05, lcapped), mean, na.rm = TRUE))%>%
  add_count()

coeff <- 1.7

ggplot(df_dual, aes(x = cluster)) +
  geom_col(aes(y = lcapped/n, fill = colonizer)) +
  geom_text(aes(y = lcapped + 0.2, label=n))+
  geom_line(aes(y = logpgdp05/coeff, group = 1), color = "black") +
  scale_y_continuous(
    name = "log settler mortality",
    sec.axis = sec_axis(~.*coeff, name = "log(GDP per capita) (2005)")
  )+
  theme_minimal()+
  labs(title="")






#Test historical model---------------------------------------------
model_test <- multinom(cluster ~  lcapped + lpd1500s,
                        data=df_model,
                        maxit=500)

z <- summary(model_test)$coefficients / summary(model_test)$standard.errors
p <- 2*(1-pnorm(abs(z),0,1))
print(p)
car::Anova(model_test,type=2)

# Make prediction grid
newdata <- expand.grid(
  lpd1500s = seq(min(df_model$lpd1500s), max(df_model$lpd1500s), length.out = 50),
  lcapped = seq(min(df_model$lcapped), max(df_model$lcapped), length.out = 50)
)

# Get predicted probabilities
preds <- predict(model_test, newdata = newdata, type = "probs")

# Convert to long format
pred_df <- cbind(newdata, preds) %>%
  tidyr::pivot_longer(cols = -c(lcapped,lpd1500s),
                      names_to = "cluster", values_to = "probability")%>%
  unique()

pred_df_max <- pred_df %>%
  group_by(lcapped,lpd1500s) %>%
  slice_max(probability, with_ties = FALSE) %>%
  ungroup()

# Plot
ggplot(pred_df_max, aes(x = lcapped, y = lpd1500s, fill = cluster)) +
  geom_tile() +
  geom_point(data=df_model, aes(x=lcapped, y=lpd1500s, fill=factor(cluster)), color="black",pch=21, alpha=0.5,size=3)+
  #facet_wrap(~colonizer) +
  theme_minimal() +
  labs(title = "Predicted Probabilities by Colonizer, Settler mortality and school enrollment (1900)")

#test accuracy
summary(model_test)
predictions <- predict(model_test, newdata = df_model)
correctness <- as.numeric(df_model$cluster) - as.numeric(predictions)==0
accuracy <- sum(correctness) / length(correctness)
print(accuracy)


predictions <- predict(multinom(cluster~1,data=df_model), newdata = df_model)
correctness <- as.numeric(df_model$cluster) - as.numeric(predictions)==0
accuracy <- sum(correctness) / length(correctness)
print(accuracy)
#Test with dls data from subnational survey article-----
load(file.choose())
dls_latest  <- dls_country %>%
  group_by(country_name) %>%
  slice_max(order_by = year_start_interviews, n = 1)%>%
  rename(Country = country_name)

df_dls <- df_model%>%
  left_join(dls_latest, by="Country")

fit <- lm(dls_index_country ~ colonizer+lcapped+prienr1900,
          data=df_dls)
summary(fit)

#Bad overlap between the two data sets
sum(is.na(df_dls$dls_index_country))/nrow(df_dls)

#Run ordinal model ----------------------------------------------------------------
library(ordinal)

fit1 <- clm(cluster~ colonizer + lcapped + protmiss + lpd1500s+
             prienr1900+ lat_abst, data=df_model)
anova(fit1)
fit2 <- update(fit1,~.-lat_abst)
anova(fit2)

fit3 <- update(fit2,~.-protmiss)
anova(fit3)


#Predict
newdata <- data.frame(x1 = c(-1, 0, 1), x2 = "A")
predict(model, newdata, type = "prob")


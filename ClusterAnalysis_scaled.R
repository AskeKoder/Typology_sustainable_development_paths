#This script analyses the clusters resulting from using lag lead imputation on scaled indicators

library(dplyr)

#Well being data
# data <- read.csv("extendedDatascaled.csv")%>%
#   select(-X)

data <- read.csv("ImputedDataLag1Lead2_maxit30_scaled.csv")%>%
  select(-c("X",".id"))%>%
  relocate(.imp, .after=last_col())

#Clusters from experiments
clusterSelection <- "Few indicators_SPI_preferred"
clusters <- readRDS("clusterVariations_laglead_scaled.RDS")%>%
  select(Country,SPI_countrycode,cluster = clusterSelection)

#Read experiment file
Experiments <- read_xlsx("IndicatorsForClusters.xlsx",sheet="Experiments")



#link experiments with included indicators
batch_scaled <- Experiments%>%
  select(c("Baseline",
           "BHN_FWB",
           "No_non_DLS",
           "No_nonDLS_limExt",
           "Few indicators_SPI_preferred",
           "Few indicators_SPI_preferred_90",
           "Few_indicators_closest_DLS_coverage",
           "Few_indicators_closest_DLS_coverage_90"
  ))%>%
  as.matrix()%>%
  na.omit()

rownames(batch_scaled) <- colnames(data)[5:(ncol(data)-1)]
rownames(batch_scaled) <- rownames(batch_scaled)[c(6,5,4,3,2,1,
                                                   10,9,8,7,
                                                   11,12,13,
                                                   18,17,16,15,14,
                                                   22,21,20,19,
                                                   23,24,25,
                                                   29,28,27,26,
                                                   33,32,31,30,
                                                   39,38,37,36,35,34,
                                                   44,43,42,41,40,
                                                   48,47,46,45,
                                                   52,51,50,49,
                                                   53,55,57,56,
                                                   54,60,58,59,
                                                   61,62,63)]

batch_scaled[,c(1:ncol(batch_scaled))] <- as.numeric(batch_scaled[,c(1:ncol(batch_scaled))])
batch_scaled

#Merge data sets
clusteredData <- data %>%
  left_join(clusters , by= c("Country","SPI_countrycode"))%>%
  data.frame()%>%
  mutate(cluster = as.factor(cluster))


#Plot -----------------------------------------------
selected_vars <- rownames(batch_scaled)[batch_scaled[, clusterSelection] == 1]
set <- clusteredData %>%
  group_by(SPI_year, cluster) %>%
  dplyr::select(c("SPI_year", "cluster",
                  selected_vars))%>%
  summarise(across(all_of(selected_vars),
                   list(
                     median = ~median(., na.rm = TRUE),
                     q1 = ~quantile(., 0.25, na.rm = TRUE),
                     q3 = ~quantile(., 0.75, na.rm = TRUE)
                   ),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")
  

#Get names of variables
variables <- selected_vars


#Convert to long format for plotting
long_median <- set %>%
  select(1,2, ends_with("_median"))%>%
  pivot_longer(cols = ends_with("_median"),
               names_to = "Variable",
               values_to = "Median",
               names_pattern = "(.*)_median")
long_q1 <- set %>%
  select(1,2, ends_with("_q1"))%>%
  pivot_longer(cols = ends_with("_q1"),
               names_to = "Variable",
               values_to = "Q1",
               names_pattern = "(.*)_q1")
long_q3 <- set %>%
  select(1,2, ends_with("_q3"))%>%
  pivot_longer(cols = ends_with("_q3"),
               names_to = "Variable",
               values_to = "Q3",
               names_pattern = "(.*)_q3")
# Combine all
long_set <- long_median %>%
  left_join(long_q1, by = c("SPI_year", "cluster", "Variable")) %>%
  left_join(long_q3, by = c("SPI_year", "cluster", "Variable")) %>%
  mutate(Variable = factor(Variable, levels = variables))


ggplot(long_set, aes(x = SPI_year, y = Median, color = cluster, fill = cluster)) +
  geom_ribbon(aes(ymin = Q1, ymax = Q3), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Variable, scales = "free_y") +
  ylim(0, 100) +
  theme_minimal() +
  # scale_fill_manual(values = cluster_colors)+
  # scale_color_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Year", y = "Cluster Median (with IQR)")+
  theme(plot.title = element_text(size = 5))


long_median_ranked <- long_median %>%
  group_by(SPI_year, Variable) %>%  # Group by year and variable (NOT group)
  mutate(Rank = case_when(
    Variable == "Share_Slums" ~ rank(-Median, ties.method = "min"),   # lower = better
    TRUE                       ~ rank(Median, ties.method = "min")   # higher = better
  ))%>%
  ungroup()

ggplot(long_median_ranked, aes(x = SPI_year, y = Rank, color = cluster, fill = cluster)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Variable, scales = "free_y") +
  ylim(0, 11) +
  theme_minimal() +
  # scale_fill_manual(values = cluster_colors)+
  # scale_color_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Year", y = "Cluster Median (with IQR)")+
  theme(plot.title = element_text(size = 5))

ggplot()+
  geom_density(data=long_median_ranked,aes(x=Rank, group=cluster, color=cluster))

factor(clusters$cluster,levels=c(10,4,1,8,9,7,3,2,5,11,6))

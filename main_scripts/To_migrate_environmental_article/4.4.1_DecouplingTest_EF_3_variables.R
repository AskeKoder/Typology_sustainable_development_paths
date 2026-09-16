rm(list = ls())

# ---- Packages ----
library(dplyr)
library(plm)
library(lmtest)
library(sandwich)
library(urca)

# =========================================================
# 0) Read data
# =========================================================
data <- read.csv("TFPdata_V2.csv")%>%
  select(-X)%>%
  filter(TFP != 0)

#data <- read.csv("TFPdata_V3.csv")

clusters <-readRDS("4_RankedClusters.RDS")%>%
  select(Country,iso3=SPI_countrycode,Cluster = Baseline)

#Combine
df <- merge(data, clusters, by = c("iso3", "Country")) %>%
  dplyr::filter(GHG > 0, GDP_PPPcap > 0) %>%  # remove zeros before log
  dplyr::mutate(
    Cluster = factor(Cluster),
    iso3     = factor(iso3),
    country  = iso3,
    year     = Year,
    lnTFP    = log(TFP),
    lnGHG     = log(GHG),
    #lnEF     = log(GHG),
    lnGDP    = log(GDP_PPPcap)
  ) 


 # df <- df %>%
 #   group_by(country) %>%
 #   arrange(year, .by_group = TRUE) %>%
 #   mutate(
 #     lnTFP = cumsum(ln.TFP._WB),        # base year lnTFP = 0
 #     TFP_index = exp(lnTFP)    # optional, level index with base = 1
 #   ) %>%
 #   ungroup()



# #Compute EF score-----------------------------------------
 world_citizen <- df%>%
   filter(year == 2020)%>%
   summarise(across(c(GHG, Biodiversity_Impact, Scarce_Water_Consumption), ~ weighted.mean(., w = Population)))

# #Normalize by world average citizen
 df<- df %>%
   #pers eq. per indicator
   mutate(nGHG = GHG/world_citizen$GHG,
          nBiodiversity_Impact = Biodiversity_Impact / world_citizen$Biodiversity_Impact,
          nScarce_Water_Consumption = Scarce_Water_Consumption / world_citizen$Scarce_Water_Consumption)%>%

   #Compute combined score
   #mutate(pers.eq.mean = rowMeans(select(.,c("nGHG", "nBiodiversity_Impact", "nScarce_Water_Consumption"))))%>%
   mutate(pers.eq.gMean = (nGHG*nBiodiversity_Impact*nScarce_Water_Consumption)^(1/3))%>%

   rename(EF = pers.eq.gMean)%>%
   mutate(lnEF = log(EF))
   #write.csv(df, "TFPdata_EF.csv")




df <- df %>%
  group_by(Cluster, iso3) %>%
  mutate(n_years = n_distinct(year)) %>%       # years per country within each cluster
  group_by(Cluster) %>%
  mutate(max_years = max(n_years)) %>%         # max years per cluster
  ungroup() %>%
  # keep:
  # - all observations from other clusters
  # - only countries with full panel in clusters 6 and 9
  filter(!(Cluster %in% c(6, 9)) | n_years == max_years) %>%
  select(-n_years, -max_years)



# =========================================================
# 0.1) Test for Unit roots
# =========================================================


testUnitRoot <- function(dg, exo = "trend") {
  pdata <- pdata.frame(dg %>% arrange(country, year),
                       index = c("country", "year"))
  vars  <- c("lnEF", "lnGDP", "lnTFP")
  
  get_p <- function(var, test_name) {
    x <- pdata[[var]]
    out <- tryCatch(
      purtest(x,
              data   = pdata,
              index  = c("country","year"),
              test   = test_name,
              pmax   = 4,
              exo    = exo,
              lags   = "AIC"),
      error = function(e) NA_real_
    )
    if (is.na(out)[1]) return(NA_real_)
    summary(out)$statistic$p.value[1]
  }
  
  ips      <- sapply(vars, get_p, test_name = "ips")
  levinlin <- sapply(vars, get_p, test_name = "levinlin")
  madwu    <- sapply(vars, get_p, test_name = "madwu")
  
  round(data.frame(IPS = ips,
                   levinlin = levinlin,
                   madwu = madwu), 5)
}

#Test for stationarity
for (g in sort(unique(df$Cluster))) {
  dg <- df %>% filter(Cluster == g)
  cat("Cluster", g,":\n")
  print(testUnitRoot(dg))
}



# Potential selection: 3,5,6,8,9,10,11


#Test if stationary after differencing
for (g in c(2,4,7,11)){
  if (g == 8) next
  dg <- df %>% filter(Cluster == g)
  if (g == 1) dg <- dg %>% filter(country != "USA")
  dg <- dg %>%
    group_by(country) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(
      across(
        c(lnEF, lnGDP, lnTFP),
        ~ . - dplyr::lag(.)
      )
    ) %>%
    ungroup()
  print(paste("Cluster", g,"differenced:"))
  print(testUnitRoot(dg))
}

# I(1) clusters:  2,3,4,7,11

# Let's affine with Pedronii or KAO test


# =========================================================
#### Kao test (evidence of cointegration @ p<0.05)
# =========================================================

for (g in c(2,4,7,11)){
  dg <- df %>% filter(Cluster == g)
  if (g == 1) dg <- dg %>% filter(country != "USA")
  
  cat("\nCluster", g, "- manual Kao test (null: no cointegration)\n")
  model_pool <- plm(lnEF ~ lnGDP + lnTFP + year, data = dg,
                    index = c("country","year"), model = "within")
  resid_panel <- residuals(model_pool)
  test <- purtest(resid_panel, test = "levinlin", exo = "none", lags = 1)
  print(summary(test))
}

# Final selection : 2,3,4,7,11


# =========================================================
# 1) Long-run (Panel DOLS ≈ FMOLS) and ECT construction
# =========================================================
estimate_longrun_cluster <- function(dg, leads_lags = 1, time_fe = TRUE) {
  pdata <- pdata.frame(dg %>% arrange(country, year), index = c("country","year"))
  
  # Build first differences for DOLS corrections (±1 lead/lag)
  pdata$DlnGDP <- diff(pdata$lnGDP)
  pdata$DlnTFP <- diff(pdata$lnTFP)
  for (k in 1:leads_lags) {
    pdata[[paste0("L",k,"_DlnGDP")]] <- lag(pdata$DlnGDP, k)
    pdata[[paste0("F",k,"_DlnGDP")]] <- lead(pdata$DlnGDP, k)
    pdata[[paste0("L",k,"_DlnTFP")]] <- lag(pdata$DlnTFP, k)
    pdata[[paste0("F",k,"_DlnTFP")]] <- lead(pdata$DlnTFP, k)
  }
  
  rhs <- c("lnGDP","lnTFP","DlnGDP","DlnTFP",
           paste0("L",1:leads_lags,"_DlnGDP"),
           paste0("F",1:leads_lags,"_DlnGDP"),
           paste0("L",1:leads_lags,"_DlnTFP"),
           paste0("F",1:leads_lags,"_DlnTFP"))
  form <- as.formula(paste("lnEF ~", paste(rhs, collapse=" + ")))
  
  fit <- plm(form, data = pdata,
             model = "within",
             #effect = if (time_fe) "twoways" else "individual",
             within=TRUE)
  
  # Robust SEs
  summ <- coeftest(fit, vcov = vcovHC(fit, method="arellano", type="HC1", cluster="group"))
  
  beta_gdp <- coef(fit)["lnGDP"]
  beta_tfp <- coef(fit)["lnTFP"]
  
  # Recover country fixed effects (α_i)
  alphas <- fixef(fit, effect = "individual")
  
  list(beta_gdp = beta_gdp,
       beta_tfp = beta_tfp,
       alphas   = alphas,
       pdata    = pdata)
}

# ======================================================
# 2) Short-run ECM estimation for each cluster
# ======================================================
estimate_ecm_cluster <- function(dg, lr, include_lags = TRUE, time_fe = TRUE) {
  dg <- dg %>% arrange(country, year)
  pdata <- pdata.frame(dg, index = c("country","year"))
  
  # Map α_i (country intercepts)
  alpha_map <- lr$alphas
  pdata$alpha_i <- alpha_map[as.character(pdata$country)]
  
  # Error-correction term (ECT)
  pdata$ECT <- pdata$lnEF - pdata$alpha_i -
    lr$beta_gdp*pdata$lnGDP - lr$beta_tfp*pdata$lnTFP
  pdata$ECT_lag <- lag(pdata$ECT, 1)
  
  # Short-run differences
  pdata$d_lnEF  <- diff(pdata$lnEF)
  pdata$d_lnGDP <- diff(pdata$lnGDP)
  pdata$d_lnTFP <- diff(pdata$lnTFP)
  
  if (include_lags) {
    pdata$d_lnEF_L1  <- lag(pdata$d_lnEF, 1)
    pdata$d_lnGDP_L1 <- lag(pdata$d_lnGDP, 1)
    pdata$d_lnTFP_L1 <- lag(pdata$d_lnTFP, 1)
    form_sr <- d_lnEF ~ ECT_lag + d_lnGDP + d_lnTFP +
      d_lnEF_L1 + d_lnGDP_L1 + d_lnTFP_L1
  } else {
    form_sr <- d_lnEF ~ ECT_lag + d_lnGDP + d_lnTFP
  }
  
  ecm <- plm(form_sr, data = pdata,
             model = "within",
             effect = if (time_fe) "twoways" else "individual")
  
  ecm_sum <- coeftest(ecm, vcov = vcovHC(ecm, method="arellano", type="HC1", cluster="group"))
  
  list(phi       = coef(ecm)["ECT_lag"],
       gamma_gdp = coef(ecm)["d_lnGDP"],
       gamma_tfp = coef(ecm)["d_lnTFP"],
       ecm_fit   = ecm,
       ecm_sum   = ecm_sum)
}

# ======================================================
# 3) Loop over clusters and assemble result tables
# ======================================================
run_all_clusters <- function(df, leads_lags = 1, include_lags = TRUE, time_fe = TRUE) {
  clusters <- sort(unique(df$Cluster))
  out_lr <- list(); out_sr <- list()
  
  for (g in clusters) {
    cat("\nProcessing cluster:", g, "\n")
    dg <- df %>% filter(Cluster == g)
    
    # Long-run FMOLS-like step
    lr <- estimate_longrun_cluster(dg, leads_lags = leads_lags, time_fe = time_fe)
    
    # Short-run ECM
    sr <- estimate_ecm_cluster(dg, lr, include_lags = include_lags, time_fe = time_fe)
    
    lr_row <- tibble(
      Cluster = g,
      `β₁ (GDP elasticity)` = round(lr$beta_gdp, 3),
      `β₂ (TFP elasticity)` = round(lr$beta_tfp, 3)
    )
    sr_row <- tibble(
      Cluster = g,
      `φ (speed of adj.)` = round(sr$phi, 3),
      `γ_GDP (short-run)` = round(sr$gamma_gdp, 3),
      `γ_TFP (short-run)` = round(sr$gamma_tfp, 3)
    )
    
    out_lr[[g]] <- lr_row
    out_sr[[g]] <- sr_row
  }
  
  list(
    long_run  = bind_rows(out_lr),
    short_run = bind_rows(out_sr)
  )
}

# ======================================================
# 4) Run it
# ======================================================
# df must contain: Cluster, country, year, lnTFP, lnGDP, lnEF
results <- run_all_clusters(df, leads_lags = 1, include_lags = TRUE, time_fe = TRUE)
results$long_run
results$short_run


df_summary <- df %>%
  group_by(Cluster) %>%
  summarize(
    count_IDs = n_distinct(iso3),
    sum_X_2020 = sum(Population[year == 2020], na.rm = TRUE)
  )
df_summary

# ---- Packages ----
library(dplyr)
library(plm)
library(lmtest)
library(sandwich)
library(urca)

# =========================================================
# 0) Read data
# =========================================================
data <- read.csv("TFPdata.csv")%>%
  select(-X)

clusters <-readRDS("4_RankedClusters.RDS")%>%
  select(Country,iso3=SPI_countrycode,Cluster = Baseline)

#Combine
df <- merge(data, clusters, by = c("iso3", "Country")) %>%
  dplyr::filter(TFP > 0, GHG > 0, GDP_PPPcap > 0) %>%  # remove zeros before log
  dplyr::mutate(
    Cluster = factor(Cluster),
    iso3     = factor(iso3),
    country  = iso3,
    year     = Year,
    lnTFP    = log(TFP),
    lnEF     = log(GHG),
    lnGDP    = log(GDP_PPPcap)
  ) %>%
  dplyr::select(Cluster, country, year, lnTFP, lnEF, lnGDP)

# =========================================================
# 0.1) Test for Unit roots
# =========================================================
testUnitRoot <- function(dg,exo="trend"){
  pdata<- pdata.frame(dg %>% arrange(country, year), index = c("country","year"))
  vars  <- c("lnEF", "lnGDP", "lnTFP")
  
  ips_test_results <- sapply(vars, function(v) {
    #Set up data for test
    y <- data.frame(split(pdata[,v], pdata$country))%>%
      # Remove rows that are all NA
      filter(if_any(everything(), ~ !is.na(.))) %>%
      # Remove columns that are all NA
      select(where(~ !all(is.na(.))))
    
    tryCatch({
      test <-  purtest(y,
                       #index = c("country","year"),
                       test = "ips",
                       pmax = 4, exo = exo,
                       lags="AIC")},
      error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
    summary(test)$statistic$p.value[1]  # extract p-value
  })
  
  levinlin_test_results <- sapply(vars, function(v) {
    #Set up data for test
    y <- data.frame(split(pdata[,v], pdata$country))%>%
      # Remove rows that are all NA
      filter(if_any(everything(), ~ !is.na(.))) %>%
      # Remove columns that are all NA
      select(where(~ !all(is.na(.))))
    tryCatch({
    test <-  purtest(y,
                     #index = c("country","year"),
                     test = "levinlin",
                     pmax = 4, exo = exo,
                     lags="AIC")},
    error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
    summary(test)$statistic$p.value[1]  # extract p-value
  })
  
  madwu_test_results <- sapply(vars, function(v) {
    #Set up data for test
    y <- data.frame(split(pdata[,v], pdata$country))%>%
      # Remove rows that are all NA
      filter(if_any(everything(), ~ !is.na(.))) %>%
      # Remove columns that are all NA
      select(where(~ !all(is.na(.))))
    
    tryCatch({
      test <-  purtest(y,
                       #index = c("country","year"),
                       test = "madwu",
                       pmax = 4, exo = exo,
                       lags="AIC")},
      error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
    summary(test)$statistic$p.value[1]  # extract p-value
  })
  
  round(data.frame("IPS" = ips_test_results,
             "levinlin" = levinlin_test_results,
             "madwu" = madwu_test_results),5)
}
for (g in sort(unique(df$Cluster))){
  if (g == 8) next
  dg <- df %>% filter(Cluster == g)
  if (g == 1) dg <- dg %>% filter(country != "USA")
  print(paste("Cluster", g,":"))
  print(testUnitRoot(dg))
}

# Potential selection: 2,3,4,5,7,11
# Let's affine with Pedronii or KAO test


# =========================================================
#### Kao test
# =========================================================

for (g in c(2,3,4,5,7,11)) {
  dg <- df %>% filter(Cluster == g)
  if (g == 1) dg <- dg %>% filter(country != "USA")
  
  cat("\nCluster", g, "- manual Kao test (null: no cointegration)\n")
  model_pool <- plm(lnEF ~ lnGDP + lnTFP, data = dg,
                    index = c("country","year"), model = "pooling")
  resid_panel <- residuals(model_pool)
  test <- purtest(resid_panel, test = "levinlin", exo = "none", lags = 1)
  print(summary(test))
}

# Final selection : 2,4,7,9


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
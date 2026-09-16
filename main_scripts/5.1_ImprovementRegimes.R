# =============================================================================
#  5.1_ImprovementRegimes.R  —  unordered "improvement regimes":
#            how fast each country improved, given where it started
#
#  This is the complementary counterpart to the LEVEL clustering of
#  3_ClusterExperimentation_ver2.R.  That pipeline clusters the 21-year
#  trajectories of the indicators themselves and produces an ORDERED ranking
#  (4_RankedClusters.RDS).  Here the pipeline is kept identical - Euclidean
#  distance, Ward.D2, NbClust c-index for k, and the Mirkin-distance ensemble
#  medoid over the 15 imputations - and ONLY the features change:
#
#     PCA on the pooled country-years   ->   dropped
#     21-year trajectory as the object  ->   one improvement rate per indicator
#                                            (OLS slope on year, 2000-2020),
#                                            residualised on the year-2000 level
#                                            with a cubic polynomial
#
#  Residualising on the starting value removes the ceiling effect, so the
#  features read "faster or slower than expected given where you started".
#  The resulting groups are therefore UNORDERED regimes, not a ranking.
#
#  PRODUCES
#     the per-imputation k and the ensemble (medoid) k
#     the paper's table: regimes x 11 indicators, mean improvement (points/yr)
#     regime size, mean development level and member countries
#     between-regime share of variance in the improvement rates, against a
#        size-matched random-partition null (200 draws)
#     regime stability: pairwise adjusted Rand index across the 15 imputations
#     adjusted Rand index against the published DLSFew_coverage partition and
#        the R2 of development level explained by regime
#     Revision_output/6_ImprovementRegimes.RDS  and  .csv   (regime per country)
#     Revision_output/6_RegimeTable_rawRates.csv             (the printed table)
#     Revision_output/6_RegimeDiagnostics.txt                (everything else printed)
#  Revision_output/ sits at the REPOSITORY ROOT, i.e. one level above
#  main_scripts/ when the script is run from there.
#
#  Run from main_scripts/ :        source("5.1_ImprovementRegimes.R")
#  or from the repository root:    source("main_scripts/5.1_ImprovementRegimes.R")
#  Requires: NbClust.  Base R otherwise; runtime is roughly one minute.
# =============================================================================
rm(list = ls())

suppressMessages(library(NbClust))   # for the c-index, as in 3_ClusterExperimentation_ver2.R

set.seed(20240517)                   # only the random-partition null is stochastic

# ---- 0. Paths (script runs from main_scripts/ or from the repository root) ---
here <- if (file.exists("2_ImputedData.csv")) "" else "main_scripts/"
if (!file.exists(paste0(here, "2_ImputedData.csv")))
  stop("Cannot find 2_ImputedData.csv - run from main_scripts/ or from the repository root.")
#Derived from `here`, not from whether the directory happens to exist: testing
#for the directory would write to the PARENT of the repository on a clean
#checkout sourced from the repository root.
outdir <- if (here == "") "../ImprovementRegimes" else "Revision_output"
dir.create(outdir, showWarnings = FALSE)

# ---- 1. Data and indicators -------------------------------------------------
data <- read.csv(paste0(here, "2_ImputedData.csv"), encoding = "UTF-8")
data <- data[, !(names(data) %in% c("X", ".id"))]

#The 11 indicators of the reference specification (Few_indicators_closest_DLS_coverage)
indicators <- c("Nutritional_deficiencies",
                "No_access_to_a_handwashing_facility",
                "Improved_sanitation",
                "Improved_water_source",
                "Prevalence_of_cooking_with_coalperbiomass",
                "Access_to_electricity",
                "Mean_years_of_schooling",
                "Internet_users",
                "Mobile_and_landline_telephone_subscriptions",
                "Universal_health_coverage",
                "Share_Slums")

#Short labels for printing
short <- c("Nutrition", "Handwashing", "Sanitation", "Water", "CleanFuels",
           "Electricity", "Schooling", "Internet", "Telephone", "UHC", "Slums")

m          <- length(unique(data$.imp))          # 15 imputations
countries  <- unique(data$Country)
isocodes   <- data$SPI_countrycode[match(countries, data$Country)]   # ASCII keys
nCountries <- length(countries)
years      <- sort(unique(data$SPI_year))
nYears     <- length(years)

cat("\n=========================================================\n")
cat("  Improvement regimes - starting-point-adjusted rates\n")
cat("=========================================================\n\n")
cat(sprintf("Countries: %d | years: %d-%d | imputations: %d | indicators: %d\n",
            nCountries, min(years), max(years), m, length(indicators)))

# ---- 2. Orientation check ---------------------------------------------------
#All 11 series are stored higher = better EXCEPT Share_Slums, which is higher =
#worse and is negated below.  Two names are misleading and are NOT negated:
#"No_access_to_a_handwashing_facility" is really the share WITH a basic
#handwashing facility, and "Prevalence_of_cooking_with_coalperbiomass" is really
#access to CLEAN fuels.  We verify this empirically rather than trusting the
#names: every series must move WITH access to electricity, except Share_Slums.
check_set <- data[data$.imp == 1, ]
orient    <- sapply(indicators, function(v) cor(check_set[, v],
                                                check_set[, "Access_to_electricity"]))
expected  <- ifelse(indicators == "Share_Slums", -1, 1)

cat("\n--- Orientation check (correlation with Access_to_electricity, imputation 1) ---\n")
print(data.frame(Indicator = short, r = round(orient, 3),
                 Expected  = ifelse(expected > 0, "higher = better", "higher = worse"),
                 row.names = NULL), right = FALSE)
if (any(sign(orient) != expected))
  stop("Orientation check FAILED for: ",
       paste(indicators[sign(orient) != expected], collapse = ", "),
       ". The sign convention assumed by this script does not hold in the data.")
cat("Orientation check passed: the two negatively-named series (handwashing,\n")
cat("cooking fuel) are stored higher = better; only Share_Slums is negated.\n")

#Apply the single sign flip
data$Share_Slums <- -data$Share_Slums

# ---- 3. Functions -----------------------------------------------------------
#Improvement rates and starting levels for one imputed set
rateFeatures <- function(set){
  #set is one imputation, balanced country x year, ordered country then year
  set <- set[order(match(set$Country, countries), set$SPI_year), ]
  if (nrow(set) != nCountries * nYears)
    stop("The imputed set is not a balanced country x year panel.")
  cy  <- years - mean(years)                       # centred time, common to all countries
  SL  <- matrix(0, nrow = nCountries, ncol = length(indicators))   # slope, points/year
  ST  <- matrix(0, nrow = nCountries, ncol = length(indicators))   # year-2000 level
  LV  <- matrix(0, nrow = nCountries, ncol = length(indicators))   # mean level 2000-2020
  rownames(SL) <- rownames(ST) <- rownames(LV) <- countries
  colnames(SL) <- colnames(ST) <- colnames(LV) <- indicators
  for (j in seq_along(indicators)){
    #Column per country, one row per year: the OLS slope is then a simple ratio
    M       <- matrix(set[, indicators[j]], nrow = nYears)
    SL[, j] <- as.vector(crossprod(cy, M)) / sum(cy^2)
    ST[, j] <- M[1, ]
    LV[, j] <- colMeans(M)
  }
  list(SL = SL, ST = ST, LV = LV)
}

#Adjusted Rand index (base R)
adjRand <- function(a, b){
  n   <- length(a)
  tab <- table(a, b)
  s   <- sum(choose(tab, 2))
  sa  <- sum(choose(rowSums(tab), 2))
  sb  <- sum(choose(colSums(tab), 2))
  e   <- sa * sb / choose(n, 2)
  (s - e) / (0.5 * (sa + sb) - e)
}

#Mirkin distance (base R): pairs together in one partition but not the other.
#This is the N10 + N01 of computePairCoefficients() used in
#3_ClusterExperimentation_ver2.R; the medoid rule is identical.
mirkinDist <- function(a, b){
  tab <- table(a, b)
  N11 <- sum(choose(tab, 2))
  N10 <- sum(choose(rowSums(tab), 2)) - N11
  N01 <- sum(choose(colSums(tab), 2)) - N11
  N10 + N01
}

#Between-group share of variance, pooled over standardised features
betweenShare <- function(Z, g){
  tot <- sum(sweep(Z, 2, colMeans(Z))^2)
  bet <- 0
  for (k in unique(g)){
    idx <- which(g == k)
    bet <- bet + length(idx) * sum((colMeans(Z[idx, , drop = FALSE]) - colMeans(Z))^2)
  }
  bet / tot
}

# ---- 4. Cluster every imputation --------------------------------------------
clusterings <- matrix(0, nrow = nCountries, ncol = m)
rownames(clusterings) <- countries
nclust <- rep(0, m)
featSL <- vector("list", m)
featST <- vector("list", m)
featZ  <- vector("list", m)   # the level-adjusted, standardised features clustered
featLV <- vector("list", m)   # mean level 2000-2020, for the development scalar

cat("\n--- Clustering the 15 imputed sets -------------------------------------\n")
for (set in 1:m){
  complete_set <- data[data$.imp == set, ]
  f  <- rateFeatures(complete_set)
  SL <- f$SL; ST <- f$ST
  featSL[[set]] <- SL; featST[[set]] <- ST; featLV[[set]] <- f$LV

  #Residualise each improvement rate on its own starting value (cubic)
  RS <- SL
  for (j in seq_along(indicators))
    RS[, j] <- resid(lm(SL[, j] ~ poly(ST[, j], 3)))

  #Standardise, Euclidean distance, and the authors' k-selection rule
  Z <- scale(RS)
  featZ[[set]] <- Z
  D <- dist(Z)
  #capture.output only silences NbClust's "only frey, mcclain, cindex ..." notice
  invisible(capture.output(
    nb <- suppressWarnings(NbClust(data = NULL, diss = D, distance = NULL,
                                   min.nc = 2, max.nc = 40,   # as in script 3
                                   index = "cindex",
                                   method = "ward.D2"))))
  nclust[set]       <- nb$Best.nc[1]
  clusterings[, set] <- cutree(hclust(D, method = "ward.D2"), k = nb$Best.nc[1])
  cat(sprintf("Set %2d finished with %2d clusters.\n", set, nclust[set]))
}

# ---- 5. Ensemble: Mirkin-distance medoid ------------------------------------
partitionQuality <- sapply(1:m, function(p)
  sum(sapply(1:m, function(x) mirkinDist(clusterings[, p], clusterings[, x]))))
best <- which.min(partitionQuality)

cat("\n--- Ensemble (medoid of the 15 partitions under the Mirkin distance) ---\n")
print(data.frame(Imputation = 1:m, k = nclust,
                 MirkinTotal = partitionQuality, row.names = NULL), right = FALSE)
cat(sprintf("\nMedoid imputation: %d   ->   k = %d regimes\n", best, nclust[best]))
cat(sprintf("Number of clusters over the 15 sets: median %.0f, range %d-%d\n",
            median(nclust), min(nclust), max(nclust)))

k       <- nclust[best]
regime0 <- clusterings[, best]
#The PARTITION comes from the medoid imputation.  Every QUANTITY reported below
#(improvement rates, starting levels, development level) is instead averaged
#over all 15 imputed sets, as elsewhere in this repository, so that the numbers
#in the paper do not depend on which single imputation supplied the medoid.
SL      <- Reduce(`+`, featSL) / m
ST      <- Reduce(`+`, featST) / m

# ---- 6. Order the regimes by development level ------------------------------
#Development level = mean of the standardised year-2000 values (higher = better
#provisioning at the start).  Regimes are UNORDERED; this only fixes the labels
#so that R1 is the highest-provisioning group.
devlevel <- rowMeans(scale(ST))
ord      <- order(tapply(devlevel, regime0, mean), decreasing = TRUE)
regime   <- match(regime0, ord)
names(regime) <- countries

# ---- 7. The table: regimes x indicators, mean improvement (points/year) ------
meanRate <- t(sapply(1:k, function(g) colMeans(SL[regime == g, , drop = FALSE])))
rownames(meanRate) <- paste0("R", 1:k)
colnames(meanRate) <- short

tab <- data.frame(Indicator = short,
                  round(t(meanRate), 2),
                  Overall   = round(colMeans(SL), 2),
                  check.names = FALSE, row.names = NULL)

cat("\n===================== IMPROVEMENT REGIMES (points/year) =====================\n")
cat("Mean annual improvement per indicator, averaged over the 15 imputed sets.\n")
cat("Regimes in columns,\n")
cat("ordered by development level in 2000 (R1 = highest). Slums is sign-flipped, so\n")
cat("a positive number is always an improvement.\n\n")
print(tab, row.names = FALSE, right = FALSE)
write.csv(tab, file.path(outdir, "6_RegimeTable_rawRates.csv"), row.names = FALSE)

# ---- 8. Regime profiles -----------------------------------------------------
cat("\n--- Regime size, development level and members ------------------------------\n")
for (g in 1:k){
  mem <- sort(countries[regime == g])
  cat(sprintf("\nR%d  n = %d   mean development level (std. 2000 value) = %+.2f\n",
              g, length(mem), mean(devlevel[regime == g])))
  cat("   ", paste(strwrap(paste(mem, collapse = ", "), width = 72), collapse = "\n    "), "\n")
}

# ---- 9. How much of the improvement variation is between regimes? -----------
#Measured on the features actually clustered: the level-adjusted, standardised
#improvement rates.  Computed inside every imputation and reported as the mean
#and range over the 15, so the figure is not read off the single imputed set
#that happened to supply the medoid partition.
obsVec  <- sapply(featZ, betweenShare, g = regime)
nullVec <- unlist(lapply(featZ, function(Z)
             replicate(200, betweenShare(Z, sample(regime)))))  # identical group sizes

cat("\n--- Between-regime share of variance in the adjusted improvement rates -------\n")
cat(sprintf("Observed                       : mean %.3f over the 15 sets (range %.3f-%.3f)\n",
            mean(obsVec), min(obsVec), max(obsVec)))
cat(sprintf("Random partitions, same sizes  : mean %.3f, 97.5th pct %.3f (%d draws: 200 per imputed set)\n",
            mean(nullVec), quantile(nullVec, 0.975), length(nullVec)))
cat(sprintf("Draws at or above the observed  : %d / %d\n",
            sum(nullVec >= mean(obsVec)), length(nullVec)))
#For reference only: the same statistic, computed the same way (inside each
#imputation, then averaged), on the RAW rates before the starting-level
#adjustment.  The gap between the two is what residualising removes.
cat(sprintf("[reference] same share on the raw, unadjusted rates : %.3f\n",
            mean(sapply(featSL, function(S) betweenShare(scale(S), regime)))))

# ---- 10. Stability across the imputations -----------------------------------
pairs <- combn(m, 2)
ari   <- apply(pairs, 2, function(p) adjRand(clusterings[, p[1]], clusterings[, p[2]]))
cat("\n--- Regime stability across the 15 imputations (adjusted Rand index) --------\n")
cat(sprintf("mean %.3f | min %.3f | max %.3f  (%d pairs)\n",
            mean(ari), min(ari), max(ari), ncol(pairs)))

# ---- 11. Relation to the published LEVEL ranking ----------------------------
rds <- if (file.exists(paste0(here, "4_RankedClusters.RDS")))
         paste0(here, "4_RankedClusters.RDS") else "4_RankedClusters.RDS"
ranked <- suppressWarnings(readRDS(rds))
mm     <- match(isocodes, ranked$SPI_countrycode)   # ISO codes, not accented names
level  <- ranked$DLSFew_coverage[mm]

cat("\n--- Improvement regimes vs the published DLSFew_coverage ranking ------------\n")
cat(sprintf("Countries matched              : %d of %d\n", sum(!is.na(mm)), nCountries))
cat(sprintf("Adjusted Rand index            : %.3f\n",
            adjRand(regime[!is.na(mm)], level[!is.na(mm)])))
#The development SCALAR reported in the SI: the mean of the eleven standardised
#indicators averaged over 2000-2020 (Share_Slums already negated above), pooled
#over the imputations.  This is a country's overall provisioning level, and is
#the variable both R2 figures in the paper are computed on.  It standardises
#within each imputation and then averages; `devlevel` above averages and then
#standardises.  The two conventions differ by less than 0.005 in the R2 and the
#choice is immaterial, but they are not interchangeable, so keep them labelled.
scalarMat <- Reduce(`+`, lapply(featLV, function(L) scale(L))) / m
scalar    <- rowMeans(scalarMat)
cat(sprintf("R2 of development scalar by regime  : %.3f\n",
            summary(lm(scalar ~ factor(regime)))$r.squared))
cat(sprintf("R2 of development scalar by ranking : %.3f  (same variable, for comparison)\n",
            summary(lm(scalar[!is.na(mm)] ~ factor(level[!is.na(mm)])))$r.squared))
cat(sprintf("  [year-2000 level instead of the 2000-2020 mean: %.3f and %.3f]\n",
            summary(lm(devlevel ~ factor(regime)))$r.squared,
            summary(lm(devlevel[!is.na(mm)] ~ factor(level[!is.na(mm)])))$r.squared))
cat("A low adjusted Rand index is the point: the regimes are not a relabelling of\n")
cat("the level ranking, they are a second, complementary partition.\n")
cat("\nCross-tabulation (rows = improvement regime, columns = level cluster):\n")
print(table(Regime = regime[!is.na(mm)], Level = level[!is.na(mm)]))

# ---- 12. Save ---------------------------------------------------------------
out <- data.frame(Country          = countries,
                  SPI_countrycode  = isocodes,
                  ImprovementRegime = as.integer(regime),
                  DevelopmentLevel = round(devlevel, 4),
                  row.names = NULL, stringsAsFactors = FALSE)
saveRDS(out, file.path(outdir, "6_ImprovementRegimes.RDS"))
write.csv(out, file.path(outdir, "6_ImprovementRegimes.csv"), row.names = FALSE)

#Everything printed above that the paper or the SI uses, as one text file, so
#that an interactive run leaves the same record as 00_RunAll.R.
diagf <- file.path(outdir, "6_RegimeDiagnostics.txt")
con <- file(diagf, open = "wt"); sink(con)
cat("Improvement regimes - diagnostics, written", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")
cat("k per imputed set and Mirkin total (medoid = smallest total):\n")
print(data.frame(Imputation = 1:m, k = nclust, MirkinTotal = partitionQuality, row.names = NULL), right = FALSE)
cat(sprintf("\nMedoid imputation %d -> k = %d regimes; k over the 15 sets: median %.0f, range %d-%d\n",
            best, k, median(nclust), min(nclust), max(nclust)))
cat(sprintf("Regime sizes: %s\n", paste(sprintf("R%d = %d", 1:k, tabulate(regime, k)), collapse = ", ")))
cat(sprintf("Between-regime share of variance, adjusted rates: %.3f (range %.3f-%.3f over the 15 sets);\n",
            mean(obsVec), min(obsVec), max(obsVec)))
cat(sprintf("  size-matched random partitions: mean %.3f, 97.5th pct %.3f, draws >= observed %d / %d\n",
            mean(nullVec), quantile(nullVec, 0.975), sum(nullVec >= mean(obsVec)), length(nullVec)))
cat(sprintf("Stability across imputations, adjusted Rand index: mean %.3f, min %.3f, max %.3f (%d pairs)\n",
            mean(ari), min(ari), max(ari), ncol(pairs)))
cat(sprintf("Adjusted Rand index vs the DLSFew_coverage ranking: %.3f\n", adjRand(regime[!is.na(mm)], level[!is.na(mm)])))
cat(sprintf("R2 of the development scalar: by regime %.3f, by ranking %.3f\n",
            summary(lm(scalar ~ factor(regime)))$r.squared,
            summary(lm(scalar[!is.na(mm)] ~ factor(level[!is.na(mm)])))$r.squared))
cat("\nCross-tabulation (rows = improvement regime, columns = level cluster):\n")
print(table(Regime = regime[!is.na(mm)], Level = level[!is.na(mm)]))
cat("\nMean annual improvement per indicator (raw points/year), regimes in columns:\n")
print(tab, row.names = FALSE, right = FALSE)
sink(); close(con)

cat("\nWritten to", normalizePath(outdir), ":\n")
cat("   6_ImprovementRegimes.RDS / .csv   regime and development level, 173 countries\n")
cat("   6_RegimeTable_rawRates.csv        the table above\n")
cat("   6_RegimeDiagnostics.txt           k per set, medoid, variance share, stability, cross-tab\n")
cat("\nNote: k is chosen by NbClust's c-index inside each imputation and the final\n")
cat("partition is the Mirkin medoid of the 15; nothing here fixes k by hand.\n")

# =============================================================================
#  Table 5 | Ordinal (proportional-odds) models of cluster membership
#            on historical covariates
#
#  Produces the two columns reported in the manuscript:
#     (1) settler mortality + primary enrolment 1900
#     (2) the same, plus coloniser identity
#  together with the type II likelihood-ratio tests and fit statistics.
#
#  Coloniser is coded as the power ruling AT INDEPENDENCE (see 4.3.6 for the
#  comparison with the longest-lasting rule and why we prefer this one).
#
#  Run from the repository root:  source("main_scripts/4.3.5_Table5.R")
#  Requires: MASS
# =============================================================================
suppressMessages(library(MASS))

# ---- 1. Estimation frame ----------------------------------------------------
#f  <- "Data_Colonial/"
#sm <- read.csv(paste0(f,"SettlerMortality.csv"), sep=";", dec=",")

ROOT <- if (dir.exists("Data_Colonial")) "." else ".."
f    <- paste0(ROOT, "/Data_Colonial/")
folder <- f
sm   <- read.csv(paste0(f, "SettlerMortality.csv"), sep=";", dec=",")

sm <- sm[!duplicated(sm$ccodealp), c("ccodealp","ajr_settmort")]   # dedupe on ISO only
names(sm) <- c("iso3","settmort")

cr <- read.csv(paste0(f,"Colonies.csv")); names(cr)[4] <- "power"
pw <- c("Belgium","France","Germany","Italy","Netherlands","Portugal","Spain","United Kingdom")
sub <- cr[cr$power %in% pw, ]
ly  <- aggregate(list(y=sub$Year), by=sub[,c("Code","power")], FUN=max)   # LAST year of rule
ly  <- ly[order(-ly$y), ]
col <- ly[!duplicated(ly$Code), c("Code","power")]; names(col) <- c("iso3","colonizer")

pe <- read.csv(paste0(f,"Prienr.csv")); names(pe)[4] <- "prienr1900"
pe <- pe[pe$Year==1900 & !pe$Code %in% c("","OWID_WRL") & !is.na(pe$Code),
         c("Code","prienr1900")]; names(pe)[1] <- "iso3"

rds <- if (file.exists("main_scripts/4_RankedClusters.RDS"))
         "main_scripts/4_RankedClusters.RDS" else "4_RankedClusters.RDS"
cl  <- readRDS(rds); names(cl)[names(cl)=="SPI_countrycode"] <- "iso3"

full <- merge(merge(sm, col, by="iso3", all=TRUE), pe, by="iso3", all=TRUE)
tb <- table(full$colonizer)
full$colonizer[full$colonizer %in% names(tb)[tb < 10]] <- "Other"   # pool on full universe
d <- merge(full, cl[,c("iso3","DLSFew_coverage")], by="iso3")
names(d)[ncol(d)] <- "cluster"
d <- d[complete.cases(d[,c("settmort","prienr1900","cluster","colonizer")]), ]
d$colonizer <- relevel(factor(d$colonizer), "United Kingdom")
d$ord <- factor(d$cluster, ordered = TRUE)

cat("Estimation sample: N =", nrow(d), " clusters:", nlevels(d$ord), "\n")
print(table(d$colonizer))

# ---- 2. The two specifications ----------------------------------------------
m1 <- polr(ord ~ settmort + prienr1900,             data=d, Hess=TRUE)   # column (1)
m2 <- polr(ord ~ settmort + prienr1900 + colonizer, data=d, Hess=TRUE)   # column (2)
m0 <- polr(ord ~ 1, data=d, Hess=TRUE)

star <- function(p) if (is.na(p)) "" else
  if (p<0.001) "***" else if (p<0.01) "**" else if (p<0.05) "*" else if (p<0.1) "." else ""
cell <- function(m, term){ ct <- coef(summary(m))
  if (!term %in% rownames(ct)) return("")
  p <- 2*pnorm(abs(ct[term,3]), lower.tail=FALSE)
  sprintf("%.3f (%.3f)%s", ct[term,1], ct[term,2], star(p)) }

rows <- c("settmort","prienr1900","colonizerFrance","colonizerOther","colonizerSpain")
labs <- c("Log settler mortality","Primary enrolment 1900 (per pp)",
          "  France","  Other","  Spain")
tab <- data.frame(Term = labs,
                  `(1)` = sapply(rows, function(r) cell(m1,r)),
                  `(2)` = sapply(rows, function(r) cell(m2,r)),
                  check.names = FALSE, stringsAsFactors = FALSE)

cat("\n================================ TABLE 5 ================================\n")
cat("Coefficient (standard error). POSITIVE = higher odds of a LOWER-ranked\n")
cat("cluster (worse provisioning). Coloniser reference category: United Kingdom.\n\n")
print(tab, row.names=FALSE, right=FALSE)

# ---- 3. Type II likelihood-ratio tests --------------------------------------
LR <- function(full, red, df){ v <- deviance(red) - deviance(full)
  p <- pchisq(v, df, lower.tail=FALSE)
  sprintf("chi2 = %6.2f (%d) %-3s p = %.4g", v, df, star(p), p) }

cat("\n--- Type II likelihood-ratio tests (the inference reported) ---\n")
cat("Column (1)\n")
cat("   settler mortality :", LR(m1, polr(ord~prienr1900,          data=d, Hess=TRUE), 1), "\n")
cat("   enrolment 1900    :", LR(m1, polr(ord~settmort,            data=d, Hess=TRUE), 1), "\n")
cat("Column (2)\n")
cat("   settler mortality :", LR(m2, polr(ord~prienr1900+colonizer,data=d, Hess=TRUE), 1), "\n")
cat("   enrolment 1900    :", LR(m2, polr(ord~settmort+colonizer,  data=d, Hess=TRUE), 1), "\n")
cat("   coloniser identity:", LR(m2, m1, 3), "\n")

# ---- 4. Fit statistics ------------------------------------------------------
cat("\n--- Fit statistics ---\n")
for (nm in c("m1","m2")){ m <- get(nm)
  cat(sprintf("   %s : N = %d | parameters = %2d | AIC = %6.1f | McFadden R2 = %.3f\n",
      ifelse(nm=="m1","(1)","(2)"), nobs(m),
      length(coef(m))+length(m$zeta), AIC(m), 1-deviance(m)/deviance(m0))) }

cat("\n--- Thresholds (cut-points between adjacent clusters), column (2) ---\n")
print(round(m2$zeta, 2))
cat("\n--- Odds ratios, column (2) ---\n")
print(round(exp(coef(m2)), 3))
cat("   enrolment per 10 percentage points: OR =",
    round(exp(10*coef(m2)["prienr1900"]), 3), "\n")

cat("\nNote: significance marks on individual coefficients are Wald-based and are\n")
cat("descriptive only. For coloniser identity, a four-level factor, the joint\n")
cat("likelihood-ratio test above is the evidence, not the per-row marks.\n")

# =============================================================================
#  Proportional-odds (parallel regression) diagnostic
#  for the historical-institutional models
#
#  The ordinal specification assumes each covariate has the same effect at every
#  cut-point of the cluster ranking. This script tests that at four resolutions
#  and shows whether the substantive conclusions depend on it.
#
#  Run from the repository root:
#      source("main_scripts/4.3.7_ProportionalOddsCheck.R")
#
#  Requires: MASS, brant        install.packages("brant")
# =============================================================================
suppressMessages({library(MASS); library(brant)})

# ---- 1. Estimation frame (at-independence coloniser coding, N = 55) ---------
ROOT <- if (dir.exists("Data_Colonial")) "." else ".."
f    <- paste0(ROOT, "/Data_Colonial/")
folder <- f
sm   <- read.csv(paste0(f, "SettlerMortality.csv"), sep=";", dec=",")



sm <- sm[!duplicated(sm$ccodealp), c("ccodealp","ajr_settmort")]
names(sm) <- c("iso3","settmort")

cr <- read.csv(paste0(f,"Colonies.csv")); names(cr)[4] <- "power"
pw <- c("Belgium","France","Germany","Italy","Netherlands","Portugal","Spain","United Kingdom")
sub <- cr[cr$power %in% pw, ]
ly  <- aggregate(list(y=sub$Year), by=sub[,c("Code","power")], FUN=max)
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
full$colonizer[full$colonizer %in% names(tb)[tb < 10]] <- "Other"
d <- merge(full, cl[,c("iso3","DLSFew_coverage")], by="iso3")
names(d)[ncol(d)] <- "cluster"
d <- d[complete.cases(d[,c("settmort","prienr1900","cluster","colonizer")]), ]
d$colonizer <- relevel(factor(d$colonizer), "United Kingdom")
cat("Estimation sample: N =", nrow(d), "countries\n")

# ---- 2. Four resolutions of the ordered outcome -----------------------------
d$ord   <- factor(d$cluster, ordered = TRUE)
d$band5 <- cut(d$cluster, c(0,3,5,7,9,11), labels=c("1-3","4-5","6-7","8-9","10-11"), ordered_result=TRUE)
d$band4 <- cut(d$cluster, c(0,4,7,9,11),   labels=c("1-4","5-7","8-9","10-11"),      ordered_result=TRUE)
d$band3 <- cut(d$cluster, c(0,5,8,11),     labels=c("1-5","6-8","9-11"),             ordered_result=TRUE)

LR <- function(a,b,df){ v <- deviance(b) - deviance(a)
  sprintf("chi2=%5.2f df=%d p=%.4f", v, df, pchisq(v, df, lower.tail=FALSE)) }

run <- function(yname, lab){
  d$Y <- d[[yname]]
  m  <- polr(Y ~ settmort + prienr1900 + colonizer, data=d, Hess=TRUE)
  m0 <- polr(Y ~ settmort + prienr1900,             data=d, Hess=TRUE)
  ms <- polr(Y ~ prienr1900 + colonizer,            data=d, Hess=TRUE)
  mp <- polr(Y ~ settmort + colonizer,              data=d, Hess=TRUE)
  cat(sprintf("\n### %-32s %2d categories, %2d parameters\n",
              lab, nlevels(d$Y), length(coef(m)) + length(m$zeta)))
  cat("   n per category    :", paste(table(d$Y), collapse=" / "), "\n")
  cat("   settler mortality :", LR(m, ms, 1), "\n")
  cat("   enrolment 1900    :", LR(m, mp, 1), "\n")
  cat("   coloniser (3 df)  :", LR(m, m0, 3), "\n")
  cat("   coefficients      :",
      paste(sprintf("%s=%.2f", names(coef(m)), coef(m)), collapse="  "), "\n")
  cat(sprintf("   AIC               : %.1f with coloniser, %.1f without\n", AIC(m), AIC(m0)))
  b <- suppressWarnings(try(brant(m), silent=TRUE))
  if (inherits(b, "try-error"))
    cat("   BRANT             : CANNOT BE COMPUTED -",
        sub(":.*", "", attr(b,"condition")$message), "\n")
  else
    cat(sprintf("   BRANT omnibus     : X2 = %.2f, df = %d, p = %.3f  -> %s\n",
        b[1,1], b[1,2], b[1,3],
        ifelse(b[1,3] < 0.05, "REJECTS proportional odds", "assumption holds")))
  invisible(m)
}

cat("\n============ PROPORTIONAL-ODDS DIAGNOSTIC ACROSS RESOLUTIONS ============\n")
cat("With 55 countries over 10 clusters the Brant test must fit nine separate\n")
cat("binary logits, several with only 2-6 observations on one side. Those\n")
cat("component fits are separated, so the test is unreliable or fails outright\n")
cat("at full resolution. Collapsing the ranking into ordered bands gives the\n")
cat("test enough data per cut while preserving the ordering.\n")

for (v in c("ord","band5","band4","band3"))
  run(v, switch(v, ord = "10 clusters (full resolution)",
                   band5 = "5 bands", band4 = "4 bands", band3 = "3 bands"))

cat("\n=========================== HOW TO READ THIS ===========================\n")
cat("If the substantive tests stay significant and the coefficients stay stable\n")
cat("across resolutions, while Brant passes wherever it can actually be\n")
cat("computed, the full-resolution rejection is a sparsity artefact rather than\n")
cat("a real violation. Report the 10-cluster model as the main specification and\n")
cat("the 4-band model as the proportional-odds robustness check.\n")

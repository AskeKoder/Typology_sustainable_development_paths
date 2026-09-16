# =============================================================================
#  Alternative historical covariates from Acemoglu, Gallego & Robinson (2014)
#
#  Tests every plausible AGR variable against the reference specification.
#  Produces SI Table S13.
#
#  Reads Data_Colonial/agr_covariates.csv, a plain-text extract of the AGR
#  replication file xcountry_data.dta (supplied alongside this script — copy it
#  into Data_Colonial/ before running). No Stata reader is required.
#
#  Run from the repository root:  source("main_scripts/4.3.8_AGRVariables.R")
#  Requires: MASS
# =============================================================================
suppressMessages(library(MASS))


ROOT <- if (dir.exists("Data_Colonial")) "." else ".."
f    <- paste0(ROOT, "/Data_Colonial/")
folder <- f
sm   <- read.csv(paste0(f, "SettlerMortality.csv"), sep=";", dec=",")


sm <- sm[!duplicated(sm$ccodealp), c("ccodealp","ajr_settmort")]; names(sm)<-c("iso3","settmort")
cr <- read.csv(paste0(f,"Colonies.csv")); names(cr)[4]<-"power"
pw <- c("Belgium","France","Germany","Italy","Netherlands","Portugal","Spain","United Kingdom")
sub <- cr[cr$power %in% pw,]
ly <- aggregate(list(y=sub$Year), by=sub[,c("Code","power")], FUN=max); ly<-ly[order(-ly$y),]
col <- ly[!duplicated(ly$Code), c("Code","power")]; names(col)<-c("iso3","colonizer")
pe <- read.csv(paste0(f,"Prienr.csv")); names(pe)[4]<-"prienr1900"
pe <- pe[pe$Year==1900 & !pe$Code %in% c("","OWID_WRL") & !is.na(pe$Code), c("Code","prienr1900")]
names(pe)[1] <- "iso3"

rds <- paste0(ROOT, "/main_scripts/4_RankedClusters.RDS")

cl <- readRDS(rds); names(cl)[names(cl)=="SPI_countrycode"] <- "iso3"
fu <- merge(merge(sm,col,by="iso3",all=TRUE), pe, by="iso3", all=TRUE)
tb <- table(fu$colonizer); fu$colonizer[fu$colonizer %in% names(tb)[tb<10]] <- "Other"
d <- merge(fu, cl[,c("iso3","DLSFew_coverage")], by="iso3"); names(d)[ncol(d)] <- "cluster"
d <- d[complete.cases(d[,c("settmort","prienr1900","cluster","colonizer")]),]
d$colonizer <- relevel(factor(d$colonizer),"United Kingdom")

agrfile <- paste0(f,"agr_covariates.csv")
if (!file.exists(agrfile))
  stop("Data_Colonial/agr_covariates.csv not found. Copy it from the revision folder (4_Code/).")
agr  <- read.csv(agrfile)
keep <- c("iso3","lcapped","protmiss","malfal94","lpd1500s","prienr1870","cath1900",
          "prot1900","musl1900","Yrsmis60","lat_abst")
d <- merge(d, agr[,keep], by="iso3", all.x=TRUE)
d$ord <- factor(d$cluster, ordered=TRUE)
cat("N =", nrow(d), "\n\n")

O <- function(fm,x) polr(fm, data=x, Hess=TRUE)
P <- function(a,b,df) pchisq(deviance(b)-deviance(a), df, lower.tail=FALSE)

cat("=== ADDED to settmort + prienr1900 + colonizer ===\n")
cat(sprintf("%-14s %4s %9s %10s %8s\n","variable","N","coef","LR p","dAIC"))
for (v in c("malfal94","lat_abst","musl1900","cath1900","protmiss","lpd1500s",
            "prot1900","Yrsmis60","prienr1870")) {
  x <- d[!is.na(d[[v]]),]; x$ord <- factor(x$cluster, ordered=TRUE); x$V <- x[[v]]
  if (nrow(x) < 25 || nlevels(x$ord) < 4) { cat(sprintf("%-14s %4d  too few observations\n", v, nrow(x))); next }
  a <- try(O(ord~settmort+prienr1900+colonizer+V, x), TRUE)
  b <- try(O(ord~settmort+prienr1900+colonizer,   x), TRUE)
  if (inherits(a,"try-error")) { cat(sprintf("%-14s %4d  did not converge\n", v, nrow(x))); next }
  cat(sprintf("%-14s %4d %9.3f %10.4f %8.1f\n", v, nrow(x), coef(a)["V"], P(a,b,1), AIC(a)-AIC(b)))
}

cat("\n=== REPLACING settler mortality ===\n")
for (v in c("settmort","lcapped","malfal94","lpd1500s")) {
  x <- d[!is.na(d[[v]]),]; x$ord <- factor(x$cluster, ordered=TRUE); x$V <- x[[v]]
  a <- O(ord~V+prienr1900+colonizer, x); b <- O(ord~prienr1900+colonizer, x)
  cat(sprintf("  %-10s N=%2d  coef=%7.3f  LR p=%.4g  AIC=%6.1f  McFadden=%.3f\n",
      v, nrow(x), coef(a)["V"], P(a,b,1), AIC(a), 1-deviance(a)/deviance(O(ord~1,x))))
}
cat("\nNote: malfal94 is measured in 1994 and partly reflects successful eradication,\n")
cat("so it is a robustness check on the disease-environment channel, not a substitute\n")
cat("for a historical covariate. See SI 04_Alternative_AGR_covariates.md.\n")

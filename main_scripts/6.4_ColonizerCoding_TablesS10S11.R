# =============================================================================
#  Coloniser identity: coding rule comparison and robustness
#
#  The manuscript assigns each country the power administering it AT INDEPENDENCE.
#  This script documents why, by comparing that rule with the longest-lasting rule
#  used in the earlier version, and reports the robustness battery for the
#  coloniser association. Produces SI Tables S10 and S11.
#
#  Run from the repository root:  source("main_scripts/4.3.6_ColonizerCoding.R")
#  Requires: MASS
# =============================================================================
suppressMessages(library(MASS))

ROOT <- if (dir.exists("Data_Colonial")) "." else ".."
f    <- paste0(ROOT, "/Data_Colonial/")
folder <- f
sm   <- read.csv(paste0(f, "SettlerMortality.csv"), sep=";", dec=",")

sm <- sm[!duplicated(sm$ccodealp), c("ccodealp","ajr_settmort")]; names(sm) <- c("iso3","settmort")
cr <- read.csv(paste0(f,"Colonies.csv")); names(cr)[4] <- "power"
pw <- c("Belgium","France","Germany","Italy","Netherlands","Portugal","Spain","United Kingdom")
sub <- cr[cr$power %in% pw, ]

# --- the two coding rules ----------------------------------------------------
cnt  <- aggregate(list(n=sub$Year), by=sub[,c("Code","power")], FUN=length)
lgst <- cnt[order(-cnt$n), ]; lgst <- lgst[!duplicated(lgst$Code), c("Code","power")]
names(lgst) <- c("iso3","col_longest")
lyr  <- aggregate(list(y=sub$Year), by=sub[,c("Code","power")], FUN=max)
last <- lyr[order(-lyr$y), ]; last <- last[!duplicated(last$Code), c("Code","power")]
names(last) <- c("iso3","col_indep")
indep <- aggregate(list(indep=sub$Year), by=list(iso3=sub$Code), FUN=max); indep$indep <- indep$indep + 1

pe <- read.csv(paste0(f,"Prienr.csv")); names(pe)[4] <- "prienr1900"
pe <- pe[pe$Year==1900 & !pe$Code %in% c("","OWID_WRL") & !is.na(pe$Code), c("Code","prienr1900")]
names(pe)[1] <- "iso3"

rds <- paste0(ROOT, "/main_scripts/4_RankedClusters.RDS")
cl  <- readRDS(rds); names(cl)[names(cl)=="SPI_countrycode"] <- "iso3"
reg <- read.csv(paste0(ROOT, "/main_scripts/extendedDataScaled.csv"))[,c("SPI_countrycode","Region")]

reg <- reg[!duplicated(reg$SPI_countrycode),]; names(reg) <- c("iso3","Region")

build <- function(colmap, nm){
  fu <- merge(merge(sm, colmap, by="iso3", all=TRUE), pe, by="iso3", all=TRUE)
  names(fu)[names(fu)==nm] <- "colonizer"
  tb <- table(fu$colonizer); fu$colonizer[fu$colonizer %in% names(tb)[tb<10]] <- "Other"
  x <- merge(merge(merge(fu, reg, by="iso3", all.x=TRUE), indep, by="iso3", all.x=TRUE),
             cl[,c("iso3","DLSFew_coverage")], by="iso3")
  names(x)[names(x)=="DLSFew_coverage"] <- "cluster"
  x <- x[complete.cases(x[,c("settmort","prienr1900","cluster","colonizer")]),]
  x$colonizer <- relevel(factor(x$colonizer),"United Kingdom")
  x$ord <- factor(x$cluster, ordered=TRUE); x }

dL <- build(lgst, "col_longest")     # longest-lasting rule
dI <- build(last, "col_indep")       # at-independence rule (manuscript)

# --- which countries change --------------------------------------------------
cmp <- merge(lgst, last, by="iso3"); cmp <- cmp[cmp$iso3 %in% dI$iso3, ]
dis <- cmp[as.character(cmp$col_longest) != as.character(cmp$col_indep), ]
dis <- merge(dis, dI[,c("iso3","cluster")], by="iso3")
cat("=== Countries where the two coding rules disagree ===\n")
print(dis[order(dis$cluster), ], row.names=FALSE)
cat(sprintf("\n%d of %d countries change (%.0f%%)\n\n", nrow(dis), nrow(cmp), 100*nrow(dis)/nrow(cmp)))

O  <- function(fm,x) polr(fm, data=x, Hess=TRUE)
pv <- function(a,b,df) pchisq(deviance(b)-deviance(a), df, lower.tail=FALSE)

report <- function(x, lab){
  k <- nlevels(x$colonizer)-1
  m2 <- O(ord~settmort+prienr1900+colonizer, x); m1 <- O(ord~settmort+prienr1900, x)
  cat(sprintf("### %s   N = %d\n", lab, nrow(x))); print(table(x$colonizer))
  cat(sprintf("  coloniser alone            : chi2=%5.2f df=%d p=%.4f\n",
      deviance(O(ord~1,x))-deviance(O(ord~colonizer,x)), k, pv(O(ord~colonizer,x),O(ord~1,x),k)))
  cat(sprintf("  + settmort + enrolment     : chi2=%5.2f df=%d p=%.4f   AIC %.1f vs %.1f\n",
      deviance(m1)-deviance(m2), k, pv(m2,m1,k), AIC(m2), AIC(m1)))
  r <- x[!is.na(x$Region),]; r$ord <- factor(r$cluster, ordered=TRUE)
  cat(sprintf("  + Region fixed effects     : p=%.4f\n",
      pv(O(ord~settmort+prienr1900+factor(Region)+colonizer,r),
         O(ord~settmort+prienr1900+factor(Region),r), k)))
  cat(sprintf("  + independence year        : p=%.4f\n",
      pv(O(ord~settmort+prienr1900+indep+colonizer,x), O(ord~settmort+prienr1900+indep,x), k)))
  for (nm in list(c("excl. settler colonies","AUS|CAN|NZL|USA"), c("excl. Spain","Spain"), c("excl. France","France"))){
    y <- if (nm[1]=="excl. settler colonies") x[!x$iso3 %in% c("AUS","CAN","NZL","USA"),] else x[x$colonizer!=nm[2],]
    y$ord <- factor(y$cluster, ordered=TRUE); y$colonizer <- droplevels(y$colonizer); kk <- nlevels(y$colonizer)-1
    if (kk < 1) next
    cat(sprintf("  %-26s : N=%2d p=%.4f\n", nm[1], nrow(y),
        pv(O(ord~settmort+prienr1900+colonizer,y), O(ord~settmort+prienr1900,y), kk))) }
  ssa <- r[r$Region=="Sub-Saharan Africa",]
  kw <- kruskal.test(cluster ~ droplevels(colonizer), data=ssa)
  cat(sprintf("  within sub-Saharan Africa  : N=%d Kruskal-Wallis chi2=%.2f df=%d p=%.3f\n",
      nrow(ssa), kw$statistic, kw$parameter, kw$p.value))
  print(rbind(n = table(droplevels(ssa$colonizer)),
              mean_cluster  = round(tapply(ssa$cluster,   droplevels(ssa$colonizer), mean), 2),
              mean_settmort = round(tapply(ssa$settmort,  droplevels(ssa$colonizer), mean), 2),
              mean_enrol    = round(tapply(ssa$prienr1900,droplevels(ssa$colonizer), mean), 2)))
  cat("\n") }

cat("=== SI TABLE S10 : the two coding rules ===\n\n")
report(dL, "LONGEST-LASTING power (earlier version)")
report(dI, "POWER AT INDEPENDENCE (manuscript)")
cat("Median independence year by coloniser, at-independence coding:\n")
print(tapply(dI$indep, dI$colonizer, median))

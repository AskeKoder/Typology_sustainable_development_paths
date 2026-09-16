# =============================================================================
#  Figure 5 | Colonial-era conditions and position along the development ranking
#
#  Replaces the multinomial predicted-probability heat map.
#    (a) the observed data, by cluster
#    (b) the fitted ordinal model for the two continuous covariates
#    (c) coloniser identity, as estimated differences from British colonies
#
#  The model drawn is the FULL specification (settler mortality + 1900 enrolment
#  + coloniser identity). This matters: the two-covariate model without coloniser
#  fails the proportional-odds (Brant) test, while the full model passes it
#  wherever it can be computed (see 4.3.7_ProportionalOddsCheck.R).
#
#  Run from the repository root:  source("main_scripts/4.3.3_ColonialFigure.R")
#  Requires: ggplot2, patchwork, MASS
# =============================================================================
library(ggplot2); library(patchwork); library(MASS)

# ---- 0. Palette (validated all-pairs: worst CVD dE 9.2, normal-vision 24.0) --
PAL <- list(uk="#2a78d6", fr="#eb6834", es="#1baf7a", other="#898781",
            line="#2a78d6", band="#9ec5f4",
            ink="#0b0b0b", ink2="#52514e", muted="#898781",
            grid="#e1e0d9", axis="#c3c2b7", surface="#fcfcfb")

# ---- 1. Estimation frame (at-independence coloniser coding) -----------------

ROOT <- if (dir.exists("Data_Colonial")) "." else ".."
f    <- paste0(ROOT, "/Data_Colonial/")
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
rds <- if (file.exists("main_scripts/4_RankedClusters.RDS")) "main_scripts/4_RankedClusters.RDS" else "4_RankedClusters.RDS"
cl <- readRDS(rds); names(cl)[names(cl)=="SPI_countrycode"] <- "iso3"
fu <- merge(merge(sm,col,by="iso3",all=TRUE), pe, by="iso3", all=TRUE)
tb <- table(fu$colonizer); fu$colonizer[fu$colonizer %in% names(tb)[tb<10]] <- "Other"
d <- merge(fu, cl[,c("iso3","DLSFew_coverage")], by="iso3"); names(d)[ncol(d)] <- "cluster"
d <- d[complete.cases(d[,c("settmort","prienr1900","cluster","colonizer")]),]
d$colonizer <- relevel(factor(d$colonizer), "United Kingdom")
d$ord <- factor(d$cluster, ordered=TRUE)
d$settler <- d$iso3 %in% c("AUS","CAN","NZL","USA")
cat("N =", nrow(d), "| clusters:", nlevels(d$ord), "\n")

# ---- 2. Fit the full ordinal model ------------------------------------------
fit <- polr(ord ~ settmort + prienr1900 + colonizer, data=d, Hess=TRUE)
print(summary(fit))

# ---- 3. Expected cluster + 95% interval, by parametric simulation ------------
lev <- as.numeric(levels(d$ord)); nb <- length(coef(fit)); nz <- length(fit$zeta)
set.seed(42); sims <- mvrnorm(4000, c(coef(fit), fit$zeta), vcov(fit))
Erank <- function(eta, zeta){ cum <- rbind(plogis(outer(zeta, eta, "-")), 1)
  colSums(rbind(cum[1,,drop=FALSE], diff(cum)) * lev) }
curve1 <- function(v, other, n=200){
  g <- seq(min(d[[v]]), max(d[[v]]), length.out=n); oth <- median(d[[other]])
  X <- cbind(if (v=="settmort") g else oth, if (v=="settmort") oth else g, 0, 0, 0)  # coloniser = UK
  fitv <- Erank(as.numeric(X %*% coef(fit)), fit$zeta)
  dr <- vapply(seq_len(nrow(sims)), function(i)
        Erank(as.numeric(X %*% sims[i,1:nb]), sims[i,(nb+1):(nb+nz)]), numeric(n))
  data.frame(x=g, fit=fitv, lo=apply(dr,1,quantile,.025), hi=apply(dr,1,quantile,.975), var=v) }
curves <- rbind(curve1("settmort","prienr1900"), curve1("prienr1900","settmort"))

# ---- 4. Long format & labels -------------------------------------------------
long <- rbind(data.frame(cluster=d$cluster, colonizer=d$colonizer, settler=d$settler, x=d$settmort,   var="settmort"),
              data.frame(cluster=d$cluster, colonizer=d$colonizer, settler=d$settler, x=d$prienr1900, var="prienr1900"))
meds <- aggregate(x ~ cluster + var, long, median)
vlab <- c(settmort="Log European settler mortality", prienr1900="Primary school enrolment, 1900 (%)")
long$var<-factor(vlab[long$var],levels=vlab); meds$var<-factor(vlab[meds$var],levels=vlab)
curves$var<-factor(vlab[curves$var],levels=vlab)
cols <- c("United Kingdom"=PAL$uk, "France"=PAL$fr, "Spain"=PAL$es, "Other"=PAL$other)

th <- theme_minimal(base_size=10) + theme(
  text=element_text(colour=PAL$ink2),
  plot.background=element_rect(fill=PAL$surface,colour=NA),
  panel.background=element_rect(fill=PAL$surface,colour=NA),
  panel.grid.major=element_line(colour=PAL$grid, linewidth=.3), panel.grid.minor=element_blank(),
  axis.line.x=element_line(colour=PAL$axis, linewidth=.4), axis.ticks=element_blank(),
  axis.text=element_text(colour=PAL$muted,size=8.5), axis.title=element_text(colour=PAL$ink2,size=9),
  strip.text=element_text(colour=PAL$ink,size=9.5,hjust=0,margin=margin(b=4)),
  plot.title=element_text(colour=PAL$ink,size=10.5,face="bold",margin=margin(b=2)),
  plot.subtitle=element_text(colour=PAL$muted,size=8.2,margin=margin(b=8)),
  legend.position="bottom", legend.title=element_blank(),
  legend.text=element_text(colour=PAL$ink2,size=8.5), legend.key.size=unit(9,"pt"),
  plot.margin=margin(6,10,4,6))
yb <- sort(unique(d$cluster))

# ---- 5. Panel a : observed data ---------------------------------------------
set.seed(7)
pA <- ggplot(long, aes(x=x, y=cluster)) +
  geom_point(aes(fill=colonizer, colour=settler), shape=21, size=2.4, stroke=.7,
             position=position_jitter(width=0, height=.16)) +
  geom_segment(data=meds, aes(x=x,xend=x,y=cluster-.38,yend=cluster+.38),
               colour=PAL$surface, linewidth=1.9, inherit.aes=FALSE) +
  geom_segment(data=meds, aes(x=x,xend=x,y=cluster-.38,yend=cluster+.38),
               colour=PAL$ink, linewidth=.75, inherit.aes=FALSE) +
  scale_fill_manual(values=cols) +
  scale_colour_manual(values=c("FALSE"=PAL$surface,"TRUE"=PAL$ink), guide="none") +
  scale_y_reverse(breaks=yb, limits=c(11.7,.4)) +
  facet_wrap(~var, scales="free_x", strip.position="bottom") +
  labs(title="a  Observed colonial-era conditions, by development cluster",
       subtitle="One point per country; vertical rule is the cluster median. Dark outline marks the four settler colonies\n(Australia, Canada, New Zealand, United States). Cluster 2 has no country in the historical sample.",
       x=NULL, y="Development cluster\n(1 = highest HDI)") +
  th + theme(strip.placement="outside")

# ---- 6. Panel b : fitted model, continuous covariates ------------------------
pB <- ggplot(curves, aes(x=x)) +
  geom_ribbon(aes(ymin=lo,ymax=hi), fill=PAL$band, alpha=.45) +
  geom_line(aes(y=fit), colour=PAL$line, linewidth=.9) +
  geom_rug(data=long, aes(x=x), inherit.aes=FALSE, sides="b",
           colour=PAL$muted, alpha=.55, length=unit(3,"pt")) +
  scale_y_reverse(breaks=yb, limits=c(11.7,.4)) +
  facet_wrap(~var, scales="free_x", strip.position="bottom") +
  labs(title="b  Fitted ordinal model: expected cluster",
       subtitle="Expected cluster with 95% interval, coloniser held at the United Kingdom and the other covariate at its median.\nTicks mark observed values; curves are not drawn beyond them. Same vertical scale as panel a.",
       x=NULL, y="Expected cluster") +
  th + theme(strip.placement="outside", legend.position="none")

# ---- 7. Panel c : coloniser effects ------------------------------------------
ct <- coef(summary(fit)); k <- grep("^colonizer", rownames(ct))
eff <- data.frame(colonizer=sub("^colonizer","",rownames(ct)[k]),
                  est=ct[k,1], lo=ct[k,1]-1.96*ct[k,2], hi=ct[k,1]+1.96*ct[k,2])
eff$colonizer <- factor(eff$colonizer, levels=rev(c("France","Other","Spain")))
pC <- ggplot(eff, aes(x=est, y=colonizer, colour=colonizer)) +
  geom_vline(xintercept=0, colour=PAL$axis, linewidth=.4) +
  geom_errorbar(aes(xmin=lo, xmax=hi), orientation="y", width=0, linewidth=.9) +
  geom_point(size=3) +
  scale_colour_manual(values=cols, guide="none") +
  labs(title="c  Coloniser identity, relative to British colonies",
       subtitle="Ordinal coefficients with 95% intervals. Negative = higher-ranked cluster than British colonies at the same\nsettler mortality and 1900 enrolment. Joint likelihood-ratio test: chi-sq = 12.45, df = 3, p = 0.006.",
       x="Change in log-odds of a lower-ranked cluster", y=NULL) +
  th + theme(panel.grid.major.y=element_blank())

#`+ plot_annotation(theme=)` sets the background of the combined canvas; each
#panel already carries it through `th`.  (The patchwork `&` operator is avoided:
#with ggplot2 >= 4.0 it fails unless patchwork >= 1.3.1 is installed.)
fig <- pA / pB / pC + plot_layout(heights=c(1.15,1,.55)) +
  plot_annotation(theme=theme(plot.background=element_rect(fill=PAL$surface, colour=NA)))
s <- paste0(ROOT, "/Figures/Article")
dir.create(s, showWarnings = FALSE, recursive = TRUE)
ggsave(paste0(s, "/Fig6_historical_covariates.png"), fig, width=8.4, height=9.6, dpi=400, bg=PAL$surface)
ggsave(paste0(s, "/Fig6_historical_covariates.pdf"), fig, width=8.4, height=9.6, bg=PAL$surface)
cat("\nWritten:", normalizePath(s), "/Fig5_historical_covariates.{png,pdf}\n")

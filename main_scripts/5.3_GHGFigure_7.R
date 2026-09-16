# =============================================================================
#  4.8_GHGFigure.R  -  Figure 7 | Emissions level and trend, by provisioning rank
#
#  TWO PANELS, ONE SHARED AXIS.  The previous version of this exhibit put mean
#  per-capita impact on a left axis and the annual rate of change on a right
#  axis of the same panel.  Two measures on two scales in one frame is the
#  commonest chart error: the reader cannot tell which axis a mark belongs to,
#  and the apparent relation between bars and points is an artefact of how the
#  two scales were aligned.  They are separated here.
#
#    (a) LEVEL - where each country stands.  Ordered, tight, and essentially
#        development level restated: the ranking adds nothing to it (p = 0.51).
#    (b) TREND - which way each country is going.  Ordered but overlapping:
#        the ranking predicts it beyond development level (p = 0.002), yet 74%
#        of the variance is within clusters, and ten of eleven clusters hold
#        countries moving in opposite directions.
#
#  Run from main_scripts/ or the repository root.  Requires ggplot2, patchwork, readxl.
# =============================================================================
library(ggplot2); library(patchwork); suppressMessages(library(readxl))
here   <- if (file.exists("2_ImputedData.csv")) "" else "main_scripts/"
root   <- if (here == "") ".." else "."
OUTD   <- file.path(root, "Figures/Article"); dir.create(OUTD, showWarnings = FALSE)
xlsx   <- file.path(root, "Data_Environmental", "GCSI_59a_Per_capita_2000-2020_05162025.xlsx")
if (!file.exists(xlsx)) stop("Cannot find ", xlsx, " - the GCSI footprint workbook used by 4.7.")

PAL <- list(fall="#2a78d6", rise="#eb6834", ink="#0b0b0b", ink2="#52514e",
            muted="#898781", grid="#e1e0d9", axis="#c3c2b7", surface="#fcfcfb")
th <- theme_minimal(base_size = 10) + theme(
  text=element_text(colour=PAL$ink2),
  plot.background=element_rect(fill=PAL$surface,colour=NA),
  panel.background=element_rect(fill=PAL$surface,colour=NA),
  panel.grid.major.x=element_line(colour=PAL$grid, linewidth=.3),
  panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
  axis.line.x=element_line(colour=PAL$axis, linewidth=.4), axis.ticks=element_blank(),
  axis.text=element_text(colour=PAL$muted,size=8.5),
  axis.title=element_text(colour=PAL$ink2,size=9),
  plot.title=element_text(colour=PAL$ink,size=10,face="bold",margin=margin(b=1)),
  plot.subtitle=element_text(colour=PAL$muted,size=8.2,margin=margin(b=7)),
  legend.position="bottom", legend.title=element_blank(),
  legend.text=element_text(colour=PAL$ink2,size=8.5), legend.key.size=unit(9,"pt"),
  plot.margin=margin(6,10,4,4))

# ---- data -------------------------------------------------------------------
#Same extraction as 4.7_GHGRegimes.R: consumption-based GHG (AR6 GWP-100),
#per capita = Sum_Value / Population, 2000-2020.
raw <- as.data.frame(read_excel(xlsx, sheet = "Sheet1"))
g   <- raw[raw$Indicator == "ghg_ar6", c("Region_acronyms","Year","Sum_Value","Population")]
names(g) <- c("iso3","year","total","pop")
g$year  <- as.integer(g$year)
g$total <- suppressWarnings(as.numeric(g$total)); g$pop <- suppressWarnings(as.numeric(g$pop))
g <- g[!is.na(g$total) & !is.na(g$pop) & g$pop > 0 & g$year>=2000 & g$year<=2020,]
g$percap <- g$total / g$pop
g <- g[g$percap > 0,]
iso <- sort(unique(g$iso3)); rate <- setNames(rep(NA_real_,length(iso)), iso)
for (i in iso) { s <- g[g$iso3==i,]; if (nrow(s) < 10) next
  rate[i] <- (exp(coef(lm(log(percap) ~ year, data=s))[2]) - 1) * 100 }
lev <- tapply(g$percap, g$iso3, mean) * 1e3          # t CO2e per person

rk_file <- if (file.exists(paste0(here,"4_RankedClusters.RDS")))
             paste0(here,"4_RankedClusters.RDS") else "4_RankedClusters.RDS"
ranked <- suppressWarnings(readRDS(rk_file))
k  <- intersect(names(rate)[!is.na(rate)], ranked$SPI_countrycode)
mm <- match(k, ranked$SPI_countrycode)
d <- data.frame(iso3=k, lev=as.numeric(lev[k]), rate=as.numeric(rate[k]),
                cl=factor(ranked$DLSFew_coverage[mm], levels=1:11))
d$cly <- factor(d$cl, levels=rev(levels(d$cl)))       # cluster 1 at the TOP
d$dir <- factor(ifelse(d$rate < 0, "Falling", "Rising"), levels=c("Falling","Rising"))
set.seed(11); d$jit <- as.numeric(d$cly) + runif(nrow(d), -.28, .28)

sm <- data.frame(cly = levels(d$cly),
                 lev = tapply(d$lev,  d$cly, median)[levels(d$cly)],
                 rate= tapply(d$rate, d$cly, mean)[levels(d$cly)])
sm$y <- match(sm$cly, levels(d$cly))

dot <- function(...) geom_point(..., shape=21, size=2.1, stroke=.45, colour=PAL$surface)
tick <- function(df, xv) list(
  geom_segment(data=df, aes(x=.data[[xv]], xend=.data[[xv]], y=y-.40, yend=y+.40),
               colour=PAL$ink, linewidth=.8, inherit.aes=FALSE))

# ---- (a) level --------------------------------------------------------------
pa <- ggplot(d, aes(lev, jit)) +
  dot(fill = PAL$muted) +
  tick(sm, "lev") +
  scale_x_log10(breaks=c(1,3,10,30), labels=c("1","3","10","30")) +
  scale_y_continuous(breaks=1:11, labels=rev(paste("Cluster", 1:11)),
                     limits=c(.4,11.6), expand=c(0,0)) +
  labs(title="a   How much each country emits",
       subtitle="Mean per-capita footprint, 2000-2020. Log scale.\nOrdered and tight: this is development level restated.",
       x="Tonnes CO2e per person, per year", y=NULL) + th

# ---- (b) trend --------------------------------------------------------------
pb <- ggplot(d, aes(rate, jit, fill=dir)) +
  annotate("segment", x=0, xend=0, y=.4, yend=11.6, colour=PAL$axis, linewidth=.5) +
  dot() + tick(sm, "rate") +
  scale_fill_manual(values=c(Falling=PAL$fall, Rising=PAL$rise)) +
  scale_x_continuous(breaks=seq(-4,5,1), labels=function(x) sprintf("%+d", x)) +
  scale_y_continuous(breaks=1:11, labels=NULL, limits=c(.4,11.6), expand=c(0,0)) +
  labs(title="b   Which way each country is going",
       subtitle="Annual change in that footprint, from a country-specific log-linear fit.\nOrdered but overlapping: ten of eleven clusters contain both.",
       x="Change per year (%)", y=NULL) + th

#The collected legend sits at the bottom via the patchwork-level theme below.
#(The patchwork `&` operator is avoided: with ggplot2 >= 4.0 it fails unless
#patchwork >= 1.3.1 is installed.)
fig <- (pa | pb) + plot_layout(widths=c(1,1.15), guides="collect")
fig <- fig + plot_annotation(
  title="Emissions level and emissions trend across the provisioning ranking",
  subtitle="Each dot is one of 153 countries; the black rule is the cluster median (a) or mean (b).",
  caption=paste("Consumption-based greenhouse gases (AR6 GWP-100), 2000-2020. Panel a: the ranking adds nothing to a cubic in development level (p = 0.51).",
                "Panel b: it does (p = 0.002), and survives region\nfixed effects and control for the starting level, but 74% of the variance in rates lies within clusters rather than between them.",
                "Twenty countries lack footprint accounts."),
  theme=theme(plot.background=element_rect(fill=PAL$surface,colour=NA),
              legend.position="bottom",
              plot.title=element_text(colour=PAL$ink,size=12,face="bold"),
              plot.subtitle=element_text(colour=PAL$muted,size=9,margin=margin(b=6)),
              plot.caption=element_text(colour=PAL$muted,size=7.6,hjust=0,margin=margin(t=8))))

W <- 9.2; H <- 5.6
ggsave(file.path(OUTD,"Fig7_GHG_level_and_trend.png"), fig, width=W, height=H, dpi=400, bg=PAL$surface)
ggsave(file.path(OUTD,"Fig7_GHG_level_and_trend.pdf"), fig, width=W, height=H, device=cairo_pdf, bg=PAL$surface)
cat(sprintf("Written: %s (%.0f x %.0f px, 400 dpi)\n",
            file.path(OUTD,"Fig7_GHG_level_and_trend.png"), W*400, H*400))
cat(sprintf("N = %d countries; %d falling, %d rising\n", nrow(d), sum(d$rate<0), sum(d$rate>=0)))

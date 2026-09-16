# =============================================================================
#  Figure 5 | Improvement regimes: what each group of countries actually moved
#
#  Six regimes of countries, 2000-2020, defined by which components of material
#  provisioning they advanced *faster or slower than their starting level would
#  predict*. The figure then shows, for each regime, the mean improvement rate
#  of the eleven indicators as a DEVIATION from the all-country mean rate for
#  that indicator, so the reader sees what is distinctive about the regime
#  rather than the global pattern ("telephones are fast everywhere").
#
#  Run from the repository's main_scripts/ directory:
#      source("<path>/4.6_RegimeFigure.R")
#
#  Requires: ggplot2, patchwork, NbClust  (all already used in this repository)
#  Writes:   <OUT>/Fig2_improvement_regimes.png   (400 dpi)
#            <OUT>/Fig2_improvement_regimes.pdf   (vector)
#            <OUT>/Figure2_caption.txt            (draft caption)
#            <OUT>/6_ImprovementRegimes.RDS       (only if the shared result is absent)
#
#  Visual language, palette (PAL) and theme (th) are copied verbatim from
#  main_scripts/4.3.3_ColonialFigure.R so that Figure 2 and Figure 5 read as a
#  set. Nothing in the repository is modified by this script.
# =============================================================================
library(ggplot2); library(patchwork)

# Output directory, resolved the same way 4.5_ImprovementRegimes.R resolves it, so
# the figure lands inside the repository whether the script is sourced from
# main_scripts/ or from the repository root.  Override by setting OUT beforehand.
if (!exists("OUT")) {
  .here <- if (file.exists("2_ImputedData.csv")) "" else "main_scripts/"
  OUT   <- if (.here == "") "../Figures/Article" else "Revision_output"
}
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- 0. Palette and theme (verbatim from 4.3.3_ColonialFigure.R) ------------
PAL <- list(uk="#2a78d6", fr="#eb6834", es="#1baf7a", other="#898781",
            line="#2a78d6", band="#9ec5f4",
            ink="#0b0b0b", ink2="#52514e", muted="#898781",
            grid="#e1e0d9", axis="#c3c2b7", surface="#fcfcfb")

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

# ---- 1. The eleven material-provisioning indicators -------------------------
#  All are stored on a 0-100 "higher = better" scale in 2_ImputedData.csv.
#  Nutritional_deficiencies and No_access_to_a_handwashing_facility have
#  deficit-sounding names but are already stored higher = better. Share_Slums
#  is the one indicator stored higher = worse, so it is negated.
IND <- c("Nutritional_deficiencies","No_access_to_a_handwashing_facility",
         "Improved_sanitation","Improved_water_source",
         "Prevalence_of_cooking_with_coalperbiomass","Access_to_electricity",
         "Mean_years_of_schooling","Internet_users",
         "Mobile_and_landline_telephone_subscriptions","Universal_health_coverage",
         "Share_Slums")
NEG <- "Share_Slums"

# =============================================================================
#  PART 1 | The regimes
# =============================================================================
rds_repo  <- "Revision_output/6_ImprovementRegimes.RDS"
rds_repo2 <- "../Revision_output/6_ImprovementRegimes.RDS"
rds_cache <- file.path(OUT, "6_ImprovementRegimes.RDS")

# ---- 1a. Mirkin distance between two partitions of the same objects ---------
#  M(P,Q) = #pairs joined in P but split in Q  +  #pairs split in P but joined
#  in Q, counted over ORDERED pairs.  This is deliberately twice the unordered
#  N10 + N01 that 4.5_ImprovementRegimes.R and 3_ClusterExperimentation_ver2.R
#  use.  The factor is a monotone rescaling and cannot change which partition is
#  the medoid, so do not "fix" one of the two to match the other:
#      M = sum_i n_i.^2 + sum_j n_.j^2 - 2 * sum_ij n_ij^2
mirkin <- function(a, b) {
  tb <- table(a, b)
  sum(rowSums(tb)^2) + sum(colSums(tb)^2) - 2 * sum(tb^2)
}

# ---- 1b. Full computation ---------------------------------------------------
compute_regimes <- function() {
  library(NbClust)
  message("Computing improvement regimes from 2_ImputedData.csv ...")
  #Resolved like every other path in this script: main_scripts/ or the root.
  .imp <- if (file.exists("2_ImputedData.csv")) "2_ImputedData.csv" else "main_scripts/2_ImputedData.csv"
  d <- read.csv(.imp, encoding = "UTF-8")
  stopifnot(all(IND %in% names(d)))
  d[[NEG]] <- -d[[NEG]]                       # only Share_Slums is negated

  imps <- sort(unique(d$.imp))
  ccs  <- sort(unique(d$SPI_countrycode))
  nC <- length(ccs); nI <- length(IND); nM <- length(imps)
  cmap <- unique(d[, c("SPI_countrycode","Country")])
  cmap <- cmap[order(cmap$SPI_countrycode), ]

  parts  <- matrix(NA_integer_, nC, nM, dimnames = list(ccs, paste0("imp", imps)))
  kvec   <- integer(nM)
  SLlist <- STlist <- LVlist <- vector("list", nM)

  for (m in seq_len(nM)) {
    dm <- d[d$.imp == imps[m], ]
    dm <- dm[order(dm$SPI_countrycode, dm$SPI_year), ]
    SL <- matrix(NA_real_, nC, nI, dimnames = list(ccs, IND))   # slope, pts/yr
    ST <- SL                                                    # value in 2000
    LV <- SL                                                    # mean level 2000-2020
    for (j in seq_len(nI)) {
      v <- dm[[IND[j]]]
      for (i in seq_len(nC)) {
        sel <- dm$SPI_countrycode == ccs[i]
        yy <- dm$SPI_year[sel]; vv <- v[sel]
        SL[i, j] <- unname(coef(lm(vv ~ yy))[2])
        ST[i, j] <- vv[which(yy == 2000)]
        LV[i, j] <- mean(vv)
      }
    }
    # Improvement net of the level it started from (cubic in the 2000 value)
    RS <- matrix(NA_real_, nC, nI, dimnames = dimnames(SL))
    for (j in seq_len(nI)) RS[, j] <- resid(lm(SL[, j] ~ poly(ST[, j], 3)))
    Z <- scale(RS)
    D <- dist(Z)
    nb <- NbClust(data = NULL, diss = D, distance = NULL,
                  min.nc = 2, max.nc = 40, index = "cindex", method = "ward.D2")
    k <- as.integer(nb$Best.nc[["Number_clusters"]])
    kvec[m]  <- k
    parts[, m] <- cutree(hclust(D, "ward.D2"), k)
    SLlist[[m]] <- SL; STlist[[m]] <- ST; LVlist[[m]] <- LV
    message(sprintf("  imputation %2d/%d : k = %d", m, nM, k))
  }

  # ---- Ensemble medoid: the partition closest to all the others -------------
  DM <- matrix(0, nM, nM)
  for (a in seq_len(nM)) for (b in seq_len(nM)) DM[a, b] <- mirkin(parts[, a], parts[, b])
  med <- which.min(rowSums(DM))
  cl  <- parts[, med]
  message(sprintf("  medoid partition = imputation %d, k = %d (k across imputations: %s)",
                  med, kvec[med], paste(kvec, collapse = ",")))

  # ---- Pooled (across-imputation) rates and levels --------------------------
  SLbar <- Reduce(`+`, SLlist) / nM
  STbar <- Reduce(`+`, STlist) / nM
  LVbar <- Reduce(`+`, LVlist) / nM
  # Development level = mean of the standardised year-2000 values (higher = better)
  level <- rowMeans(scale(STbar)); names(level) <- ccs

  # ---- Order regimes by mean development level, descending ------------------
  ordg   <- order(tapply(level, cl, mean), decreasing = TRUE)
  regime <- match(cl, ordg); names(regime) <- names(cl)

  list(regime = regime, iso3 = ccs, country = cmap, level = level,
       SL = SLbar, ST = STbar, LV = LVbar, partitions = parts, k = kvec,
       medoid = as.integer(med), mirkin = DM, indicators = IND)
}

# ---- 1c. Rates: computed once, then cached (they are needed either way) -----
#  The cache is stamped with the modification time of the input file and with the
#  indicator list, and is discarded if either has moved, so the figure can never
#  be drawn from rates that no longer correspond to the data.
rates_cache <- file.path(OUT, ".regime_rates_cache.RDS")
.stamp <- list(mtime = file.mtime(if (file.exists("2_ImputedData.csv"))
                                    "2_ImputedData.csv" else "main_scripts/2_ImputedData.csv"),
               ind = IND)
FULL <- NULL
if (file.exists(rates_cache)) {
  .cached <- readRDS(rates_cache)
  if (identical(.cached$stamp, .stamp)) FULL <- .cached
  else message("Cached improvement rates are stale (input file or indicator list changed); recomputing.")
}
if (is.null(FULL)) {
  FULL <- compute_regimes()
  FULL$stamp <- .stamp
  saveRDS(FULL, rates_cache)
}
SLbar <- FULL$SL; cmap <- FULL$country

# ---- 1d. Regimes: prefer the analysis script's shared result ----------------
read_shared <- function(obj) {
  # Accept the plausible shapes of Revision_output/6_ImprovementRegimes.RDS
  if (is.data.frame(obj)) {
    gcol <- grep("regime", names(obj), ignore.case = TRUE)[1]
    icol <- grep("iso3|countrycode", names(obj), ignore.case = TRUE)[1]
    lcol <- grep("level", names(obj), ignore.case = TRUE)[1]
    if (!is.na(gcol) && !is.na(icol)) {
      v <- as.integer(obj[[gcol]]); names(v) <- as.character(obj[[icol]])
      lv <- if (!is.na(lcol)) setNames(as.numeric(obj[[lcol]]), names(v)) else NULL
      return(list(regime = v, level = lv))
    }
  }
  if (is.list(obj) && !is.null(obj$regime) && !is.null(names(obj$regime)))
    return(list(regime = obj$regime, level = obj$level))
  if (is.atomic(obj) && !is.null(names(obj)))
    return(list(regime = setNames(as.integer(obj), names(obj)), level = NULL))
  NULL
}

SRC <- NULL; SRC_FROM <- ""
for (p in c(rds_repo, rds_repo2)) {
  if (is.null(SRC) && file.exists(p)) {
    cand <- read_shared(readRDS(p))
    if (!is.null(cand)) { SRC <- cand; SRC_FROM <- normalizePath(p) }
    else warning("Found ", p, " but could not read a regime vector from it.")
  }
}
if (!is.null(SRC)) {
  cat("USING EXISTING RESULT:", SRC_FROM,
      "\n  -> the regimes were NOT recomputed; only the improvement rates were rebuilt",
      "\n     from 2_ImputedData.csv so that the figure can be drawn.\n")
  agree <- table(FULL$regime[names(SRC$regime)], SRC$regime)
  cat("  agreement with an independent recomputation:",
      sprintf("%.1f%% of countries in matching groups\n",
              100 * sum(apply(agree, 2, max)) / sum(agree)))
} else {
  cat("No Revision_output/6_ImprovementRegimes.RDS found - regimes computed here.\n")
  SRC <- list(regime = FULL$regime, level = FULL$level)
  saveRDS(FULL, rds_cache); cat("Saved:", rds_cache, "\n")
}

regime <- SRC$regime
level  <- if (!is.null(SRC$level)) SRC$level else FULL$level[names(regime)]
SLbar  <- SLbar[names(regime), , drop = FALSE]
G <- max(regime, na.rm = TRUE)

# =============================================================================
#  PART 2 | Regime profiles: deviation from the all-country mean rate
# =============================================================================
#  IMPORTANT.  Two quantities are available and they are NOT interchangeable:
#
#    raw      = regime mean improvement rate minus the all-country mean rate.
#               Interpretable in index points per year.  This is what Table 4 reports.
#    adjusted = regime mean of the rate AFTER residualising on a cubic in the
#               country's own year-2000 level.  This is the space the regimes were
#               actually formed in, and it is what the figure must show, because on
#               the raw scale a regime can look as though it advanced a component
#               fast when in fact it advanced it exactly as fast as its starting
#               position predicted.  Eleven of the sixty-six cells change sign
#               between the two, so the choice is not cosmetic.
#
#  The figure plots ADJUSTED.  Table 4 reports RAW.  The caption says so.
gmean <- colMeans(SLbar[, IND, drop = FALSE])                 # all-country mean rate
prof  <- t(sapply(seq_len(G), function(g)
            colMeans(SLbar[regime == g, IND, drop = FALSE]))) # regime mean rate
dev_raw <- sweep(prof, 2, gmean)                              # raw, for reference only

STbar <- FULL$ST[names(regime), , drop = FALSE]
RSbar <- SLbar[, IND, drop = FALSE]
for (j in IND) RSbar[, j] <- resid(lm(SLbar[, j] ~ poly(STbar[, j], 3)))
dev   <- t(sapply(seq_len(G), function(g)
            colMeans(RSbar[regime == g, , drop = FALSE])))    # the plotted quantity
colnames(dev) <- IND

cat("\n--- Sign disagreement between the raw and the adjusted profile ---\n")
cat(sprintf("%d of %d cells change sign; correlation across cells %.3f\n",
            sum(sign(dev_raw) != sign(dev)), length(dev),
            cor(as.vector(dev_raw), as.vector(dev))))
rownames(prof) <- rownames(dev) <- paste0("R", seq_len(G))

cat("\n--- All-country mean improvement rate (points per year, 2000-2020) ---\n")
print(round(sort(gmean, decreasing = TRUE), 3))
cat("\n--- Regime mean rate minus all-country mean rate (points per year) ---\n")
print(round(dev, 2))

# ---- Regime size, development level, and four example countries -------------
#  Examples are chosen deterministically: the countries at the 12th, 38th, 62nd
#  and 88th percentile of the regime's own development-level distribution, so
#  they span the regime rather than describing only its richest members.
SHORTEN <- c("Congo, Democratic Republic of" = "DR Congo",
             "Congo, Republic of"            = "Congo",
             "Central African Republic"      = "C. African Rep.",
             "Gambia, The"                   = "Gambia",
             "Republic of North Macedonia"   = "North Macedonia",
             "Korea, Democratic Republic of" = "North Korea",
             "Korea, Republic of"            = "South Korea",
             "Bosnia and Herzegovina"        = "Bosnia & Herz.",
             "United Arab Emirates"          = "UAE",
             "Sao Tome and Principe"         = "Sao Tome",
             "West Bank and Gaza"            = "Palestine",
             "Trinidad and Tobago"           = "Trinidad",
             "Papua New Guinea"              = "Papua N. Guinea",
             "Dominican Republic"            = "Dominican Rep.")
examples <- character(G); nvec <- integer(G); lvl <- numeric(G)
for (g in seq_len(G)) {
  cs <- names(regime)[regime == g]
  cs <- cs[order(-level[cs])]
  q  <- unique(round(quantile(seq_along(cs), c(.12,.38,.62,.88), type = 1)))
  nm <- cmap$Country[match(cs[q], cmap$SPI_countrycode)]
  nm <- unname(ifelse(nm %in% names(SHORTEN), SHORTEN[nm], nm))
  examples[g] <- paste(nm, collapse = "|")   # "|" marks the only legal break points
  nvec[g] <- length(cs); lvl[g] <- mean(level[cs])
}

# ---- Regime names, derived from the computed profiles -----------------------
#  Named by a deterministic rule on the deviation matrix, so that the labels
#  follow the data rather than an assumed ordering of the groups.
name_regimes <- function(dev) {
  #Names are derived from the ADJUSTED profile - the space the regimes were formed
  #in - not from the raw rates.  On the raw scale a regime can appear to advance a
  #component quickly when it advanced it exactly as fast as its starting position
  #predicted; eleven of the sixty-six cells change sign between the two scales, so
  #naming on the raw profile would misdescribe four of the six regimes.
  if (nrow(dev) != 6)
    stop("name_regimes() is written for exactly 6 regimes but received ", nrow(dev),
         ". Add or remove a naming rule, and adjust the 3 + 3 panel layout below.")
  left <- seq_len(nrow(dev)); nm <- character(nrow(dev)); why <- character(nrow(dev))
  take <- function(g, label, reason) { nm[g] <<- label; why[g] <<- reason; left <<- setdiff(left, g) }
  pick <- function(v) left[which.max(v[left])]

  # 1  above expectation on every front
  g <- pick(rowSums(dev > 0))
  take(g, "Advancing on every front",
       sprintf("%d of 11 indicators above expectation", sum(dev[g, ] > 0)))

  # 2  below expectation on almost every front
  g <- left[which.min(rowMeans(dev)[left])]
  take(g, "Falling behind on every front",
       sprintf("%d of 11 below expectation; mean deviation %+.2f",
               sum(dev[g, ] < 0), mean(dev[g, ])))

  # 3  the electricity outlier
  g <- pick(dev[, "Access_to_electricity"])
  take(g, "Electrification surge",
       sprintf("largest electricity deviation (%+.2f pts/yr adjusted)", dev[g, "Access_to_electricity"]))

  # 4  the regime that is not distinctive at all: it moved as its starting
  #    position predicted, on every indicator.  This is the reference group and
  #    must be named as such rather than on whichever indicator is marginally
  #    positive, which on the raw scale produced the misleading "connectivity only".
  g <- left[which.min(rowMeans(abs(dev))[left])]
  take(g, "Moving as expected",
       sprintf("largest deviation of any indicator only %+.2f; mean |deviation| %.2f",
               dev[g, which.max(abs(dev[g, ]))], mean(abs(dev[g, ]))))

  # 5  of the two that remain, the one whose shortfall is in food, fuel and health
  fh <- dev[, "Nutritional_deficiencies"] + dev[, "No_access_to_a_handwashing_facility"] +
        dev[, "Prevalence_of_cooking_with_coalperbiomass"] + dev[, "Universal_health_coverage"]
  g <- left[which.min(fh[left])]
  take(g, "Food, fuel and health lag",
       sprintf("nutrition %+.2f, handwashing %+.2f, clean fuels %+.2f, health %+.2f",
               dev[g, "Nutritional_deficiencies"], dev[g, "No_access_to_a_handwashing_facility"],
               dev[g, "Prevalence_of_cooking_with_coalperbiomass"], dev[g, "Universal_health_coverage"]))

  # 6  leftover: the shortfall is in water and energy
  take(left[1], "Water and energy lag",
       sprintf("water %+.2f, electricity %+.2f, housing %+.2f",
               dev[left[1], "Improved_water_source"], dev[left[1], "Access_to_electricity"],
               dev[left[1], "Share_Slums"]))
  list(name = nm, why = why)
}
NM    <- name_regimes(dev)
NAMES <- NM$name

cat("\n--- Regimes (ordered by mean development level, descending) ---\n")
for (g in seq_len(G))
  cat(sprintf("R%d  %-30s n = %3d  level = %+5.2f  |  %s\n     (%s)\n",
              g, NAMES[g], nvec[g], lvl[g], gsub("|", ", ", examples[g], fixed=TRUE), NM$why[g]))

# =============================================================================
#  PART 3 | The figure
# =============================================================================
# ---- Indicator order: grouped by decent-living-standards theme --------------
#  Grouping by theme rather than by speed is the point of the exhibit: the
#  question is which *domains* of provisioning a regime advanced, and a reader
#  can then see at a glance that, say, all of water-sanitation-hygiene moved
#  together. Ordering by speed would instead re-impose the global gradient that
#  the deviation is designed to remove.
LAB <- c(Nutritional_deficiencies                    = "Nutrition",
         Improved_water_source                       = "Drinking water",
         Improved_sanitation                         = "Sanitation",
         No_access_to_a_handwashing_facility         = "Handwashing",
         Access_to_electricity                       = "Electricity",
         Prevalence_of_cooking_with_coalperbiomass   = "Clean cooking",
         Share_Slums                                 = "Housing, non-slum",
         Universal_health_coverage                   = "Health coverage",
         Mean_years_of_schooling                     = "Schooling",
         Mobile_and_landline_telephone_subscriptions = "Telephone",
         Internet_users                              = "Internet")
THEME <- c(Nutritional_deficiencies                    = "Food",
           Improved_water_source                       = "Water &\nsanitation",
           Improved_sanitation                         = "Water &\nsanitation",
           No_access_to_a_handwashing_facility         = "Water &\nsanitation",
           Access_to_electricity                       = "Energy",
           Prevalence_of_cooking_with_coalperbiomass   = "Energy",
           Share_Slums                                 = "Shelter",
           Universal_health_coverage                   = "Health",
           Mean_years_of_schooling                     = "Education",
           Mobile_and_landline_telephone_subscriptions = "Communication",
           Internet_users                              = "Communication")
ORD  <- names(LAB)
THLV <- unique(THEME[ORD])

# ---- Long frame -------------------------------------------------------------
pd <- expand.grid(ind = ORD, g = seq_len(G), stringsAsFactors = FALSE)
pd$value <- dev[cbind(pd$g, match(pd$ind, colnames(dev)))]
pd$lab   <- factor(LAB[pd$ind], levels = rev(LAB[ORD]))
pd$theme <- factor(THEME[pd$ind], levels = THLV)
pd$sign  <- factor(ifelse(pd$value >= 0, "Faster than its starting position predicts",
                                         "Slower than its starting position predicts"),
                   levels = c("Faster than its starting position predicts",
                              "Slower than its starting position predicts"))
# ---- Panel headers ----------------------------------------------------------
#  Three fixed lines: regime name / size and level / example countries. The list
#  of countries is packed into exactly two lines that never break inside a
#  country name, and a name that will not fit is dropped rather than clipped.
wrap_list <- function(s, width = 31) {
  it <- strsplit(s, "|", fixed = TRUE)[[1]]
  repeat {
    # every split into two lines that both fit; take the most evenly balanced
    cuts <- seq_len(length(it) - 1)
    a <- vapply(cuts, function(i) paste(it[seq_len(i)],  collapse = ", "), "")
    b <- vapply(cuts, function(i) paste(it[-seq_len(i)], collapse = ", "), "")
    ok <- nchar(a) <= width & nchar(b) <= width
    if (any(ok)) {
      best <- cuts[ok][which.min(abs(nchar(a[ok]) - nchar(b[ok])))]
      return(paste0(a[best], "\n", b[best]))
    }
    if (length(it) <= 2) return(paste0(it[1], "\n", paste(it[-1], collapse = ", ")))
    it <- utils::head(it, -1)                # drop the last example and retry
  }
}
hdr <- sprintf("R%d  %s\nn = %d,  mean level %+.2f\n%s",
               seq_len(G), NAMES, nvec, lvl,
               vapply(examples, wrap_list, "", USE.NAMES = FALSE))
pd$panel <- factor(hdr[pd$g], levels = hdr)

SIGN_COL <- c("Faster than its starting position predicts" = PAL$uk,
              "Slower than its starting position predicts" = PAL$fr)
SIGN_SHP <- c("Faster than its starting position predicts" = 21,   # filled disc
              "Slower than its starting position predicts" = 21)   # ring (surface fill)
SIGN_FIL <- c("Faster than its starting position predicts" = PAL$uk,
              "Slower than its starting position predicts" = PAL$surface)

xr <- range(pd$value); xlim <- c(min(-1.1, xr[1] - .24), max(1.1, xr[2] + .24))

make_row <- function(gs, xaxis) {
  sub <- pd[pd$g %in% gs, ]
  sub$panel <- droplevels(sub$panel)
  ggplot(sub, aes(x = value, y = lab)) +
    geom_vline(xintercept = 0, colour = PAL$axis, linewidth = .45) +
    geom_segment(aes(x = 0, xend = value, yend = lab, colour = sign),
                 linewidth = .62, lineend = "butt", show.legend = FALSE) +
    geom_point(aes(colour = sign, fill = sign, shape = sign),
               size = 1.85, stroke = .62) +
    scale_colour_manual(values = SIGN_COL, drop = FALSE) +
    scale_fill_manual(values   = SIGN_FIL, drop = FALSE) +
    scale_shape_manual(values  = SIGN_SHP, drop = FALSE) +
    scale_x_continuous(limits = xlim, breaks = seq(-2, 3, 1), expand = c(0, 0)) +
    facet_grid(theme ~ panel, scales = "free_y", space = "free_y", switch = "y") +
    labs(x = if (xaxis) "Improvement rate net of starting level (index points per year; 0 = as expected)" else NULL,
         y = NULL) +
    th +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = PAL$grid, linewidth = .28),
          panel.spacing.x = unit(9, "pt"), panel.spacing.y = unit(3.2, "pt"),
          axis.line.x = element_blank(),
          axis.text.y = element_text(colour = PAL$ink2, size = 7.6, hjust = 1),
          axis.text.x = element_text(colour = PAL$muted, size = 7.4),
          axis.title.x = element_text(colour = PAL$ink2, size = 8.2, margin = margin(t = 5)),
          strip.placement = "outside",
          strip.text.x = element_text(colour = PAL$ink, size = 8.1, hjust = 0,
                                      lineheight = 1.22, margin = margin(b = 4, t = 1)),
          strip.text.y.left = element_text(colour = PAL$muted, size = 7.1, angle = 0,
                                           hjust = 0, lineheight = 1.05,
                                           margin = margin(r = 5)),
          legend.position = if (xaxis) "none" else "top",
          legend.justification = "left", legend.location = "plot",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(3, "pt"),
          legend.text = element_text(colour = PAL$ink2, size = 8),
          plot.margin = if (xaxis) margin(9, 8, 2, 2) else margin(2, 8, 2, 2))
}

# 3 + 3 panels; name_regimes() above has already guaranteed G == 6.
rowA <- make_row(1:3, FALSE)
rowB <- make_row(4:6, TRUE)

# Re-flowed at a conservative width so the line never overruns the page in
# either device: cairo_pdf sets slightly wider than the PNG device.
NOTE <- paste(strwrap(paste(
  "Each dot is one of eleven material-provisioning indicators (0-100, higher is better). Improvement rates are ordinary-least-squares slopes on year over 2000-2020,",
  "averaged over the 15 multiply-imputed datasets. Each rate is then residualised on a cubic in the country's own year-2000 level of the same indicator, so zero",
  "means a regime advanced that component exactly as fast as countries starting from the same position typically did; the plotted value is the regime mean of that",
  "residual. This is the space in which the regimes were formed. Table 4 reports the same rates unadjusted, in raw points per year: the two scales disagree in sign",
  "for 11 of 66 cells, so they are not interchangeable. Indicators appear in the same order in every panel, grouped by decent-living-standards domain, and all six",
  "panels share one horizontal scale. Regimes are ordered by mean development level; the example countries are those at the 12th, 38th, 62nd and 88th percentile of",
  "each regime's own level distribution."
), width = 158), collapse = "\n")

fig <- (rowA / rowB) +
  plot_layout(heights = c(1, 1.06)) +
  plot_annotation(
    title = "Six improvement regimes, 2000-2020",
    subtitle = "What each group advanced faster, and slower, than its starting position predicts",
    caption = NOTE,
    theme = th + theme(
      plot.title = element_text(colour = PAL$ink, size = 11, face = "bold", margin = margin(b = 2)),
      plot.subtitle = element_text(colour = PAL$muted, size = 8.6, margin = margin(b = 6)),
      plot.caption = element_text(colour = PAL$muted, size = 6.5, hjust = 0,
                                  lineheight = 1.3, margin = margin(t = 9)),
      plot.margin = margin(9, 11, 8, 9)))
#Every panel carries plot.background through `th`, and so does the annotation
#theme above; the patchwork `&` operator is not needed (and fails with
#ggplot2 >= 4.0 unless patchwork >= 1.3.1 is installed).

W <- 7.5; H <- 7.3; DPI <- 400
f_png <- file.path(OUT, "Fig5_improvement_regimes.png")
f_pdf <- file.path(OUT, "Fig5_improvement_regimes.pdf")
ggsave(f_png, fig, width = W, height = H, dpi = DPI, bg = PAL$surface)
ggsave(f_pdf, fig, width = W, height = H, bg = PAL$surface, device = cairo_pdf)

# ---- Verify the effective resolution from the PNG header (no extra packages) -
png_dim <- function(f) {                      # IHDR: width and height, big-endian
  con <- file(f, "rb"); on.exit(close(con))
  readBin(con, "raw", 16)
  c(w = readBin(con, "integer", 1, 4, endian = "big"),
    h = readBin(con, "integer", 1, 4, endian = "big"))
}
sz <- png_dim(f_png)
cat(sprintf("\nWritten: %s  (%d x %d px)\n", f_png, sz["w"], sz["h"]))
cat(sprintf("Effective resolution: %.0f x %.0f dpi at %.2f x %.2f in  (Elsevier minimum 300 dpi: %s)\n",
            sz["w"]/W, sz["h"]/H, W, H, if (sz["w"]/W >= 300) "PASS" else "FAIL"))
cat(sprintf("Written: %s  (vector, %.1f KB)\n", f_pdf, file.size(f_pdf)/1024))

# ---- Draft caption, in the voice of the Figure 5 caption --------------------
cap <- paste0(
"Figure 5 | Six improvement regimes of material provisioning, 2000-2020.\n\n",
"Each panel is one regime; each dot is one of eleven decent-living indicators. Improvement\n",
"rates are ordinary-least-squares slopes on year over 2000-2020, in index points per year,\n",
"averaged over the fifteen multiply-imputed datasets. Because a country already close to\n",
"the ceiling on a service cannot improve on it as fast as one starting from a low base,\n",
"each rate is residualised on a cubic polynomial in that country's own year-2000 value of\n",
"the same indicator; the plotted quantity is the regime mean of that residual. Zero\n",
"therefore means a regime advanced that component exactly as fast as countries starting\n",
"from the same position typically did, filled dots to the right mark faster-than-expected\n",
"progress, and open dots to the left mark slower-than-expected progress. This is the space\n",
"in which the regimes were formed: they are the groups obtained by Ward minimum-variance\n",
"clustering on Euclidean distances between these standardised residuals, with the number\n",
"of groups chosen by the C-index and the ensemble partition taken as the Mirkin medoid of\n",
"the fifteen imputation-specific partitions. Table 4 reports the same improvement rates\n",
"unadjusted, in raw points per year; the raw and adjusted quantities disagree in sign for\n",
"eleven of the sixty-six regime-indicator cells, so the table and this figure are not\n",
"alternative displays of one quantity and should be read together. The eleven indicators\n",
"appear in the same order in every panel, grouped by decent-living-standards domain, and\n",
"all six panels share one horizontal scale. Regimes are ordered by mean development level,\n",
"the regime mean of the country's average standardised year-2000 value across the eleven\n",
"indicators, shown in each header with the regime's size and four example countries drawn\n",
"from the 12th, 38th, 62nd and 88th percentile of that regime's own level distribution.\n",
"N = ",
length(regime), " countries.\n\n",
paste(sprintf("R%d  %-26s n = %3d   mean level %+5.2f   e.g. %s",
              seq_len(G), NAMES, nvec, lvl,
              gsub("|", ", ", examples, fixed = TRUE)), collapse = "\n"), "\n")
writeLines(cap, file.path(OUT, "Figure2_caption.txt"))
cat(sprintf("Written: %s\n", file.path(OUT, "Figure2_caption.txt")))

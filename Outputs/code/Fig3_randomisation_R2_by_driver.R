## ============================================================
## Variance explained by each candidate driver, against a randomised null (Fig. 3b)
## For every driver: adjusted R2 of an AR(1) GAMM fitted alone, compared with the R2
## obtained when that driver is randomly SHIFTED in time within each collapse phase
## (wrap-around), which keeps its serial structure and destroys only its alignment
## with occupancy. B = 1000 randomisations per driver.
## Run from Outputs/code: Rscript Fig3_randomisation_R2_by_driver.R
## ============================================================
suppressMessages({library(tidyverse); library(mgcv); library(nlme)})
set.seed(42)
B <- 1000
ALIGNED <- "/Users/sarah/Documents/Postdoc/Canada Research/Pacific herring/CLAUDE CODE/PAPER/Aligned_2026_08_06"

dm <- read_csv(file.path(ALIGNED, "Outputs/derived/annual_series.csv"), show_col_types = FALSE) %>%
  filter(Year %in% 1951:2024, !is.na(occupancy), !is.na(prop_old_wt), !is.na(mean_SST),
         !is.na(ONI_ann), !is.na(mean_PDO), !is.na(value_per_ton)) %>%
  arrange(Year) %>%
  mutate(phase = factor(if_else(Year < 1984, "Pre", "Post")))
cat("n =", nrow(dm), "years | B =", B, "\n")

DRIVERS <- tribble(
  ~var,                ~label,                                    ~class,
  "prop_old_wt",       "Proportion of old fish",                   "Selectivity",
  "mean_SST",          "Mean SST",                                 "Climate",
  "ONI_ann",           "Oceanic Nino Index",                       "Climate",
  "mean_PDO",          "Pacific Decadal Oscillation",              "Climate",
  "catch_total_t",     "Total catch",                              "Amount / value / abundance",
  "old_removed_t",     "Tonnage of old fish removed",              "Amount / value / abundance",
  "old_removed_frac",  "Population fraction of old fish removed",  "Amount / value / abundance",
  "value_per_ton",     "Value per tonne (real 2024 CAD)",          "Amount / value / abundance",
  "SSB",               "Spawning stock biomass",                   "Amount / value / abundance")

r2_of <- function(d, vars) {
  rhs <- paste(sprintf("s(%s, k = 5)", vars), collapse = " + ")
  m <- tryCatch(gamm(as.formula(paste("logit_occ ~", rhs)), data = d,
                     correlation = corAR1(form = ~ Year), method = "REML"),
                error = function(e) NULL)
  if (is.null(m)) NA_real_ else summary(m$gam)$r.sq
}
shift_within <- function(x, g) { out <- x
  for (lv in unique(g)) { i <- which(g == lv); n <- length(i); k <- sample.int(n, 1)
    out[i] <- x[i][((seq_len(n) - 1 + k) %% n) + 1] }; out }

one_driver <- function(var, label, class) {
  obs  <- r2_of(dm, var)
  null <- map_dbl(seq_len(B), function(i) {
    d <- dm; d[[var]] <- shift_within(dm[[var]], dm$phase); r2_of(d, var) })
  null <- null[!is.na(null)]
  tibble(Driver = label, Class = class,
         `adj R2` = round(obs, 3),
         `null R2 median` = round(median(null), 3),
         `null R2 95th pct` = round(quantile(null, 0.95), 3),
         `Empirical P` = signif((1 + sum(null >= obs)) / (1 + length(null)), 3),
         B = length(null))
}
res <- pmap_dfr(list(DRIVERS$var, DRIVERS$label, DRIVERS$class), one_driver)

## best two-term model, null shifts both predictors independently
obs_best <- r2_of(dm, c("prop_old_wt", "mean_SST"))
null_best <- map_dbl(seq_len(B), function(i) {
  d <- dm
  d$prop_old_wt <- shift_within(dm$prop_old_wt, dm$phase)
  d$mean_SST    <- shift_within(dm$mean_SST,    dm$phase)
  r2_of(d, c("prop_old_wt", "mean_SST")) })
null_best <- null_best[!is.na(null_best)]
res <- bind_rows(res,
  tibble(Driver = "Old fish + SST (best model)", Class = "Best model",
         `adj R2` = round(obs_best, 3),
         `null R2 median` = round(median(null_best), 3),
         `null R2 95th pct` = round(quantile(null_best, 0.95), 3),
         `Empirical P` = signif((1 + sum(null_best >= obs_best)) / (1 + length(null_best)), 3),
         B = length(null_best))) %>%
  arrange(desc(`adj R2`))

print(as.data.frame(res), row.names = FALSE)
write_csv(res, file.path(ALIGNED, "Outputs/tables/randomisation_R2_by_driver.csv"))
cat("\nsaved: Outputs/tables/randomisation_R2_by_driver.csv\n")

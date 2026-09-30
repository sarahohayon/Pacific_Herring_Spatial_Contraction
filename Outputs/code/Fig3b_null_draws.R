## =============================================================================
## Null distributions of variance explained, for Figure 3b and Table S8b
## Each driver is fitted alone (AR(1) GAMM); its null comes from 1,000 shifts of the driver
## in time within collapse phase. Every draw is kept so the null can be plotted.
## Run from Outputs/code:  Rscript Fig3b_null_draws.R, then Fig3b_null_draws_extra.R
## =============================================================================
suppressMessages({library(tidyverse); library(mgcv); library(nlme)})
set.seed(42)
B <- 1000
ALIGNED <- normalizePath("../..")   # repository root; run from Outputs/code

dm <- read_csv(file.path(ALIGNED, "Outputs/derived/annual_series.csv"), show_col_types = FALSE) %>%
  filter(Year %in% 1951:2024, !is.na(occupancy), !is.na(prop_old_wt), !is.na(mean_SST),
         !is.na(ONI_ann), !is.na(mean_PDO)) %>%
  arrange(Year) %>% mutate(phase = factor(if_else(Year < 1984, "Pre", "Post")))

DRIVERS <- tribble(
  ~var,            ~label,                   ~class,
  "prop_old_wt",   "Old fish removed",        "Selectivity",
  "catch_total_t", "Total catch",             "Other",
  "mean_SST",      "Mean SST",                "Other",
  "SSB",           "Spawning stock biomass",  "Other")

r2_of <- function(d, v) {
  m <- tryCatch(gamm(as.formula(sprintf("logit_occ ~ s(%s, k = 5)", v)), data = d,
                     correlation = corAR1(form = ~ Year), method = "REML"),
                error = function(e) NULL)
  if (is.null(m)) NA_real_ else summary(m$gam)$r.sq
}
shift_within <- function(x, g) { out <- x
  for (lv in unique(g)) { i <- which(g == lv); n <- length(i); k <- sample.int(n, 1)
    out[i] <- x[i][((seq_len(n) - 1 + k) %% n) + 1] }; out }

out <- pmap_dfr(list(DRIVERS$var, DRIVERS$label, DRIVERS$class), function(v, lab, cls) {
  obs  <- r2_of(dm, v)
  null <- map_dbl(seq_len(B), function(i) { d <- dm; d[[v]] <- shift_within(dm[[v]], dm$phase); r2_of(d, v) })
  null <- null[!is.na(null)]
  cat(sprintf("%-24s observed %.3f | null median %.3f | P = %.3f\n", lab, obs, median(null),
              (1 + sum(null >= obs)) / (1 + length(null))))
  tibble(driver = lab, class = cls, observed = obs, null_r2 = null,
         p_emp = (1 + sum(null >= obs)) / (1 + length(null)))
})
write_csv(out, file.path(ALIGNED, "Outputs/derived/fig3b_null_draws.csv"))
cat("\nsaved: Outputs/derived/fig3b_null_draws.csv (", nrow(out), "rows )\n")

## =============================================================================
## Null draws for the remaining drivers of Figure 3b and Table S8b (ONI, PDO and the two
## removal metrics). Same scheme as Fig3b_null_draws.R; the draws are appended to its file.
## Run after Fig3b_null_draws.R.
## =============================================================================
suppressMessages({library(tidyverse); library(mgcv); library(nlme)})
B <- 1000
ALIGNED <- normalizePath("../..")   # repository root; run from Outputs/code
f_out <- file.path(ALIGNED, "Outputs/derived/fig3b_null_draws.csv")
have <- read_csv(f_out, show_col_types = FALSE)
dm <- read_csv(file.path(ALIGNED, "Outputs/derived/annual_series.csv"), show_col_types = FALSE) %>%
  filter(Year %in% 1951:2024, !is.na(occupancy), !is.na(prop_old_wt), !is.na(mean_SST),
         !is.na(ONI_ann), !is.na(mean_PDO)) %>%
  arrange(Year) %>% mutate(phase = factor(if_else(Year < 1984, "Pre", "Post")))
EXTRA <- tribble(
  ~var,               ~label,                                     ~class,
  "ONI_ann",          "Oceanic Nino Index",                        "Climate",
  "mean_PDO",         "Pacific Decadal Oscillation",               "Climate",
  "old_removed_t",    "Tonnage of old fish removed",               "Fishery",
  "old_removed_frac", "Population fraction of old fish removed",   "Fishery")
r2_of <- function(d, v) {
  m <- tryCatch(gamm(as.formula(sprintf("logit_occ ~ s(%s, k = 5)", v)), data = d,
                     correlation = corAR1(form = ~ Year), method = "REML"), error = function(e) NULL)
  if (is.null(m)) NA_real_ else summary(m$gam)$r.sq }
shift_within <- function(x, g) { out <- x
  for (lv in unique(g)) { i <- which(g == lv); n <- length(i); k <- sample.int(n, 1)
    out[i] <- x[i][((seq_len(n) - 1 + k) %% n) + 1] }; out }
new <- pmap_dfr(list(EXTRA$var, EXTRA$label, EXTRA$class), function(v, lab, cls) {
  set.seed(1000 + match(v, EXTRA$var))          # per-driver seed, reproducible on its own
  obs  <- r2_of(dm, v)
  null <- map_dbl(seq_len(B), function(i) { d <- dm; d[[v]] <- shift_within(dm[[v]], dm$phase); r2_of(d, v) })
  null <- null[!is.na(null)]
  p <- (1 + sum(null >= obs)) / (1 + length(null))
  cat(sprintf("%-42s observed %7.3f | null median %6.3f | P = %.3f\n", lab, obs, median(null), p))
  tibble(driver = lab, class = cls, observed = obs, null_r2 = null, p_emp = p) })
write_csv(bind_rows(have, new), f_out)
cat("\nfig3b_null_draws.csv now has", n_distinct(bind_rows(have, new)$driver), "drivers\n")

## =============================================================================
## Does SST add more to the older-fish model than chance would? (Table S8c)
## SST is shifted in time within collapse phase while older fish is left unchanged; each time
## older fish + SST is refitted and the gain in adjusted R2 is recorded (1,000 times).
## Run from Outputs/code:  Rscript Fig3b_sst_increment_test.R
## =============================================================================
suppressMessages({library(tidyverse); library(mgcv); library(nlme)})
set.seed(42); B <- 1000
ALIGNED <- normalizePath("../..")   # repository root; run from Outputs/code
dm <- read_csv(file.path(ALIGNED, "Outputs/derived/annual_series.csv"), show_col_types = FALSE) %>%
  filter(Year %in% 1951:2024, !is.na(occupancy), !is.na(prop_old_wt), !is.na(mean_SST)) %>%
  arrange(Year) %>% mutate(phase = factor(if_else(Year < 1984, "Pre", "Post")))
r2 <- function(d, rhs) {
  m <- tryCatch(gamm(as.formula(paste("logit_occ ~", rhs)), data = d,
                     correlation = corAR1(form = ~ Year), method = "REML"), error = function(e) NULL)
  if (is.null(m)) NA_real_ else summary(m$gam)$r.sq
}
shift_within <- function(x, g) { out <- x
  for (lv in unique(g)) { i <- which(g == lv); n <- length(i); k <- sample.int(n, 1)
    out[i] <- x[i][((seq_len(n) - 1 + k) %% n) + 1] }; out }

r2_old  <- r2(dm, "s(prop_old_wt, k = 5)")
r2_both <- r2(dm, "s(prop_old_wt, k = 5) + s(mean_SST, k = 5)")
obs_inc <- r2_both - r2_old
cat(sprintf("observed: prop_old %.3f | prop_old + SST %.3f | increment %.3f\n", r2_old, r2_both, obs_inc))

null_inc <- map_dbl(seq_len(B), function(i) {
  d <- dm; d$mean_SST <- shift_within(dm$mean_SST, dm$phase)
  r2(d, "s(prop_old_wt, k = 5) + s(mean_SST, k = 5)") - r2_old })
null_inc <- null_inc[!is.na(null_inc)]
p <- (1 + sum(null_inc >= obs_inc)) / (1 + length(null_inc))
cat(sprintf("null increment: median %.3f, 95th pct %.3f (B = %d)\n",
            median(null_inc), quantile(null_inc, 0.95), length(null_inc)))
cat(sprintf("SST increment beyond old fish: P = %.3f\n", p))
write_csv(tibble(r2_old, r2_both, obs_inc, null_med = median(null_inc),
                 null_p95 = quantile(null_inc, 0.95), p_emp = p, B = length(null_inc)),
          file.path(ALIGNED, "Outputs/tables/randomisation_SST_increment.csv"))

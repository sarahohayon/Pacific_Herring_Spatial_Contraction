## =============================================================================
## Randomization test of the slope of occupancy on the proportion of older fish (Table S8a)
## Is the observed slope stronger than when the two series are paired at random?
##   shift    moves the older-fish series in time with wrap-around, within era or phase.
##            It keeps the series' own autocorrelation, so this is the test reported.
##   shuffle  free permutation; it ignores autocorrelation and is shown for comparison.
## Models: GLS with AR(1) errors, without and with collapse phase as a covariate.
## Also fits the gear-era model reported in the Results.
## Run from Outputs/code:  Rscript Fig3_randomisation_test.R
## =============================================================================
suppressMessages({library(tidyverse); library(nlme)})
set.seed(42)
B <- 1000
ALIGNED <- normalizePath("../..")   # repository root; run from Outputs/code

dm <- read_csv(file.path(ALIGNED, "Outputs/derived/annual_series.csv"), show_col_types = FALSE) %>%
  filter(Year %in% 1951:2024, !is.na(occupancy), !is.na(prop_old_wt), !is.na(mean_SST)) %>%
  arrange(Year) %>%
  mutate(era   = factor(if_else(Year <= 1965, "Reduction seine", "Roe gillnet")),
         phase = factor(if_else(Year < 1984, "Pre-collapse", "Post-collapse"),
                        levels = c("Pre-collapse", "Post-collapse")))
cat(sprintf("n = %d years | %d randomisations per test\n", nrow(dm), B))

## gear era as covariate (Methods; Results "It held likewise after accounting for gear era")
m_era <- gls(logit_occ ~ era + prop_old_wt, data = dm, correlation = corAR1(form = ~ Year), method = "ML")
ci_era <- intervals(m_era, which = "coef")$coef["prop_old_wt", ]
cat(sprintf("Gear-era model: slope = %.2f (95%% CI %.2f to %.2f), P = %.3g\n",
            ci_era[2], ci_era[1], ci_era[3], summary(m_era)$tTable["prop_old_wt", 4]))

stat_of <- function(x, model) {
  d <- dm; d$prop_old_wt <- x
  f <- if (model == "basic") logit_occ ~ prop_old_wt else logit_occ ~ phase + prop_old_wt
  m <- tryCatch(gls(f, data = d, correlation = corAR1(form = ~ Year), method = "ML"),
                error = function(e) NULL)
  if (is.null(m)) return(c(slope = NA_real_, t = NA_real_))
  tt <- summary(m)$tTable
  c(slope = unname(tt["prop_old_wt", 1]), t = unname(tt["prop_old_wt", 3]))
}
shuffle_within <- function(x, g) { out <- x
  for (lv in unique(g)) { i <- which(g == lv); out[i] <- sample(x[i]) }; out }
shift_within <- function(x, g) { out <- x
  for (lv in unique(g)) { i <- which(g == lv); n <- length(i); k <- sample.int(n, 1)
    out[i] <- x[i][((seq_len(n) - 1 + k) %% n) + 1] }; out }

run_test <- function(model, group, scheme) {
  g   <- dm[[group]]
  obs <- stat_of(dm$prop_old_wt, model)
  fun <- if (scheme == "shuffle") shuffle_within else shift_within
  null <- map_dfr(seq_len(B), ~ as_tibble_row(stat_of(fun(dm$prop_old_wt, g), model))) %>% drop_na()
  p_emp <- (1 + sum(null$t <= obs[["t"]])) / (1 + nrow(null))      # one-sided: more negative than observed
  tibble(Model = model, `Randomised within` = group, Scheme = scheme,
         `Observed slope` = round(obs[["slope"]], 2), `Observed t` = round(obs[["t"]], 2),
         `Null slope mean` = round(mean(null$slope), 2), `Null slope SD` = round(sd(null$slope), 2),
         `Null t 5th pct` = round(quantile(null$t, 0.05), 2),
         `n randomisations` = nrow(null),
         `Empirical P` = signif(p_emp, 3))
}

res <- bind_rows(
  run_test("basic", "era",   "shuffle"),
  run_test("basic", "era",   "shift"),
  run_test("phase", "phase", "shuffle"),
  run_test("phase", "phase", "shift"))
print(as.data.frame(res), row.names = FALSE)
write_csv(res, file.path(ALIGNED, "Outputs/tables/randomisation_test_occupancy_prop_old.csv"))
cat("\nsaved: Outputs/tables/randomisation_test_occupancy_prop_old.csv\n")

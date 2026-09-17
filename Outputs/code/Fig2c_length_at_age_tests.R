## ============================================================
## Figure 2c - statistical tests for length-at-age (Results 2.2, Fig. 2c caption)
## Test-fishery fish, statistical areas 13, 14, 15, 17, 18, ages 2-8, from Data/Biosample_sog.csv
##  1  1980 -> 1983 drop: linear mixed model, length ~ age class * year + (1 | sample);
##     LRT for the drop and for the age x year interaction (is the drop larger in older fish?);
##     robustness: statistical area and month added as fixed effects
##  2  post-collapse decline: annual mean length-at-age 1984-2024 (ages 2-7), GLS with AR(1) within age class
## Run from Outputs/code: Rscript Fig2c_length_at_age_tests.R
## ============================================================
suppressMessages({library(tidyverse); library(lme4); library(nlme)})
BASE <- "/Users/sarah/Documents/Postdoc/Canada Research/Pacific herring/CLAUDE CODE"
b <- read_csv(file.path(BASE, "Data/Biosample_sog.csv"), show_col_types = FALSE, guess_max = 100000) %>%
  mutate(Year = as.integer(Year), age = suppressWarnings(as.numeric(age)),
         length_mm = suppressWarnings(as.numeric(length_mm)),
         stat_area = suppressWarnings(as.numeric(stat_area)),
         sample_id = paste(Year, season, sample_number)) %>%
  filter(stock_assessment_region == "Strait of Georgia", stat_area %in% c(13,14,15,17,18),
         source == "Test Fishery", !is.na(age), !is.na(length_mm), age %in% 2:8)

## ---- Test 1: the 1980 -> 1983 drop, and whether it was larger in older fish ----
d1 <- b %>% filter(Year %in% c(1980, 1983)) %>%
  mutate(yr = factor(Year), age_f = factor(age), len_cm = length_mm / 10)
cat("Fish:", nrow(d1), "| samples:", n_distinct(d1$sample_id), "\n")
print(d1 %>% count(Year, age) %>% pivot_wider(names_from = age, values_from = n))
m_int  <- lmer(len_cm ~ age_f * yr + (1 | sample_id), data = d1, REML = FALSE)
m_add  <- lmer(len_cm ~ age_f + yr + (1 | sample_id), data = d1, REML = FALSE)
m_none <- lmer(len_cm ~ age_f + (1 | sample_id), data = d1, REML = FALSE)
lr_drop <- anova(m_none, m_add); lr_int <- anova(m_add, m_int)
cat(sprintf("\nDrop 1980->1983 (all ages): LRT chi2 = %.1f, df = %d, P = %.2g\n", lr_drop$Chisq[2], lr_drop$Df[2], lr_drop$`Pr(>Chisq)`[2]))
cat(sprintf("Drop differs among ages (age x year): LRT chi2 = %.1f, df = %d, P = %.2g\n", lr_int$Chisq[2], lr_int$Df[2], lr_int$`Pr(>Chisq)`[2]))
## linear version: does the drop grow with age?
d1 <- d1 %>% mutate(post = as.numeric(yr == "1983"), age0 = age - 2)
m_lin <- lmer(len_cm ~ age_f + post + post:age0 + (1 | sample_id), data = d1, REML = FALSE)
cf <- summary(m_lin)$coefficients
cat(sprintf("Drop at age 2 = %.2f cm; extra drop per year of age = %.2f cm (SE %.2f, t = %.1f)\n",
            cf["post", 1], cf["post:age0", 1], cf["post:age0", 2], cf["post:age0", 3]))
cat(sprintf("LRT for the age-dependence (linear): P = %.2g\n", anova(m_add, m_lin)$`Pr(>Chisq)`[2]))
## per-age drop with 95% CI from the interaction model
nd <- expand_grid(age_f = levels(d1$age_f), yr = levels(d1$yr))
X <- model.matrix(~ age_f * yr, nd %>% mutate(age_f = factor(age_f, levels(d1$age_f)), yr = factor(yr, levels(d1$yr))))
beta <- fixef(m_int); V <- as.matrix(vcov(m_int))
per_age <- map_dfr(levels(d1$age_f), function(a) {
  L <- X[nd$age_f == a & nd$yr == "1983", ] - X[nd$age_f == a & nd$yr == "1980", ]
  est <- sum(L * beta); se <- sqrt(drop(t(L) %*% V %*% L))
  tibble(age = a, drop_cm = round(est, 2), lo = round(est - 1.96*se, 2), hi = round(est + 1.96*se, 2))
})
print(per_age)

## robustness: were 1980 and 1983 sampled in different areas or months?
print(d1 %>% count(Year, stat_area) %>% group_by(Year) %>% mutate(pct = round(100*n/sum(n))) %>% select(-n) %>%
        pivot_wider(names_from = stat_area, values_from = pct, values_fill = 0))
print(d1 %>% count(Year, month) %>% group_by(Year) %>% mutate(pct = round(100*n/sum(n))) %>% select(-n) %>%
        pivot_wider(names_from = month, values_from = pct, values_fill = 0))
d1r <- d1 %>% mutate(area_f = factor(stat_area), month_f = factor(month))
r_add <- lmer(len_cm ~ age_f + yr + area_f + month_f + (1 | sample_id), data = d1r, REML = FALSE)
r_int <- lmer(len_cm ~ age_f * yr + area_f + month_f + (1 | sample_id), data = d1r, REML = FALSE)
r_lin <- lmer(len_cm ~ age_f + post + post:age0 + area_f + month_f + (1 | sample_id), data = d1r, REML = FALSE)
lr_r <- anova(r_add, r_int); cr <- summary(r_lin)$coefficients
cat(sprintf("Robustness (+ area + month): age x year LRT chi2 = %.1f, P = %.2g | extra drop per year of age = %.2f cm (SE %.2f)\n",
            lr_r$Chisq[2], lr_r$`Pr(>Chisq)`[2], cr["post:age0", 1], cr["post:age0", 2]))


## ---- Test 2: progressive decline after the collapse (annual means 1984-2024, AR(1)) ----
d2 <- b %>% filter(Year >= 1984) %>% group_by(Year, age) %>%
  summarise(len_cm = mean(length_mm)/10, n = n(), .groups = "drop") %>% filter(n >= 30, age <= 7) %>%
  mutate(age_f = factor(age))
g_common <- gls(len_cm ~ age_f + Year, data = d2, correlation = corAR1(form = ~ Year | age_f), method = "ML")
g_byage  <- gls(len_cm ~ age_f + age_f:Year, data = d2, correlation = corAR1(form = ~ Year | age_f), method = "ML")
tt <- summary(g_common)$tTable
cat(sprintf("\nPost-1984 decline, common slope: %.3f cm/yr = %.2f cm/decade (SE %.3f, P = %.2g), AR1 phi = %.2f\n",
            tt["Year", 1], 10*tt["Year", 1], tt["Year", 2], tt["Year", 4],
            coef(g_common$modelStruct$corStruct, unconstrained = FALSE)))
tb <- summary(g_byage)$tTable; tb <- tb[grepl(":Year", rownames(tb)), ]
print(round(cbind(cm_per_decade = 10*tb[, 1], P = tb[, 4]), 3))
cat(sprintf("Slopes differ among ages: LRT P = %.2g\n", anova(g_common, g_byage)$`p-value`[2]))

## ---- Rate comparison: 1980-83 drop per year vs post-1984 rate, per age ----
rate <- per_age %>% mutate(cm_per_yr_1980_83 = round(drop_cm/3, 2))
print(rate)

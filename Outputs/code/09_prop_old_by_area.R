## =============================================================================
## 09_prop_old_by_area.R
## Proportion of older fish within single statistical areas.
##
## 1. Does the 1972-1980 rise in older fish in the roe-gillnet catch hold within single
##    areas, or does it come from samples shifting between areas?
## 2. Did Area 14 keep more older fish than the peripheral areas in the test fishery
##    (Results, "The core also retained its older fish")?
##
## Proportions are computed per sample (at least MIN_FISH aged fish), then averaged per year.
## Run from Outputs/code:  Rscript 09_prop_old_by_area.R
## =============================================================================

source("00_setup.R")
suppressPackageStartupMessages(library(lme4))
## showtext renders at 96 dpi unless told otherwise, which shrinks text in ggsave(dpi = 300)
if (nzchar(HERRING_FONT)) showtext::showtext_opts(dpi = 300)

bio <- read_csv(F_BIO, show_col_types = FALSE, name_repair = "minimal",
                guess_max = 100000)
names(bio) <- BIO_COLS

bio <- bio %>%
  mutate(age       = suppressWarnings(as.numeric(age)),
         weight_g  = suppressWarnings(as.numeric(weight_g)),
         stat_area = suppressWarnings(as.numeric(stat_area)),
         Year      = recode_two_digit_year(year),
         gear_grp  = gear_group(str_trim(gear)),
         src = case_when(
           source == "Reduction Fishery"                    ~ "Reduction",
           source == "Test Fishery"                         ~ "Test fishery",
           source == "Roe Fishery" & gear_grp == "Gillnet"  ~ "Roe gillnet",
           source == "Roe Fishery" & gear_grp == "Seine"    ~ "Roe seine",
           TRUE ~ NA_character_)) %>%
  filter(stat_area %in% BIO_AREAS, !is.na(src), !is.na(age), age >= 1, !is.na(weight_g),
         keep_tf(source, Year, str_trim(gear)))

## ---- per-sample, then per area-year ----------------------------------------
samples <- bio %>%
  group_by(src, stat_area, Year, season, sample_number) %>%
  summarise(n_fish = n(),
            po_wt  = sum(weight_g[age >= OLD_AGE]) / sum(weight_g),
            po_n   = mean(age >= OLD_AGE), .groups = "drop") %>%
  filter(n_fish >= MIN_FISH)

annual_area <- samples %>%
  group_by(src, stat_area, Year) %>%
  summarise(n_samples = n(), n_aged = sum(n_fish),
            prop_old_wt = mean(po_wt), prop_old_n = mean(po_n), .groups = "drop")

annual_pooled <- samples %>%
  group_by(src, Year) %>%
  summarise(n_samples = n(), prop_old_wt = mean(po_wt), prop_old_n = mean(po_n),
            .groups = "drop") %>%
  mutate(stat_area = "Pooled (5 areas)")

save_tab(annual_area, "prop_old_by_area_annual")

## ---- Q1: roe-gillnet rise, 1972-1980, by area --------------------------------
cat("\n================ Q1  Roe gillnet: proportion old (by weight) by area ================\n")
gn <- annual_area %>% filter(src == "Roe gillnet")

cat("\nYear-by-area values 1972-1990 (prop_old_wt [n samples]); '.' = no sample\n")
gn %>% filter(Year <= 1990) %>%
  mutate(cell = sprintf("%.2f [%d]", prop_old_wt, n_samples)) %>%
  select(Year, stat_area, cell) %>%
  pivot_wider(names_from = stat_area, values_from = cell, values_fill = ".",
              names_prefix = "Area ") %>%
  left_join(annual_pooled %>% filter(src == "Roe gillnet") %>%
              transmute(Year, Pooled = sprintf("%.2f", prop_old_wt)), by = "Year") %>%
  arrange(Year) %>% print(n = 30, width = 200)

rise <- gn %>% filter(Year >= 1972, Year <= 1980) %>%
  group_by(stat_area) %>%
  summarise(n_years   = n(),
            years     = paste(sort(Year), collapse = ","),
            mean_7274 = mean(prop_old_wt[Year <= 1974]),
            mean_7880 = mean(prop_old_wt[Year >= 1978]),
            slope_yr  = if (n() >= 3) coef(lm(prop_old_wt ~ Year))[2] else NA_real_,
            P_slope   = if (n() >= 4) summary(lm(prop_old_wt ~ Year))$coefficients[2, 4] else NA_real_,
            .groups = "drop")
cat("\nRise 1972-1980 within each area (weight-based)\n")
print(rise %>% mutate(across(where(is.double), ~ round(.x, 3))), width = 200)
save_tab(rise, "prop_old_gillnet_rise_1972_1980_by_area")

## composition-controlled trend: year slope with area fixed effects
m_fe  <- lm(prop_old_wt ~ Year + factor(stat_area), data = gn %>% filter(Year >= 1972, Year <= 1980))
m_raw <- lm(prop_old_wt ~ Year, data = annual_pooled %>% filter(src == "Roe gillnet", Year >= 1972, Year <= 1980))
cat(sprintf("\n1972-1980 slope, pooled series:              %+.3f per yr (P %s)\n",
            coef(m_raw)[2], fmt_p(summary(m_raw)$coefficients[2, 4])))
cat(sprintf("1972-1980 slope, within areas (area FE):     %+.3f per yr (P %s)\n",
            coef(m_fe)[2], fmt_p(summary(m_fe)$coefficients[2, 4])))

## does the post-1981 pooled series simply equal Area 14?
post <- annual_pooled %>% filter(src == "Roe gillnet", Year >= 1981) %>%
  select(Year, pooled = prop_old_wt) %>%
  left_join(gn %>% filter(stat_area == 14) %>% select(Year, a14 = prop_old_wt), by = "Year") %>%
  left_join(gn %>% filter(stat_area == 17) %>% select(Year, a17 = prop_old_wt), by = "Year")
cat(sprintf("\n1981-2024: pooled vs Area-14-only r = %.2f (n = %d yrs); pooled vs Area-17-only r = %.2f (n = %d yrs)\n",
            cor(post$pooled, post$a14, use = "complete.obs"), sum(complete.cases(post[, c("pooled","a14")])),
            cor(post$pooled, post$a17, use = "complete.obs"), sum(complete.cases(post[, c("pooled","a17")]))))

## ---- Q2: test fishery, did Area 14 retain old fish? --------------------------
cat("\n================ Q2  Test fishery: proportion old (by count) by area ================\n")
tf <- annual_area %>% filter(src == "Test fishery", Year >= 1975) %>%
  mutate(phase = factor(PHASE3(Year), levels = c("Pre-collapse","Collapse","Post-collapse")),
         group = ifelse(stat_area == 14, "Refuge (14)", "Periphery"))

tf_phase <- tf %>% group_by(stat_area, phase) %>%
  summarise(n_years = n(), mean = mean(prop_old_n), sd = sd(prop_old_n), .groups = "drop") %>%
  mutate(cell = sprintf("%.2f ± %.2f (%d)", mean, sd, n_years))
cat("\nMean ± SD (n years) by area and phase\n")
tf_phase %>% select(stat_area, phase, cell) %>%
  pivot_wider(names_from = phase, values_from = cell) %>% print(width = 200)
save_tab(tf_phase %>% select(-cell), "prop_old_testfishery_by_area_phase")

## change from pre-collapse, per area
chg <- tf_phase %>% select(stat_area, phase, mean) %>%
  pivot_wider(names_from = phase, values_from = mean) %>%
  mutate(`Post vs pre (%)` = round(100 * (`Post-collapse` / `Pre-collapse` - 1)))
cat("\nPost-collapse relative to pre-collapse, by area\n"); print(chg %>% mutate(across(where(is.double), ~ round(.x, 3))))

## rank of Area 14 among areas sampled in the same year
rk <- tf %>% group_by(Year) %>% filter(n() >= 3, any(stat_area == 14)) %>%
  mutate(rank14 = rank(-prop_old_n)[stat_area == 14][1], n_areas = n()) %>%
  summarise(phase = first(phase), rank14 = first(rank14), n_areas = first(n_areas), .groups = "drop")
cat("\nYears in which Area 14 had the HIGHEST prop_old among >= 3 sampled areas:\n")
print(rk %>% group_by(phase) %>% summarise(years = n(), top = sum(rank14 == 1),
                                            bottom = sum(rank14 == n_areas), mean_rank = round(mean(rank14), 2)))

## area x phase model on annual area values
m_tf <- lm(prop_old_n ~ group * phase, data = tf)
cat("\nlm(prop_old_n ~ group * phase), annual area values\n")
print(anova(m_tf))
emm <- tf %>% group_by(group, phase) %>% summarise(mean = round(mean(prop_old_n), 3), n = n(), .groups = "drop")
print(emm %>% pivot_wider(names_from = group, values_from = c(mean, n)))

## ---- Q2b: sample-level GLMM (the test to report) -----------------------------
## The annual-means model above has little power (8 pre-collapse Area 14 years), so here
## each test-fishery sample is an observation, with a year random effect for the shared
## annual age structure and an observation-level random effect for overdispersion.
tfs <- bio %>%
  filter(src == "Test fishery", Year >= 1975) %>%
  group_by(stat_area, Year, season, sample_number) %>%
  summarise(n = n(), n_old = sum(age >= OLD_AGE), .groups = "drop") %>%
  filter(n >= MIN_FISH) %>%
  mutate(group = factor(ifelse(stat_area == 14, "Refuge14", "Periphery"),
                        levels = c("Periphery", "Refuge14")),
         phase = factor(ifelse(Year <= 1983, "Pre", "Post"), levels = c("Pre", "Post")),
         sid = factor(row_number()), yr = factor(Year), area = factor(stat_area))

ctl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
g1 <- glmer(cbind(n_old, n - n_old) ~ group * phase + (1 | yr) + (1 | sid),
            family = binomial, data = tfs, control = ctl)
g0 <- update(g1, . ~ group + phase + (1 | yr) + (1 | sid))
g2 <- update(g1, . ~ . + (1 | area))
lrt <- anova(g0, g1)

contrast_or <- function(m, post) {
  b <- fixef(m); V <- vcov(m); L <- c(0, 1, 0, as.numeric(post))
  est <- sum(L * b); se <- sqrt(drop(t(L) %*% V %*% L))
  tibble(OR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se),
         P = 2 * pnorm(-abs(est / se)))
}
glmm_tab <- bind_rows(
  contrast_or(g1, FALSE) %>% mutate(Contrast = "Area 14 vs periphery, pre-collapse (<=1983)"),
  contrast_or(g1, TRUE)  %>% mutate(Contrast = "Area 14 vs periphery, post-collapse (>=1984)"),
  contrast_or(g2, TRUE)  %>% mutate(Contrast = "Post-collapse, + area random effect"))
for (a in c(13, 15, 17, 18)) {
  d  <- tfs %>% filter(phase == "Post", stat_area %in% c(14, a))
  m  <- glmer(cbind(n_old, n - n_old) ~ group + (1 | yr) + (1 | sid),
              family = binomial, data = d, control = ctl)
  co <- summary(m)$coefficients["groupRefuge14", ]
  glmm_tab <- bind_rows(glmm_tab, tibble(
    Contrast = sprintf("Post-collapse, Area 14 vs Area %d (%d vs %d samples)", a,
                       sum(d$stat_area == 14), sum(d$stat_area == a)),
    OR = exp(co[1]), lo = exp(co[1] - 1.96 * co[2]), hi = exp(co[1] + 1.96 * co[2]), P = co[4]))
}
glmm_tab <- glmm_tab %>% select(Contrast, everything())

cat(sprintf("\n================ Q2b  Sample-level GLMM (%d samples, %s fish) ================\n",
            nrow(tfs), format(sum(tfs$n), big.mark = ",")))
cat(sprintf("Area type x phase interaction: LRT chi2 = %.2f, df = 1, P %s\n",
            lrt$Chisq[2], fmt_p(lrt$`Pr(>Chisq)`[2])))
print(glmm_tab %>% mutate(across(c(OR, lo, hi), ~ round(.x, 2)), P = signif(P, 2)), width = 200)

paired <- tfs %>% filter(phase == "Post") %>%
  group_by(Year, group) %>% summarise(p = sum(n_old) / sum(n), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = p) %>% filter(!is.na(Refuge14), !is.na(Periphery))
cat(sprintf("Post-collapse years with both sampled: %d; Area 14 higher in %d; mean difference %+.3f (paired t P %s)\n",
            nrow(paired), sum(paired$Refuge14 > paired$Periphery),
            mean(paired$Refuge14 - paired$Periphery),
            fmt_p(t.test(paired$Refuge14, paired$Periphery, paired = TRUE)$p.value)))
save_tab(glmm_tab, "prop_old_testfishery_area14_glmm")

## ---- figure ------------------------------------------------------------------
plot_df <- annual_area %>%
  filter(src %in% c("Roe gillnet", "Test fishery"), Year >= 1972) %>%
  mutate(area = factor(stat_area),
         value = ifelse(src == "Roe gillnet", prop_old_wt, prop_old_n),
         panel = ifelse(src == "Roe gillnet",
                        "a  Roe-gillnet catch (proportion age ≥5 by weight)",
                        "b  Test fishery (proportion age ≥5 by count)"))

p <- ggplot(plot_df, aes(Year, value, colour = area)) +
  band_postcoll() +
  geom_point(aes(size = n_samples), alpha = 0.45) +
  geom_smooth(method = "loess", span = 0.5, se = FALSE, linewidth = 1,
              data = plot_df %>% group_by(panel, area) %>% filter(n() >= 8)) +
  geom_vline(xintercept = 1981, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = area_cols, name = "Statistical area") +
  scale_size_area(max_size = 4, name = "Samples") +
  scale_x_continuous(limits = c(1972, 2025), breaks = seq(1975, 2025, 10)) +
  labs(x = "Year", y = "Proportion of old fish",
       caption = "Dashed line: 1981 Area Licensing. Grey: post-collapse (≥1984). Lines: loess for areas with ≥8 sampled years.") +
  theme_herring()
save_fig(p, "prop_old_by_area", width = 9, height = 8)
cat("\nWrote Outputs/figures/prop_old_by_area.png and tables prop_old_*_by_area*.csv\n")

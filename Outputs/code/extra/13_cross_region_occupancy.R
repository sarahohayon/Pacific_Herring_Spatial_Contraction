## ---------------------------------------------------------------------------
## Does the Strait of Georgia pattern hold elsewhere in British Columbia?
##
## The same two series, built the same way, for every herring stock assessment
## region: spawning occupancy from the spawn index, and the proportion of old
## fish in the most size-selective gear of each era from the biosample database.
##
## Rules are the paper's, not the earlier draft's:
##   presence  a surveyed spawn record, including records without a biomass
##             estimate, excluding those DFO classes as Incomplete (keep_surveyed)
##   sections  historically used = spawn recorded in at least MIN_ACTIVE years
##   samples   at least MIN_FISH aged fish, proportions averaged across samples
##   models    AR(1) GAMM of logit(occupancy) on a smooth of the selectivity index
##
## The regional biosample files live outside this repository, in the DFO
## biosample database folder; BIO_DB below points at them.
## ---------------------------------------------------------------------------
source("00_setup.R")
suppressMessages({library(data.table); library(mgcv); library(nlme); library(strucchange)
                  library(tidyr); library(patchwork); library(scales)})

BIO_DB <- "/Users/sarah/Documents/Postdoc/Canada Research/Pacific herring/CLAUDE CODE/Data/Herring Biosample Database"
REG <- c(HG = "Haida_Gwaii.csv", PRD = "Prince_Rupert_District.csv", CC = "Central_Coast.csv",
         SoG = "Strait_of_Georgia.csv", WCVI = "West_Coast_Vancouver_Island.csv",
         A27 = "Area_27.csv", A2W = "Area_2W.csv")
REG_NAME <- c(HG = "Haida Gwaii", PRD = "Prince Rupert District", CC = "Central Coast",
              SoG = "Strait of Georgia", WCVI = "West Coast Vancouver Island",
              A27 = "Area 27", A2W = "Area 2W")

## ---- 1  occupancy, every region -------------------------------------------
sp <- read_csv(F_SPAWN, show_col_types = FALSE, guess_max = 1e5) %>% keep_surveyed()
## records without a Region label are assigned by the region their statistical area mostly belongs to
lut <- sp %>% filter(!is.na(Region)) %>% count(StatisticalArea, Region) %>%
  group_by(StatisticalArea) %>% slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(StatisticalArea, Region_fill = Region)
sp <- sp %>% left_join(lut, by = "StatisticalArea") %>%
  mutate(reg = coalesce(Region, Region_fill)) %>%
  filter(!is.na(reg), !is.na(Section), Year %in% 1951:2024) %>%
  ## the Strait is defined exactly as in the main analysis: the sections DFO does
  ## not assign to the SoG stock are dropped, so this row reproduces the paper
  filter(!(reg == "SoG" & Section %in% c(131, 133, 134, 136, 293)))

occ_region <- function(rg) {
  d <- sp %>% filter(reg == rg)
  hist_sec <- d %>% distinct(Section, Year) %>% count(Section) %>%
    filter(n >= MIN_ACTIVE) %>% pull(Section)
  if (length(hist_sec) < 3) return(NULL)
  d %>% filter(Section %in% hist_sec) %>% distinct(Section, Year) %>%
    count(Year, name = "n_active") %>%
    complete(Year = 1951:2024, fill = list(n_active = 0)) %>%
    mutate(occupancy = n_active / length(hist_sec), reg = rg, N_all = length(hist_sec))
}
occ <- bind_rows(lapply(names(REG), occ_region))

## ---- 2  age structure, every region ---------------------------------------
## selectivity index: the dominant size-selective gear of each era, reduction
## seine before the closure and roe gillnet after it, as in the main analysis
age_region <- function(rg) {
  f <- file.path(BIO_DB, REG[rg])
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, select = c(2, 3, 15, 17, 22, 26),
              col.names = c("sample", "year", "source", "gear", "weight_g", "age"),
              showProgress = FALSE)
  d <- as_tibble(dt) %>%
    mutate(age = suppressWarnings(as.numeric(age)),
           weight_g = suppressWarnings(as.numeric(weight_g)),
           Year = suppressWarnings(as.numeric(year)),
           Year = ifelse(Year < 30, Year + 2000, ifelse(Year < 100, Year + 1900, Year))) %>%
    filter(!is.na(age), age >= 1, !is.na(Year))

  sel <- d %>%
    mutate(strat = case_when(
      source == "Reduction Fishery" & grepl("seine", gear, ignore.case = TRUE) ~ "reduction",
      source == "Roe Fishery"       & grepl("gill",  gear, ignore.case = TRUE) ~ "gillnet",
      TRUE ~ NA_character_)) %>%
    filter(!is.na(strat))
  ## sample-level proportions, then the unweighted annual mean
  samp <- sel %>% group_by(Year, strat, sample) %>%
    summarise(n = n(), old_n = sum(age >= OLD_AGE),
              old_w = sum(weight_g[age >= OLD_AGE], na.rm = TRUE),
              tot_w = sum(weight_g, na.rm = TRUE), .groups = "drop") %>%
    filter(n >= MIN_FISH)
  ann <- samp %>% group_by(Year, strat) %>%
    summarise(prop_old_n = mean(old_n / n),
              prop_old_wt = ifelse(sum(tot_w) > 0, mean(old_w / pmax(tot_w, 1e-9)), NA_real_),
              n_samples = n(), .groups = "drop")
  ## splice: gillnet where it exists, reduction seine before it
  spliced <- ann %>% group_by(Year) %>%
    summarise(prop_old = ifelse(any(strat == "gillnet"),
                                prop_old_n[strat == "gillnet"][1],
                                prop_old_n[strat == "reduction"][1]),
              source = ifelse(any(strat == "gillnet"), "gillnet", "reduction"),
              n_samples = sum(n_samples), .groups = "drop") %>%
    filter(!is.na(prop_old))
  ## test fishery, where it exists: the population's own age structure
  tf <- d %>% filter(source == "Test Fishery") %>%
    group_by(Year, sample) %>% summarise(n = n(), old_n = sum(age >= OLD_AGE), .groups = "drop") %>%
    filter(n >= MIN_FISH) %>% group_by(Year) %>%
    summarise(pop_old = mean(old_n / n), .groups = "drop")
  spliced %>% left_join(tf, by = "Year") %>% mutate(reg = rg)
}
ages <- bind_rows(lapply(names(REG), age_region))

dat <- occ %>% left_join(ages, by = c("Year", "reg")) %>%
  mutate(region = REG_NAME[reg],
         logit_occ = qlogis(pmin(pmax(occupancy, 0.001), 0.999)))
write_csv(dat, file.path(DERIVED, "cross_region_series.csv"))

## ---- 3  per-region statistics ---------------------------------------------
one_region <- function(rg) {
  d <- dat %>% filter(reg == rg) %>% arrange(Year)
  o <- d %>% filter(!is.na(occupancy))
  ## structural break in occupancy
  br <- tryCatch({
    sf <- sctest(o$occupancy ~ 1, type = "supF")
    bp <- breakpoints(o$occupancy ~ 1, breaks = 1)
    yr <- o$Year[bp$breakpoints]
    yr <- if (length(yr) == 0 || all(is.na(yr))) NA_integer_ else as.integer(yr[1])
    list(supF = unname(sf$statistic)[1], p = sf$p.value[1], yr = yr)
  }, error = function(e) list(supF = NA_real_, p = NA_real_, yr = NA_integer_))
  ## occupancy before and after that break
  pre  <- if (is.na(br$yr)) NA_real_ else mean(o$occupancy[o$Year <= br$yr], na.rm = TRUE)
  post <- if (is.na(br$yr)) NA_real_ else mean(o$occupancy[o$Year >  br$yr], na.rm = TRUE)
  ## occupancy on the selectivity index
  m <- d %>% filter(!is.na(logit_occ), !is.na(prop_old))
  fit <- tryCatch(gamm(logit_occ ~ s(prop_old, k = 5), data = m,
                       correlation = corAR1(form = ~ Year), method = "REML"),
                  error = function(e) NULL)
  st <- if (!is.null(fit)) summary(fit$gam)$s.table else NULL
  tibble(region = unname(REG_NAME[rg]), N_sections = first(d$N_all), n_years = nrow(m),
         supF = round(br$supF), break_year = br$yr,
         occ_pre = round(100*pre), occ_post = round(100*post),
         change_pct = round(100*(post - pre)/pre),
         edf = if (is.null(st)) NA else round(st[1, "edf"], 2),
         P = if (is.null(st)) NA else signif(st[1, "p-value"], 3),
         adj_R2 = if (is.null(fit)) NA else round(summary(fit$gam)$r.sq, 3))
}
res <- bind_rows(lapply(names(REG), one_region))
print(as.data.frame(res), row.names = FALSE)
save_tab(res, "cross_region_results")
cat("\nwritten: Outputs/tables/cross_region_results.csv and derived/cross_region_series.csv\n")

## ---- 4  pooled across regions ---------------------------------------------
z <- function(x) as.numeric(scale(x))
pooled <- dat %>% filter(!is.na(logit_occ), !is.na(prop_old)) %>%
  group_by(reg) %>% filter(n() >= 15) %>% mutate(prop_old_s = z(prop_old)) %>% ungroup() %>%
  mutate(south = ifelse(reg %in% c("SoG", "WCVI"), "southern", "northern"))

## region enters as a random intercept: "south" is nested inside it, so putting
## both in as fixed effects makes the design singular
m_add <- lme(logit_occ ~ prop_old_s, random = ~ 1 | region, data = pooled,
             correlation = corAR1(form = ~ Year | region), method = "ML")
m_int <- lme(logit_occ ~ prop_old_s * south, random = ~ 1 | region, data = pooled,
             correlation = corAR1(form = ~ Year | region), method = "ML")
sa <- summary(m_add)$tTable; si <- summary(m_int)$tTable
cat(sprintf("\nPOOLED  slope on prop_old = %+.2f (P = %.3g), n = %d over %d regions\n",
            sa["prop_old_s","Value"], sa["prop_old_s","p-value"], nrow(pooled), n_distinct(pooled$reg)))
cat(sprintf("SOUTH vs NORTH  interaction P = %.3g | northern slope %+.2f | southern slope %+.2f\n",
            si[grep(":", rownames(si))[1], "p-value"], si["prop_old_s","Value"],
            si["prop_old_s","Value"] + si[grep(":", rownames(si))[1], "Value"]))

## ---- 5  figures ------------------------------------------------------------
ord <- res %>% arrange(desc(adj_R2)) %>% pull(region)
pd <- dat %>% filter(!is.na(occupancy)) %>% mutate(region = factor(region, levels = ord))
brk <- res %>% filter(!is.na(break_year)) %>% mutate(region = factor(region, levels = ord))

p_ts <- ggplot(pd, aes(Year)) +
  geom_vline(data = brk, aes(xintercept = break_year), linetype = "dashed", colour = "grey45") +
  geom_line(aes(y = 100*occupancy), colour = "grey20", linewidth = 0.9) +
  geom_line(aes(y = 100*prop_old), colour = "#B23A34", linewidth = 0.9) +
  facet_wrap(~ region, ncol = 2) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Year", y = "Spawning occupancy (black) and old fish in the selective gear (red)") +
  theme_science + theme(strip.background = element_blank(),
                        strip.text = element_text(size = 15, face = "bold"),
                        axis.title = element_text(size = 15), axis.text = element_text(size = 13))
save_fig(p_ts, "CrossRegion_timeseries", width = 12, height = 12)

p_rel <- ggplot(pd %>% filter(!is.na(prop_old)), aes(100*prop_old, 100*occupancy)) +
  geom_point(size = 1.8, colour = "grey35", alpha = 0.75) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE,
              colour = "#B23A34", fill = "#B23A34", alpha = 0.18, linewidth = 1.2) +
  facet_wrap(~ region, ncol = 2, scales = "free_x") +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Old fish in the most size-selective gear", y = "Spawning occupancy") +
  theme_science + theme(strip.background = element_blank(),
                        strip.text = element_text(size = 15, face = "bold"),
                        axis.title = element_text(size = 15), axis.text = element_text(size = 13))
save_fig(p_rel, "CrossRegion_relationship", width = 12, height = 12)
cat("figures: CrossRegion_timeseries | CrossRegion_relationship\n")

## ---------------------------------------------------------------------------
## The same analysis at two-region resolution: north (CC, HG, PRD) against
## south (SoG, WCVI). Sections and biosamples are pooled before anything is
## computed, so each group is treated as one stock rather than as an average of
## its parts. Areas 27 and 2W are left out: too few usable age years to model,
## and Area 27 sits equidistant between WCVI and the Central Coast.
## ---------------------------------------------------------------------------
GRP <- list(North = c("CC", "HG", "PRD"), South = c("SoG", "WCVI"))

occ_group <- function(gname) {
  d <- sp %>% filter(reg %in% GRP[[gname]])
  hist_sec <- d %>% distinct(reg, Section, Year) %>% count(reg, Section) %>%
    filter(n >= MIN_ACTIVE) %>% select(reg, Section)
  d %>% semi_join(hist_sec, by = c("reg", "Section")) %>%
    distinct(reg, Section, Year) %>% count(Year, name = "n_active") %>%
    complete(Year = 1951:2024, fill = list(n_active = 0)) %>%
    mutate(occupancy = n_active / nrow(hist_sec), grp = gname, N_all = nrow(hist_sec))
}
age_group <- function(gname) {
  ages %>% filter(reg %in% GRP[[gname]]) %>% group_by(Year) %>%
    summarise(prop_old = weighted.mean(prop_old, n_samples, na.rm = TRUE),
              n_samples = sum(n_samples), n_regions = n_distinct(reg), .groups = "drop") %>%
    mutate(grp = gname)
}
gdat <- bind_rows(lapply(names(GRP), occ_group)) %>%
  left_join(bind_rows(lapply(names(GRP), age_group)), by = c("Year", "grp")) %>%
  mutate(logit_occ = qlogis(pmin(pmax(occupancy, 0.001), 0.999)))
write_csv(gdat, file.path(DERIVED, "north_south_series.csv"))

one_group <- function(gname) {
  d <- gdat %>% filter(grp == gname) %>% arrange(Year)
  o <- d %>% filter(!is.na(occupancy))
  sf <- sctest(o$occupancy ~ 1, type = "supF")
  bp <- breakpoints(o$occupancy ~ 1, breaks = 1)
  yr <- o$Year[bp$breakpoints]; yr <- if (length(yr) == 0) NA_integer_ else as.integer(yr[1])
  m <- d %>% filter(!is.na(logit_occ), !is.na(prop_old))
  fit <- gamm(logit_occ ~ s(prop_old, k = 5), data = m,
              correlation = corAR1(form = ~ Year), method = "REML")
  st <- summary(fit$gam)$s.table
  lin <- summary(gls(logit_occ ~ prop_old, data = m, correlation = corAR1(form = ~ Year)))$tTable
  tibble(group = gname, N_sections = first(d$N_all), n_years = nrow(m),
         supF = round(unname(sf$statistic)), break_year = yr,
         occ_pre = round(100*mean(o$occupancy[o$Year <= yr], na.rm = TRUE)),
         occ_post = round(100*mean(o$occupancy[o$Year > yr], na.rm = TRUE)),
         edf = round(st[1, "edf"], 2), P = signif(st[1, "p-value"], 3),
         adj_R2 = round(summary(fit$gam)$r.sq, 3),
         linear_slope = round(lin["prop_old", "Value"], 2),
         linear_P = signif(lin["prop_old", "p-value"], 3))
}
gres <- bind_rows(lapply(names(GRP), one_group))
print(as.data.frame(gres), row.names = FALSE)
save_tab(gres, "north_south_results")

p_ns <- ggplot(gdat %>% filter(!is.na(occupancy)), aes(Year)) +
  geom_line(aes(y = 100*occupancy), colour = "grey20", linewidth = 1.1) +
  geom_line(aes(y = 100*prop_old), colour = "#B23A34", linewidth = 1.1) +
  geom_vline(data = gres %>% filter(!is.na(break_year)) %>% rename(grp = group),
             aes(xintercept = break_year), linetype = "dashed", colour = "grey45") +
  facet_wrap(~ grp, ncol = 1) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "Year", y = "Spawning occupancy (black) and old fish in the selective gear (red)") +
  theme_science + theme(strip.background = element_blank(),
                        strip.text = element_text(size = 17, face = "bold"),
                        axis.title = element_text(size = 16), axis.text = element_text(size = 14))
save_fig(p_ns, "NorthSouth_occupancy", width = 12, height = 9)
cat("\nwritten: north_south_results, north_south_series, NorthSouth_occupancy\n")

## ---------------------------------------------------------------------------
## The population's own age structure: proportion aged 5 and older in the test
## fishery, which samples the spawning aggregation rather than the catch, for the
## five major stocks. Each region gets its own change-point analysis.
## ---------------------------------------------------------------------------
MAIN <- c("SoG", "WCVI", "CC", "PRD", "HG")

tf_dat <- dat %>% filter(reg %in% MAIN, !is.na(pop_old)) %>%
  mutate(region = factor(REG_NAME[reg], levels = unname(REG_NAME[MAIN]))) %>%
  arrange(region, Year)

tf_break <- function(rg) {
  d <- tf_dat %>% filter(reg == rg) %>% arrange(Year)
  ## a change point needs enough years on both sides of it
  if (nrow(d) < 20)
    return(tibble(region = unname(REG_NAME[rg]), n_years = nrow(d),
                  first = min(d$Year), last = max(d$Year),
                  supF = NA_real_, P = NA_real_, break_year = NA_integer_, CI = NA_character_,
                  before = NA_real_, after = NA_real_))
  sf <- tryCatch(sctest(d$pop_old ~ 1, type = "supF"),
                 error = function(e) list(statistic = NA_real_, p.value = NA_real_))
  ## breakpoints() returns the BIC-preferred model, which can be "no break" even
  ## when supF is significant, so the single-break solution is extracted explicitly:
  ## supF tests whether a break exists, this estimates where it is
  bp0 <- tryCatch(breakpoints(d$pop_old ~ 1, h = max(5/nrow(d), 0.15)), error = function(e) NULL)
  bp  <- if (is.null(bp0)) NULL else tryCatch(breakpoints(bp0, breaks = 1), error = function(e) NULL)
  ci <- if (is.null(bp)) NULL else tryCatch(confint(bp0, breaks = 1)$confint, error = function(e) NULL)
  yr <- if (is.null(bp)) NA_integer_ else d$Year[bp$breakpoints]
  yr <- if (length(yr) == 0 || all(is.na(yr))) NA_integer_ else as.integer(yr[1])
  lo <- if (is.null(ci)) NA else d$Year[max(1, ci[1, 1])]
  hi <- if (is.null(ci)) NA else d$Year[min(nrow(d), ci[1, 3])]
  tibble(region = unname(REG_NAME[rg]), n_years = nrow(d),
         first = min(d$Year), last = max(d$Year),
         supF = round(unname(sf$statistic)), P = signif(sf$p.value, 3),
         break_year = yr, CI = if (is.na(yr)) NA else paste0(lo, "-", hi),
         before = if (is.na(yr)) NA else round(100*mean(d$pop_old[d$Year <= yr], na.rm = TRUE)),
         after  = if (is.na(yr)) NA else round(100*mean(d$pop_old[d$Year >  yr], na.rm = TRUE))) %>%
    mutate(change_pct = round(100*(after - before)/before))
}
tf_res <- bind_rows(lapply(MAIN, tf_break))
print(as.data.frame(tf_res), row.names = FALSE)
save_tab(tf_res, "test_fishery_breakpoints_by_region")

seg <- tf_res %>% filter(!is.na(break_year)) %>%
  mutate(region = factor(region, levels = unname(REG_NAME[MAIN])))
lvl <- bind_rows(
  seg %>% transmute(region, x = -Inf, xend = break_year, y = before/100),
  seg %>% transmute(region, x = break_year, xend = Inf, y = after/100))

p_tf <- ggplot(tf_dat, aes(Year, pop_old)) +
  geom_vline(data = seg, aes(xintercept = break_year), linetype = "dashed", colour = "grey40") +
  geom_line(colour = "grey30", linewidth = 0.9) +
  geom_point(size = 1.6, colour = "grey30") +
  geom_segment(data = lvl, aes(x = x, xend = xend, y = y, yend = y),
               colour = "#B23A34", linewidth = 1.3, inherit.aes = FALSE) +
  facet_wrap(~ region, ncol = 1, scales = "free_y") +
  scale_y_continuous(labels = function(x) paste0(100*x, "%")) +
  scale_x_continuous(limits = c(1975, 2024), breaks = seq(1980, 2020, 10)) +
  labs(x = "Year", y = "Old fish (age 5 and older) in the test fishery") +
  theme_science + theme(strip.background = element_blank(),
                        strip.text = element_text(size = 15, face = "bold"),
                        axis.title = element_text(size = 16), axis.text = element_text(size = 14))
save_fig(p_tf, "TestFishery_propold_breakpoints", width = 11, height = 14)
cat("\nwritten: test_fishery_breakpoints_by_region, TestFishery_propold_breakpoints\n")

## ---------------------------------------------------------------------------
## Change points by the paper's method, so the regions are comparable with Fig S3:
## segmented regression on binomial proportions weighted by the number of aged
## fish, breakpoint started at 1986, significance from a Davies test. The earlier
## level-shift test above answers a different question (a change in mean, not in
## trend) and is kept only as a cross-check.
## ---------------------------------------------------------------------------
suppressMessages(library(segmented))

tf_n <- bind_rows(lapply(MAIN, function(rg) {
  f <- file.path(BIO_DB, REG[rg])
  dt <- fread(f, select = c(2, 3, 15, 26), col.names = c("sample", "year", "source", "age"),
              showProgress = FALSE)
  as_tibble(dt) %>%
    mutate(age = suppressWarnings(as.numeric(age)),
           Year = suppressWarnings(as.numeric(year)),
           Year = ifelse(Year < 30, Year + 2000, ifelse(Year < 100, Year + 1900, Year))) %>%
    filter(!is.na(age), age >= 1, !is.na(Year), source == "Test Fishery") %>%
    ## 1975 is excluded, as in the main analysis: its samples are February only,
    ## with no March component, unlike every later year. Including it moves the
    ## Strait's segmented change point from 1985.1 to 1978.7, so the exclusion is
    ## applied to every region for consistency rather than to the Strait alone.
    filter(Year >= 1976) %>%
    group_by(Year, sample) %>% summarise(n = n(), old = sum(age >= OLD_AGE), .groups = "drop") %>%
    filter(n >= MIN_FISH) %>%
    group_by(Year) %>% summarise(prop_old = mean(old / n), n_aged = sum(n), .groups = "drop") %>%
    mutate(reg = rg)
}))

seg_region <- function(rg) {
  d <- tf_n %>% filter(reg == rg) %>% arrange(Year)
  if (nrow(d) < 15) return(tibble(region = unname(REG_NAME[rg]), n_years = nrow(d)))
  m <- glm(prop_old ~ Year, family = binomial, data = d, weights = n_aged)
  sg <- try(segmented(m, seg.Z = ~ Year, psi = list(Year = 1986)), silent = TRUE)
  if (inherits(sg, "try-error")) return(tibble(region = unname(REG_NAME[rg]), n_years = nrow(d)))
  ps <- summary(sg)$psi
  dv <- tryCatch(davies.test(m, ~ Year), error = function(e) NULL)
  tibble(region = unname(REG_NAME[rg]), n_years = nrow(d),
         break_year = round(ps[1, "Est."], 1), SE = round(ps[1, "St.Err"], 2),
         CI = sprintf("%.1f-%.1f", ps[1, "Est."] - 1.96*ps[1, "St.Err"],
                      ps[1, "Est."] + 1.96*ps[1, "St.Err"]),
         davies_P = if (is.null(dv)) NA else signif(dv$p.value, 3))
}
seg_res <- bind_rows(lapply(MAIN, seg_region))
print(as.data.frame(seg_res), row.names = FALSE)
save_tab(seg_res, "test_fishery_segmented_by_region")
cat("\nwritten: test_fishery_segmented_by_region\n")

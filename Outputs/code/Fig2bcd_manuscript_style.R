## =============================================================================
## Figure 2 b, c, d
## b  proportion of older fish by sampling source; quasibinomial GAM with a separate trend
##    for each source (7 fixed degrees of freedom), the model reported in the Results
## c  mean length-at-age: reduction fishery before the closure, test fishery after it
## d  proportion of older fish in the test fishery, with its change point; the 1984-1986
##    collapse years are shaded
## Panel a is rebuilt only so that all four panels share the same widths.
## Panels are saved to Figures/Main/Figure 2.
## Run from Outputs/code:  Rscript Fig2bcd_manuscript_style.R
## =============================================================================

suppressMessages({
  library(tidyverse); library(mgcv); library(segmented); library(patchwork)
  library(scales); library(showtext); library(sysfonts)
})

ALIGNED <- normalizePath("../..")   # repository root; run from Outputs/code
RAW     <- file.path(ALIGNED, "Raw_data")
OUT     <- file.path(ALIGNED, "Figures/Main/Figure 2")
setup <- new.env(); sys.source("00_setup.R", envir = setup)   # biosample column names and year recode

BIO_AREAS <- c(13, 14, 15, 17, 18)
MIN_FISH  <- 30

## ---- style ----
font_add(family = "helvetica", regular = "/System/Library/Fonts/Helvetica.ttc")
showtext_auto()

break_year  <- 1984
shade_start <- break_year
shade_end   <- 2025
MOR_ALPHA   <- 0.6

scale_x_decades <- scale_x_continuous(
  limits = c(1950, 2025),
  breaks = seq(1950, 2020, by = 10),
  expand = expansion(mult = c(0.01, 0.02)))

shade_df  <- tibble(xmin = shade_start, xmax = shade_end, ymin = -Inf, ymax = Inf)
shade_age <- shade_df

theme_science <- theme_classic(base_size = 20, base_family = "helvetica") +
  theme(
    axis.title.x      = element_text(size = 20, margin = margin(t = 14)),
    axis.title.y      = element_text(size = 20, margin = margin(r = 14)),
    axis.text         = element_text(size = 20, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.6),
    axis.line         = element_line(linewidth = 0.5),
    axis.ticks.length = unit(5, "pt"),
    plot.title        = element_text(size = 24, face = "plain"),
    legend.text       = element_text(size = 16),
    legend.title      = element_text(size = 16),
    plot.margin       = margin(15, 15, 15, 15, "pt"))

common_panel <- list(scale_x_decades, theme_science, theme(legend.position = "none"))

fishery_gear_cols <- c(
  "Reduction/Winter/Seine"   = "#2F6F8F",
  "Food & bait/Winter/Seine" = "#A7CDE3",
  "Roe/Spring/Seine"         = "#A52A4A",
  "Roe/Spring/Gillnet"       = "#7B1B38")
SOURCE_COLS <- c(
  "Reduction Fishery" = "#2F6F8F",
  "Food & Bait"       = "#9BC5E8",
  "Roe Seine"         = "#A52A4A",
  "Roe Gillnet"       = "#7B1B38",
  "Test Fishery"      = "#4A4A4A")
age_palette <- setNames(
  c("#2D6A4F","#40916C","#74C69D","#B7E4C7","#FFD166","#F4A261","#E76F51","#C0392B"),
  paste0("Age ", 2:9))

## test fishery = all months, seine sets, from 1977 (as 00_setup.R keep_tf)
tf_ok <- function(source, year, gear)
  source != "Test Fishery" | (as.numeric(year) >= 1977 & gear %in% c("Seine", "Other seine"))
## the biosample extract has no usable header row: names are assigned from 00_setup.R
Biosample_sog <- read_csv(file.path(RAW, "Biosample_Strait_of_Georgia.csv"), show_col_types = FALSE,
                          name_repair = "minimal", guess_max = 100000) %>%
  setNames(setup$BIO_COLS) %>%
  mutate(Year = setup$recode_two_digit_year(year))

## ============================================================
## a - catch by fishery / season / gear (unchanged, alignment only)
## ============================================================
SOG_herring_catch_all_fishing_types <- read.csv(file.path(RAW, "SOG_herring_catch_all_fishing_types.csv")) %>%
  mutate(Season = as.numeric(Season), New_Year = floor(Season / 10) + 1)
SOG_tac <- read.csv(file.path(RAW, "SOG_TAC.csv")) %>% dplyr::select(TAC_metric_tons, New_Year)

df_fs <- SOG_herring_catch_all_fishing_types %>%
  mutate(catch_mt = suppressWarnings(readr::parse_number(as.character(Sum.of.Catch.metric.tons))),
         cat = case_when(
           Fishery == "Reduction fishery"             ~ "Reduction/Winter/Seine",
           Fishery == "Food and bait and special use" ~ "Food & bait/Winter/Seine",
           Fishery == "Roe fishery" & grepl("gill", Gear, ignore.case = TRUE) ~ "Roe/Spring/Gillnet",
           Fishery == "Roe fishery"                   ~ "Roe/Spring/Seine",
           TRUE ~ NA_character_)) %>%
  filter(!is.na(cat)) %>%
  group_by(New_Year, cat) %>%
  summarise(sum_catch = sum(catch_mt, na.rm = TRUE), .groups = "drop") %>%
  filter(New_Year >= 1951) %>%
  tidyr::complete(New_Year, cat, fill = list(sum_catch = 0))
df_fs$cat <- factor(df_fs$cat, levels = names(fishery_gear_cols))

p_fishery <- ggplot(df_fs, aes(x = New_Year, y = sum_catch, fill = cat)) +
  geom_rect(data = shade_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.45) +
  annotate("rect", xmin = 1968, xmax = 1971, ymin = 0, ymax = Inf, fill = "#f4a3a3", alpha = MOR_ALPHA) +
  geom_area(alpha = 0.95, linewidth = 0, position = "stack") +
  geom_line(data = filter(SOG_tac, New_Year >= 1951, New_Year <= 2025),
            aes(x = New_Year, y = TAC_metric_tons, colour = "Total allowable catch (TAC)"),
            linewidth = 1.2, lineend = "round", inherit.aes = FALSE) +
  scale_fill_manual(values = fishery_gear_cols, name = "Fishery / Season / Gear") +
  scale_color_manual(values = c("Total allowable catch (TAC)" = "black"), name = NULL) +
  guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2)) +
  scale_y_continuous(labels = function(x) scales::label_number(big.mark = ",")(x / 1000),
                     expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = "", y = "Catch (x1000 tonnes)", title = NULL) +
  common_panel + theme(legend.position = "inside", legend.position.inside = c(0.97, 0.97),
                       legend.justification = c(1, 1), legend.text = element_text(size = 12),
                       legend.title = element_text(size = 13), legend.background = element_blank(),
                       legend.key = element_blank())

## ============================================================
## b - proportion old fish (age >= 5) by source
## ============================================================
age_prop <- Biosample_sog %>%
  mutate(Year = as.integer(Year),
         age  = suppressWarnings(readr::parse_number(as.character(age))),
         gear = as.character(gear),
         stat_area = suppressWarnings(readr::parse_number(as.character(stat_area))),
         source_clean = case_when(
           source == "Reduction Fishery" ~ "Reduction Fishery",
           grepl("Food|Bait", source) ~ "Food & Bait",
           source == "Roe Fishery" & grepl("seine", gear, ignore.case = TRUE) ~ "Roe Seine",
           source == "Roe Fishery" & grepl("gill",  gear, ignore.case = TRUE) ~ "Roe Gillnet",
           source == "Test Fishery" ~ "Test Fishery",
           TRUE ~ NA_character_),
         old = age >= 5,
         source_clean = factor(source_clean, levels = names(SOURCE_COLS))) %>%
  filter(stock_assessment_region == "Strait of Georgia", stat_area %in% BIO_AREAS,
         tf_ok(source, Year, gear),
         !is.na(Year), !is.na(age), age >= 1, !is.na(source_clean)) %>%
  group_by(Year, source_clean) %>%
  summarise(n_old = sum(old), n = n(), n_young = n - n_old,
            prop_old = n_old / n, .groups = "drop") %>%
  filter(n >= MIN_FISH)

print(age_prop %>% group_by(source_clean) %>%
        summarise(years = paste(range(Year), collapse = "-"), n_years = n(),
                  mean = round(mean(prop_old), 2), sd = round(sd(prop_old), 2)))

m_source <- gam(cbind(n_old, n_young) ~ source_clean + s(Year, by = source_clean, k = 8, fx = TRUE),
                family = quasibinomial, data = age_prop)
cat(sprintf("b: dispersion = %.1f\n", summary(m_source)$scale))

pred_df <- age_prop %>%
  arrange(source_clean, Year) %>%
  group_by(source_clean) %>%
  mutate(run = cumsum(c(1, diff(Year) > 2))) %>%
  group_by(source_clean, run) %>%
  summarise(min_year = min(Year), max_year = max(Year), n_years = n(), .groups = "drop") %>%
  filter(n_years >= 3) %>%
  rowwise() %>% mutate(Year = list(seq(min_year, max_year))) %>%
  tidyr::unnest(Year) %>% dplyr::select(Year, source_clean, run)
pr <- predict(m_source, newdata = pred_df, type = "link", se.fit = TRUE)
pred_df <- pred_df %>%
  mutate(fit = plogis(pr$fit),
         lwr = plogis(pr$fit - 1.96 * pr$se.fit),
         upr = plogis(pr$fit + 1.96 * pr$se.fit))

p_propold <- ggplot() +
  geom_rect(data = shade_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.45) +
  annotate("rect", xmin = 1968, xmax = 1971, ymin = 0, ymax = Inf, fill = "#f4a3a3", alpha = MOR_ALPHA) +
  geom_ribbon(data = pred_df, aes(Year, ymin = lwr, ymax = upr, fill = source_clean,
                                  group = interaction(source_clean, run)),
              alpha = 0.12, colour = NA) +
  geom_line(data = pred_df, aes(Year, fit, colour = source_clean,
                                group = interaction(source_clean, run)), linewidth = 1.2) +
  geom_point(data = age_prop, aes(Year, prop_old, colour = source_clean, size = n), alpha = 0.55) +
  scale_colour_manual(values = SOURCE_COLS, name = "Sample source") +
  scale_fill_manual(values = SOURCE_COLS, guide = "none") +
  scale_size_continuous(range = c(1.2, 6), guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     labels = scales::label_number(accuracy = 0.01),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "Proportion of old fish (age ≥ 5)", title = NULL) +
  common_panel +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.title = element_text(size = 13),
        legend.text = element_text(size = 12), legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill = scales::alpha("white", 0.75), colour = NA))

## ============================================================
## c - mean length at age
## ============================================================
AGES <- 2:9; LINE_AGES <- 2:9; MIN_N <- 30
bio <- Biosample_sog %>%
  mutate(Year = as.integer(Year),
         age = suppressWarnings(readr::parse_number(as.character(age))),
         length_mm = suppressWarnings(readr::parse_number(as.character(length_mm))),
         stat_area = suppressWarnings(readr::parse_number(as.character(stat_area)))) %>%
  filter(stock_assessment_region == "Strait of Georgia", stat_area %in% BIO_AREAS,
         source %in% c("Reduction Fishery","Test Fishery"), tf_ok(source, Year, gear),
         !is.na(Year), !is.na(age), !is.na(length_mm), age %in% AGES)

len_yr <- bio %>%
  group_by(Year, age, source) %>%
  summarise(mean_len_cm = mean(length_mm)/10, se_len_cm = sd(length_mm)/10/sqrt(n()),
            n = n(), .groups = "drop") %>%
  filter(n >= MIN_N) %>%
  mutate(age_lbl = factor(paste0("Age ", age), levels = paste0("Age ", AGES)),
         data_source = ifelse(source == "Reduction Fishery", "Reduction fishery", "Test fishery")) %>%
  arrange(age, data_source, Year) %>%
  group_by(age, data_source) %>%
  mutate(run = cumsum(c(1, diff(Year) > 2))) %>%   # break lines across sampling gaps
  ungroup()

len_lines <- len_yr %>% filter(age %in% LINE_AGES)

## 1980 -> 1983 change and later trend, for the caption
chg <- len_yr %>% filter(data_source == "Test fishery", Year %in% c(1980, 1983)) %>%
  group_by(age) %>% filter(n_distinct(Year) == 2) %>%
  summarise(d_1980_83 = round(mean_len_cm[Year == 1983] - mean_len_cm[Year == 1980], 2), .groups = "drop")
trd <- len_yr %>% filter(data_source == "Test fishery", Year >= 1984, age %in% LINE_AGES) %>%
  group_by(age) %>% summarise(cm_per_decade = round(10 * coef(lm(mean_len_cm ~ Year))[2], 2), .groups = "drop")
print(full_join(chg, trd, by = "age"))

len_labels <- len_yr %>%
  filter(Year > 1973) %>%
  group_by(age_lbl) %>% slice_min(Year, n = 1, with_ties = FALSE) %>% ungroup()
len_labels_end <- len_lines %>% filter(data_source == "Test fishery") %>%
  group_by(age_lbl) %>% slice_max(Year, n = 1, with_ties = FALSE) %>% ungroup() %>%
  filter(Year >= 2020)   # end labels only for age classes sampled to the present
## keep the end labels at least 0.55 cm apart so they do not overlap
len_labels_end <- len_labels_end %>% arrange(mean_len_cm) %>%
  mutate(y_lab = accumulate(mean_len_cm, ~ max(.y, .x + 0.55)))
hdr <- len_labels %>% filter(Year == min(Year)) %>% slice_max(mean_len_cm, n = 1, with_ties = FALSE)

p_length <- ggplot(len_yr, aes(Year, mean_len_cm, colour = age_lbl, linetype = data_source)) +
  geom_rect(data = shade_age, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.45) +
  annotate("rect", xmin = 1968, xmax = 1971, ymin = -Inf, ymax = Inf, fill = "#f4a3a3", alpha = MOR_ALPHA) +
  annotate("rect", xmin = 1980, xmax = 1984, ymin = -Inf, ymax = Inf, fill = "#F4A261", alpha = 0.15) +  # 1980-83 drop
  geom_ribbon(data = len_lines,
              aes(ymin = mean_len_cm - se_len_cm, ymax = mean_len_cm + se_len_cm,
                  fill = age_lbl, group = interaction(age_lbl, data_source, run)),
              alpha = 0.12, colour = NA) +
  geom_line(data = len_lines, aes(group = interaction(age_lbl, data_source, run)),
            linewidth = 0.9, lineend = "round", alpha = 0.9) +
  geom_point(aes(size = n), alpha = 0.65) +
  geom_text(data = len_labels, aes(Year, mean_len_cm, label = age, colour = age_lbl),
            hjust = 1, nudge_x = -0.6, size = 6, fontface = "bold", show.legend = FALSE) +
  geom_text(data = len_labels_end, aes(Year, y_lab, label = age, colour = age_lbl),
            hjust = 0, nudge_x = 0.6, size = 5, fontface = "bold", show.legend = FALSE) +
  annotate("text", x = hdr$Year - 1.2, y = hdr$mean_len_cm + 0.8, label = "Age",   # centred over the number column
           hjust = 0.5, size = 5.5, fontface = "bold", colour = "grey25") +
  scale_colour_manual(values = age_palette, guide = "none") +
  scale_fill_manual(values = age_palette, guide = "none") +
  scale_linetype_manual(values = c("Reduction fishery" = "dashed", "Test fishery" = "solid"),
                        name = "Sample source") +
  scale_size_continuous(range = c(1.2, 7), guide = "none") +
  scale_y_continuous(breaks = seq(12, 24, 2), limits = c(12, 25),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = NULL, y = "Mean length (cm)", title = NULL) +
  common_panel +
  theme(legend.position = "inside", legend.position.inside = c(0.70, 0.98),
        legend.justification = c(0, 1), legend.text = element_text(size = 12),
        legend.title = element_text(size = 13), legend.key.width = unit(1.4, "cm"),
        legend.background = element_rect(fill = scales::alpha("white", 0.8), colour = NA),
        legend.key = element_blank())

## ============================================================
## d - proportion age >= 5, test fishery: change point + phase means
## ============================================================
sog_age <- Biosample_sog %>%
  mutate(age = suppressWarnings(as.numeric(age))) %>%
  filter(Year >= 1977, gear %in% c("Seine","Other seine"),   # all months, from 1977
         source == "Test Fishery", Year >= 1970, !is.na(age)) %>%
  group_by(Year) %>%
  summarise(n_fish = n(), n_age5 = sum(age >= 5), .groups = "drop") %>%
  filter(n_fish >= 30) %>%
  mutate(prop_age5 = n_age5 / n_fish,
         phase = factor(case_when(Year <= 1983 ~ "Pre-collapse",
                                  Year <= 1986 ~ "Collapse",
                                  TRUE         ~ "Post-collapse"),
                        levels = c("Pre-collapse", "Collapse", "Post-collapse")))

m_age5   <- glm(cbind(n_age5, n_fish - n_age5) ~ Year, family = binomial, data = sog_age)
seg_age5 <- segmented(m_age5, seg.Z = ~Year, psi = list(Year = 1986))
sog_age$seg_fit <- fitted(seg_age5, type = "response")
ci_brk <- confint(seg_age5)
cat(sprintf("d: break = %.1f (95%% CI %.1f-%.1f); Davies test in 04 / Table of change points\n",
            seg_age5$psi[1, "Est."], ci_brk[1, 2], ci_brk[1, 3]))

phase_means <- sog_age %>% group_by(phase) %>%
  summarise(m = mean(prop_age5), se = sd(prop_age5)/sqrt(n()), sd = sd(prop_age5),
            x0 = min(Year), x1 = max(Year), n = n(), .groups = "drop")
print(phase_means)
## value labels: pre below its band, collapse to the left of its segment, post above its band
phase_lab <- phase_means %>%
  mutate(lab_x = case_when(phase == "Pre-collapse" ~ (x0 + x1)/2, phase == "Collapse" ~ x0 - 0.4, TRUE ~ 1992),
         lab_y = case_when(phase == "Pre-collapse" ~ m - se - 0.02, phase == "Collapse" ~ m, TRUE ~ m + se + 0.02),
         hj    = if_else(phase == "Collapse", 1, 0.5))
a <- summary(aov(prop_age5 ~ phase, data = sog_age))[[1]]
cat(sprintf("d: ANOVA F(%d,%d) = %.2f, P = %.4f\n", a$Df[1], a$Df[2], a$`F value`[1], a$`Pr(>F)`[1]))

age_bg_full <- list(
  geom_rect(data = shade_age, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.45),
  annotate("rect", xmin = 1968, xmax = 1971, ymin = -Inf, ymax = Inf, fill = "#f4a3a3", alpha = MOR_ALPHA),
  annotate("rect", xmin = 1984, xmax = 1987, ymin = -Inf, ymax = Inf, fill = "grey60", alpha = 0.25))

p_age5 <- ggplot(sog_age, aes(Year, prop_age5)) +
  age_bg_full +
  geom_line(colour = "grey35", linewidth = 0.9, lineend = "round") +
  geom_point(aes(size = n_fish), colour = "grey35", alpha = 0.75, shape = 16) +
  geom_line(aes(y = seg_fit), colour = "#8B1A1A", linewidth = 1.1, lineend = "round") +
  annotate("text", x = 1951, y = max(sog_age$prop_age5), label = "Test fishery\nbegins 1975",
           hjust = 0, vjust = 1, size = 5, colour = "grey40") +
  common_panel + scale_size_area(max_size = 4, guide = "none") +
  labs(x = "Year", y = "Proportion of old fish (age ≥ 5)")

## ---- align + save b, c, d ----
plots_aligned <- align_patches(p_fishery, p_propold, p_length, p_age5)
ggsave(file.path(OUT, "Fig2b_propold_by_source.pdf"),   plot = plots_aligned[[2]], width = 9, height = 6)
ggsave(file.path(OUT, "Fig2c_length_at_age.pdf"),       plot = plots_aligned[[3]], width = 9, height = 6)
ggsave(file.path(OUT, "Fig2d_propold_testfishery.pdf"), plot = plots_aligned[[4]], width = 9, height = 6)
cat("Saved b, c, d to", OUT, "\n")

## =============================================================================
## Figure 1 b, c, d
## b  spawning occupancy (19 historically used sections)
## c  spawning occupancy against SSB before and after the collapse (hysteresis)
## d  latitudinal extent of spawning
## Reads the series written by 01_data_preparation.Rmd.
## Run from Outputs/code:  Rscript Fig1_bcd_SoG_pipeline.R
## =============================================================================

suppressMessages({
  library(tidyverse); library(mgcv); library(scales)
  library(patchwork); library(showtext); library(sysfonts)
})

ALIGNED <- normalizePath("../..")   # repository root; run from Outputs/code
DERIVED <- file.path(ALIGNED, "Outputs/derived")
OUT     <- file.path(ALIGNED, "Figures/Main/Figure 1/SoG_sum_2026_09_15")
dir.create(OUT, showWarnings = FALSE)

YEARS <- 1951:2024

## ---- style ----
cluster_cols_9 <- c(
  "13" = "#5B8C5A", "15" = "#7EB5D6", "14" = "#E07B3F", "17" = "#3A9E8F",
  "18" = "#3461A8", "19" = "#8B6BAE", "16" = "#C03B3B", "28" = "#2D6A4F",
  "29" = "#1A2F6B")

font_add(family = "helvetica", regular = "/System/Library/Fonts/Helvetica.ttc")
showtext_auto()

break_year  <- 1984
shade_start <- break_year
shade_end   <- 2025

scale_x_decades <- scale_x_continuous(
  limits = c(1950, 2025),
  breaks = seq(1950, 2020, by = 10),
  expand = expansion(mult = c(0.01, 0.02)))

shade_df <- tibble(xmin = shade_start, xmax = shade_end, ymin = -Inf, ymax = Inf)

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

## ---- Data from the pipeline ----
annual <- read_csv(file.path(DERIVED, "annual_series.csv"), show_col_types = FALSE) %>%
  filter(Year %in% YEARS, !is.na(occupancy))

site_year <- read_csv(file.path(DERIVED, "site_year_spawn.csv"), show_col_types = FALSE) %>%
  filter(Year %in% YEARS)   # every recorded spawn; records without a biomass estimate have spawn_bio_sum = 0

N_all <- round(unique(annual$n_active / annual$occupancy)[1])
cat(sprintf("Occupancy denominator: %d sections | years %d-%d\n",
            N_all, min(annual$Year), max(annual$Year)))

## ============================================================
## Panel b - spawning occupancy
## ============================================================
m1_occupancy <- annual %>% select(Year, n_active, occupancy)

p_occupancy <- ggplot(m1_occupancy, aes(Year, occupancy)) +
  geom_rect(data = shade_df, aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf),
            fill = "grey85", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(linewidth = 0.9, lineend = "round") +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 10),
    se = TRUE,
    color = "darkblue",
    fill = "grey75",
    linewidth = 1
  )+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)),
                     labels = scales::percent_format(accuracy=1)) +
  labs(x="Year", y="Spawning occupancy", title=NULL) +
  common_panel

## ============================================================
## Panel c - biomass-occupancy hysteresis (SSB = sum of all records)
## ============================================================
post_lab <- sprintf("Post-collapse (%d-%d)", break_year, max(annual$Year))
pre_lab  <- sprintf("Pre-collapse (%d-%d)", min(annual$Year), break_year - 1)

d <- annual %>%
  transmute(Year, occupancy, SSB, SSBk = SSB/1000,
            phase = factor(if_else(Year < break_year, pre_lab, post_lab),
                           levels = c(pre_lab, post_lab))) %>%
  arrange(Year)

## statistics check (should match 03_biomass_occupancy_decoupling)
mh <- lm(occupancy ~ scale(SSB) + phase, data = d)
print(round(summary(mh)$coefficients, 4))
pre      <- d %>% filter(phase == pre_lab)
post_all <- d %>% filter(phase == post_lab)
ov       <- range(pre$SSB)
post_in  <- post_all %>% filter(SSB >= ov[1], SSB <= ov[2])
cat(sprintf("Matched SSB: pre %.0f%% (n=%d) vs post %.0f%% (n=%d)\n",
            100*mean(pre$occupancy), nrow(pre), 100*mean(post_in$occupancy), nrow(post_in)))
cat(sprintf("Occupancy: pre %.1f +/- %.1f sections (%.0f%%), post %.1f +/- %.1f (%.0f%%)\n",
            mean(pre$occupancy*N_all), sd(pre$occupancy*N_all), 100*mean(pre$occupancy),
            mean(post_all$occupancy*N_all), sd(post_all$occupancy*N_all), 100*mean(post_all$occupancy)))
cat(sprintf("SSB: min %.0f kt (%d), max %.0f kt (%d), pre max %.0f kt, means %.0f vs %.0f kt\n",
            min(d$SSBk), d$Year[which.min(d$SSB)], max(d$SSBk), d$Year[which.max(d$SSB)],
            max(pre$SSBk), mean(pre$SSBk), mean(post_all$SSBk)))

lev <- levels(d$phase)
phase_fill <- setNames(c("grey85", "grey35"), lev)
phase_edge <- setNames(c("grey45", "grey15"), lev)

hull_fill  <- setNames(c("grey80", "grey30"), lev)

pre_hull  <- pre[chull(pre$SSBk, pre$occupancy), ]
post_hull <- post_all[chull(post_all$SSBk, post_all$occupancy), ]

pd <- ggplot(d, aes(x = SSBk, y = occupancy)) +
  geom_polygon(data = pre_hull, aes(x = SSBk, y = occupancy),
               inherit.aes = FALSE, fill = hull_fill[[1]], alpha = 0.28) +
  geom_polygon(data = post_hull, aes(x = SSBk, y = occupancy),
               inherit.aes = FALSE, fill = hull_fill[[2]], alpha = 0.22) +
  geom_path(colour = "grey60", linewidth = 0.4, alpha = 0.55) +
  geom_point(aes(fill = phase, colour = phase), shape = 21, size = 3.4, stroke = 0.5) +
  scale_fill_manual(values = phase_fill, name = NULL) +
  scale_colour_manual(values = phase_edge, name = NULL, guide = "none") +
  scale_x_continuous(
    "SSB (×1000 tonnes)",
    labels = scales::label_number(big.mark = ","),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_continuous(
    "Spawning occupancy",
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  labs(title = NULL) +
  theme_science +
  theme(
    legend.text          = element_text(size = 16),
    legend.position      = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background    = element_rect(fill = scales::alpha("white", 0.6), colour = NA),
    legend.key           = element_blank(),
    panel.grid           = element_blank()
  )

## ============================================================
## Panel d - spawning latitudinal extent
## Bubbles = occupied location-years (SoG), size = sum of all spawn records
## ============================================================
spawning_sog <- site_year %>%
  mutate(Spawning_biomass = spawn_bio_sum,
         StatisticalArea  = factor(StatisticalArea))

lat_year_q <- spawning_sog %>%
  group_by(Year) %>%
  summarise(lat05 = quantile(Latitude, 0.05, na.rm = TRUE),
            lat95 = quantile(Latitude, 0.95, na.rm = TRUE), .groups = "drop")

## Southern (5%) edge: automatic smoothing gives a straight line on these data (edf 1, also with
## REML), so the wiggliness is fixed at 11 df, the edf of the Aug version, for display (fx = TRUE).
## The pre/post shift is the same either way (48.78 -> 49.02 degN).
g05 <- gam(lat05 ~ s(Year, k = 12, fx = TRUE), data = lat_year_q)
g95 <- gam(lat95 ~ s(Year, k = 12, fx = TRUE), data = lat_year_q)  # top edge: fixed 11 df (auto edf was 7.9), display only

year_grid <- data.frame(Year = seq(min(lat_year_q$Year), max(lat_year_q$Year), length.out = 200))
pred05 <- predict(g05, newdata = year_grid, se.fit = TRUE)
pred95 <- predict(g95, newdata = year_grid, se.fit = TRUE)

quant_col <- "#696969"
pred_05 <- tibble(Year = year_grid$Year, fit = pred05$fit, se = pred05$se.fit)
pred_95 <- tibble(Year = year_grid$Year, fit = pred95$fit, se = pred95$se.fit)

p2 <- ggplot(spawning_sog, aes(x = Year, y = Latitude)) +
  annotate("rect", xmin = shade_start, xmax = shade_end, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.45) +
  geom_point(data = filter(spawning_sog, Spawning_biomass > 0),
             aes(size = Spawning_biomass, color = StatisticalArea), alpha = 0.65) +
  # recorded spawn without a biomass estimate (e.g. Area 28, Howe Sound): smallest dot size
  geom_point(data = filter(spawning_sog, Spawning_biomass == 0),
             aes(color = StatisticalArea), size = 0.6, alpha = 0.65) +
  geom_ribbon(data = pred_05, inherit.aes = FALSE,
              aes(x = Year, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = quant_col, alpha = 0.10) +
  geom_ribbon(data = pred_95, inherit.aes = FALSE,
              aes(x = Year, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = quant_col, alpha = 0.10) +
  geom_line(data = pred_05, inherit.aes = FALSE, aes(x = Year, y = fit),
            color = quant_col, linewidth = 0.9) +
  geom_line(data = pred_95, inherit.aes = FALSE, aes(x = Year, y = fit),
            color = quant_col, linewidth = 0.9) +
  scale_size_continuous(name = "Spawning biomass (t)", range = c(0.6, 5.5),
                        labels = scales::label_number(big.mark = ",")) +
  scale_color_manual(values = cluster_cols_9, name = "Statistical area") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = "Year", y = "Latitude", title = NULL) +
  common_panel

## ---- Save (same size as the Rmd: 9 x 6 in, aligned panels) ----
plots_aligned <- align_patches(p_occupancy, pd, p2)
ggsave(file.path(OUT, "p_occupancy.pdf"),      plot = plots_aligned[[1]], width = 9, height = 6)
ggsave(file.path(OUT, "hystersis.pdf"),        plot = plots_aligned[[2]], width = 9, height = 6)
ggsave(file.path(OUT, "SOG_quantile_gam.pdf"), plot = plots_aligned[[3]], width = 9, height = 6)
cat("Saved to", OUT, "\n")

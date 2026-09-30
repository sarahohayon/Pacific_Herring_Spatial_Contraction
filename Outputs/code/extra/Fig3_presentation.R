## ============================================================
## Figure 3 for the talk: the two series on their own, at projector size.
##   - the model prediction and its confidence band are dropped; they sat on top
##     of the two trend lines, which are the message
##   - panel b becomes an inset: observed variance explained against its own
##     randomised null, five drivers only, large type
## Run from Outputs/code: Rscript Fig3_presentation.R
## ============================================================
suppressMessages({library(tidyverse); library(scales); library(patchwork)
  library(showtext); library(sysfonts)})
font_add(family = "helvetica", regular = "/System/Library/Fonts/Helvetica.ttc")
showtext_auto(); showtext_opts(dpi = 170)
ALIGNED <- "../.."
OUT <- file.path("..", "figures")

col_occ <- "grey20"; col_old <- "#B23A34"; col_sst <- "#2A9D8F"; col_ssb <- "#1F3B73"

theme_talk <- theme_classic(base_size = 26, base_family = "helvetica") +
  theme(axis.title.x = element_text(size = 26, margin = margin(t = 12)),
        axis.title.y = element_text(size = 26, margin = margin(r = 12)),
        axis.text    = element_text(size = 24, colour = "black"),
        axis.ticks   = element_line(linewidth = 0.7),
        axis.line    = element_line(linewidth = 0.6),
        plot.margin  = margin(12, 12, 12, 12, "pt"))

annual <- read_csv(file.path(ALIGNED, "Outputs/derived/annual_series.csv"), show_col_types = FALSE) %>%
  filter(Year %in% 1951:2024)
sel <- annual %>% filter(!is.na(prop_old_wt)) %>%
  mutate(gear = ifelse(Year <= 1971, "Reduction seine", "Roe gillnet"))

LAB_OCC <- "Spawning occupancy"
LAB_OLD <- "Old fish in the catch (age 5+)"
LEV <- c(LAB_OCC, LAB_OLD)

p_main <- ggplot() +
  annotate("rect", xmin = 1984, xmax = 2024, ymin = 0, ymax = 1, fill = "grey78", alpha = 0.5) +
  annotate("rect", xmin = 1968, xmax = 1971, ymin = 0, ymax = 1, fill = "#f4a3a3", alpha = 0.45) +
  annotate("text", x = 2024, y = 0.075, hjust = 1, size = 8.5, colour = "grey15",
           label = "Adjusted R\u00b2 = 0.54") +
  annotate("text", x = 2024, y = 0.015, hjust = 1, size = 8.5, colour = "grey15",
           label = "P < 0.001") +
  geom_line(data = annual, aes(Year, occupancy, colour = LAB_OCC), linewidth = 1.8, lineend = "round") +
  geom_point(data = annual, aes(Year, occupancy, colour = LAB_OCC), size = 3.0, shape = 16) +
  geom_line(data = sel, aes(Year, prop_old_wt, group = gear, colour = LAB_OLD),
            linewidth = 1.8, lineend = "round") +
  geom_point(data = sel, aes(Year, prop_old_wt, colour = LAB_OLD, shape = gear), size = 3.6) +
  scale_colour_manual(values = setNames(c(col_occ, col_old), LEV), breaks = LEV, name = NULL) +
  scale_shape_manual(values = c("Reduction seine" = 17, "Roe gillnet" = 16), guide = "none") +
  guides(colour = guide_legend(override.aes = list(linewidth = 2.4, size = 4))) +
  scale_y_continuous(NULL, labels = percent_format(accuracy = 1), limits = c(0, 1),
                     breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.02, 0.06))) +
  scale_x_continuous("Year", breaks = seq(1950, 2020, 10), limits = c(1951, 2024),
                     expand = expansion(mult = c(0.01, 0.02))) +
  theme_talk +
  theme(legend.position = "top", legend.justification = "left",
        legend.text = element_text(size = 24), legend.key.width = unit(2.6, "lines"),
        legend.margin = margin(b = 4), legend.key = element_blank())

## ---- inset: observed variance explained against its randomised null ----------
nd <- read_csv(file.path(ALIGNED, "Outputs/derived/fig3b_null_draws.csv"), show_col_types = FALSE) %>%
  mutate(driver = recode(driver,
           "Old fish removed" = "Old fish", "Total catch" = "Catch",
           "Mean SST" = "SST", "Pacific Decadal Oscillation" = "PDO",
           "Oceanic Nino Index" = "ONI", "Spawning stock biomass" = "SSB")) %>%
  filter(driver %in% c("Old fish", "Catch", "SST", "PDO", "ONI", "SSB")) %>%
  mutate(class = case_when(driver %in% c("Old fish", "Catch") ~ "Fishery",
                           driver %in% c("SST", "PDO", "ONI") ~ "Climate",
                           TRUE ~ "SSB"))
ord <- nd %>% distinct(driver, observed) %>% arrange(observed) %>% pull(driver)
nd <- nd %>% mutate(driver = factor(driver, levels = ord))
CLASS_COLS <- c("Fishery" = col_old, "Climate" = col_sst, "SSB" = col_ssb)

dens <- nd %>% group_by(driver, class) %>%
  summarise(d = list(density(null_r2, n = 256)), .groups = "drop") %>%
  mutate(out = map(d, ~tibble(x = .x$x, h = .x$y/max(.x$y)))) %>%
  select(-d) %>% unnest(out) %>%
  mutate(yi = as.numeric(driver), ymin = yi - 0.05, ymax = yi - 0.05 + 0.62*h)
pts <- nd %>% distinct(driver, class, observed) %>% mutate(yi = as.numeric(driver))

p_inset <- ggplot() +
  geom_ribbon(data = dens, aes(x = x, ymin = ymin, ymax = ymax, group = driver),
              fill = "grey84", colour = "grey58", linewidth = 0.3) +
  geom_segment(data = pts, aes(x = observed, xend = observed, y = yi - 0.12, yend = yi + 0.5,
                               colour = class), linewidth = 2.2, lineend = "round") +
  scale_colour_manual(values = CLASS_COLS, name = NULL) +
  scale_y_continuous(NULL, breaks = seq_along(ord), labels = ord,
                     expand = expansion(add = c(0.2, 0.6))) +
  scale_x_continuous("Variance explained", limits = c(-0.45, 0.72), expand = c(0, 0)) +
  theme_classic(base_size = 22, base_family = "helvetica") +
  theme(axis.text.y = element_text(size = 22, colour = "black"),
        axis.text.x = element_text(size = 19, colour = "black"),
        axis.title.x = element_text(size = 22, margin = margin(t = 8)),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.background = element_blank(),
        plot.margin = margin(6, 8, 6, 6, "pt"))

## the two series on their own, for a slide that builds
ggsave(file.path(OUT, "Fig3_talk_trends_only.pdf"), p_main, width = 16, height = 9)
ggsave(file.path(OUT, "Fig3_talk_trends_only.png"), p_main, width = 16, height = 9, dpi = 200)

## side by side: b as a small companion, no data hidden
side <- p_main + p_inset + plot_layout(widths = c(1, 0.42))

p_box <- p_inset + theme(plot.background = element_rect(fill = alpha("white", 0.93), colour = "grey60"),
                         plot.margin = margin(8, 10, 8, 8, "pt"))
over <- p_main + inset_element(p_box, left = 0.015, bottom = 0.03, right = 0.46, top = 0.46)
ggsave(file.path(OUT, "Fig3_talk_inset.png"), over, width = 16, height = 9, dpi = 170)
ggsave(file.path(OUT, "Fig3_talk.pdf"), side, width = 20, height = 9)
ggsave(file.path(OUT, "Fig3_talk.png"), side, width = 20, height = 9, dpi = 170)

cat("written: Outputs/figures/Fig3_talk*, Fig3_talk_trends_only*\n")

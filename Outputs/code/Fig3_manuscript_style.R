## =============================================================================
## Figure 3
## a  spawning occupancy and the proportion of older fish in the dominant gear of each era,
##    with occupancy predicted from older fish (AR(1) GAMM, 95% band). Lines break over the
##    1968-1971 closure, where the index switches from reduction seines to roe gillnets.
## b  variance in occupancy explained by each candidate driver, against its own null from
##    1,000 randomizations (draws written by Fig3b_null_draws.R)
## Panels are saved separately and as one combined figure to Figures/Main/Figure 3.
## Run from Outputs/code:  Rscript Fig3_manuscript_style.R
## =============================================================================

suppressMessages({
  library(tidyverse); library(mgcv); library(nlme); library(patchwork)
  library(scales); library(showtext); library(sysfonts)
})

ALIGNED <- normalizePath("../..")   # repository root; run from Outputs/code
OUT     <- file.path(ALIGNED, "Figures/Main/Figure 3")
dir.create(OUT, showWarnings = FALSE)

font_add(family = "helvetica", regular = "/System/Library/Fonts/Helvetica.ttc")
font_add(family = "tagbold", regular = "/System/Library/Fonts/Supplemental/Arial Bold.ttf")   # bold panel letters
showtext_auto()

col_occ <- "#555555"; col_old <- "#C0392B"; col_sst <- "#2A9D8F"; col_null <- "grey55"
col_fit <- "darkblue"   # model prediction, as the GAM smooths in Figure 1
MOR_ALPHA <- 0.6

scale_x_decades <- scale_x_continuous(
  limits = c(1950, 2025), breaks = seq(1950, 2020, by = 10),
  expand = expansion(mult = c(0.01, 0.02)))
shade_df <- tibble(xmin = 1984, xmax = 2025, ymin = -Inf, ymax = Inf)

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
bands <- list(
  geom_rect(data = shade_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.45),
  annotate("rect", xmin = 1968, xmax = 1971, ymin = -Inf, ymax = Inf,
           fill = "#f4a3a3", alpha = MOR_ALPHA))

annual <- read_csv(file.path(ALIGNED, "Outputs/derived/annual_series.csv"), show_col_types = FALSE) %>%
  filter(Year %in% 1951:2024)
dm <- annual %>%
  filter(!is.na(occupancy), !is.na(prop_old_wt), !is.na(mean_SST), !is.na(ONI_ann), !is.na(mean_PDO)) %>%
  arrange(Year)
cat(sprintf("Model frame: %d years (%d-%d)\n", nrow(dm), min(dm$Year), max(dm$Year)))

## ============================================================
## a - occupancy and the selectivity index, one 0-100% scale, gears kept apart
## ============================================================
era_of <- function(y) if_else(y <= 1965, "Reduction seine", "Roe gillnet")
sel <- annual %>% filter(!is.na(prop_old_wt)) %>%
  mutate(gear = factor(era_of(Year), levels = c("Reduction seine", "Roe gillnet")))

p3a <- ggplot() +
  bands +
  geom_line(data = annual, aes(Year, occupancy), colour = col_occ, linewidth = 0.9, lineend = "round") +
  geom_point(data = annual, aes(Year, occupancy), colour = col_occ, size = 1.8, alpha = 0.8) +
  geom_line(data = sel, aes(Year, prop_old_wt, group = gear), colour = col_old,
            linewidth = 0.9, lineend = "round") +
  geom_point(data = sel, aes(Year, prop_old_wt, shape = gear), colour = col_old, size = 2.4) +
  scale_shape_manual(values = c("Reduction seine" = 17, "Roe gillnet" = 16), name = NULL) +
  scale_y_continuous("Spawning occupancy", labels = percent_format(accuracy = 1),
                     limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     sec.axis = sec_axis(~ ., name = "Proportion of old fish removed",
                                         labels = percent_format(accuracy = 1),
                                         breaks = seq(0, 1, 0.25))) +
  labs(x = NULL) +
  common_panel +
  theme(axis.title.y.right = element_text(colour = col_old, margin = margin(l = 14)),
        axis.text.y.right   = element_text(colour = col_old),
        axis.title.y.left   = element_text(colour = col_occ),
        axis.text.y.left    = element_text(colour = col_occ),
        legend.position = "inside", legend.position.inside = c(0.98, 0.02),
        legend.justification = c(1, 0), legend.text = element_text(size = 14),
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.key = element_blank())

## ============================================================
## models: AR(1) GAMMs, as in 05_occupancy_drivers_GAMM.Rmd
## ============================================================
ctrl <- corAR1(form = ~ Year)
m_old  <- gamm(logit_occ ~ s(prop_old_wt, k = 5), data = dm, correlation = ctrl, method = "REML")
m_sst  <- gamm(logit_occ ~ s(mean_SST, k = 5),    data = dm, correlation = ctrl, method = "REML")
m_null <- gls(logit_occ ~ 1, data = dm, correlation = ctrl, method = "REML")
st <- summary(m_old$gam)$s.table
cat(sprintf("prop_old alone: edf %.2f, F %.2f, P %.3g, adj R2 %.2f | SST alone: adj R2 %.2f\n",
            st[1, "edf"], st[1, "F"], st[1, "p-value"], summary(m_old$gam)$r.sq,
            summary(m_sst$gam)$r.sq))

## within-era check: does occupancy track old fish inside the roe era alone?
roe <- dm %>% filter(Year >= 1972)
red <- dm %>% filter(Year <= 1965)
for (dd in list(list(roe, "roe gillnet era 1972-2024"), list(red, "reduction era 1951-1965"))) {
  d0 <- dd[[1]]
  if (nrow(d0) > 5) {
    g <- gls(logit_occ ~ prop_old_wt, data = d0, correlation = corAR1(form = ~ Year), method = "ML")
    tt <- summary(g)$tTable
    cat(sprintf("within %s (n = %d): slope = %.2f (SE %.2f), P = %.3f\n",
                dd[[2]], nrow(d0), tt[2, 1], tt[2, 2], tt[2, 4]))
  }
}

## ============================================================
## b - partial effect of the proportion of old fish, with partial residuals
## ============================================================
grid_b <- tibble(prop_old_wt = seq(min(dm$prop_old_wt), max(dm$prop_old_wt), length.out = 200))
pr_b   <- predict(m_old$gam, newdata = grid_b, se.fit = TRUE)
b_line <- grid_b %>% mutate(fit = plogis(pr_b$fit),
                            lwr = plogis(pr_b$fit - 1.96 * pr_b$se.fit),
                            upr = plogis(pr_b$fit + 1.96 * pr_b$se.fit))
b_pts  <- dm %>% mutate(part = plogis(fitted(m_old$gam) + residuals(m_old$gam)),
                        gear = factor(era_of(Year), levels = c("Reduction seine", "Roe gillnet")))

p3b <- ggplot() +
  geom_ribbon(data = b_line, aes(prop_old_wt, ymin = lwr, ymax = upr), fill = col_old, alpha = 0.15) +
  geom_line(data = b_line, aes(prop_old_wt, fit), colour = col_old, linewidth = 1.3) +
  geom_point(data = b_pts, aes(prop_old_wt, part, shape = gear), colour = col_occ,
             size = 2.6, alpha = 0.8) +
  scale_shape_manual(values = c("Reduction seine" = 17, "Roe gillnet" = 16), name = NULL) +
  scale_x_continuous("Proportion of old fish removed", labels = percent_format(accuracy = 1)) +
  scale_y_continuous("Spawning occupancy", labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.05, 0.12))) +
  theme_science +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98),
        legend.justification = c(1, 1), legend.text = element_text(size = 14),
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.key = element_blank())

## ============================================================
## c - observed occupancy and what each model predicts
## ============================================================
fit_line <- function(m, label) dm %>%
  transmute(Year, fit = plogis(if (inherits(m, "gls")) fitted(m) else fitted(m$gam)), model = label)
MODELS <- c("Old fish", "Sea-surface temperature", "No driver (mean)")
c_fits <- bind_rows(fit_line(m_old, MODELS[1]), fit_line(m_sst, MODELS[2]), fit_line(m_null, MODELS[3])) %>%
  mutate(model = factor(model, levels = MODELS)) %>%
  arrange(model, Year) %>% group_by(model) %>%
  mutate(run = cumsum(c(1, diff(Year) > 1))) %>% ungroup()   # no fitted line over the 1966-1971 gap
mod_cols <- setNames(c(col_old, col_sst, col_null), MODELS)

p3c <- ggplot() +
  bands +
  geom_line(data = annual, aes(Year, occupancy), colour = col_occ, linewidth = 0.9, lineend = "round") +
  geom_point(data = annual, aes(Year, occupancy), colour = col_occ, size = 1.6, alpha = 0.7) +
  geom_line(data = c_fits, aes(Year, fit, colour = model, group = interaction(model, run)),
            linewidth = 1.2, lineend = "round") +
  scale_colour_manual(values = mod_cols, name = NULL) +
  scale_y_continuous("Spawning occupancy", labels = percent_format(accuracy = 1),
                     limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Year") +
  common_panel +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98),
        legend.justification = c(1, 1), legend.text = element_text(size = 14),
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.key = element_blank())

## ============================================================
## d - model support (dAIC), best model at zero
## ============================================================
tabS2 <- read_csv(file.path(ALIGNED, "Outputs/tables/TableS2_model_selection.csv"), show_col_types = FALSE)
lab <- c(prop_old = "Old fish", SST = "SST", ONI = "ONI", PDO = "PDO", null = "No driver")
pretty <- function(x) {
  x <- gsub("prop_old", "Old fish", x); x <- gsub("\\*", " × ", x); x <- gsub("\\+", " + ", x)
  gsub("null", "No driver", x)
}
d_tab <- tabS2 %>%
  mutate(Model = pretty(Model),
         Class = recode(Class, fishing = "Old fish", climate = "Climate",
                        `fish+clim` = "Old fish + climate", interaction = "Interaction", null = "No driver"),
         Model = fct_reorder(Model, -dAIC))
class_cols <- c("Old fish" = col_old, "Climate" = col_sst, "Old fish + climate" = "#7B5B9D",
                "Interaction" = "grey35", "No driver" = col_null)

p3d <- ggplot(d_tab, aes(dAIC, Model, colour = Class)) +
  geom_segment(aes(x = 0, xend = dAIC, yend = Model), linewidth = 0.6) +
  geom_point(size = 3.4) +
  geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.5) +
  scale_colour_manual(values = class_cols, name = NULL) +
  scale_x_continuous("ΔAIC (0 = best-supported model)", expand = expansion(mult = c(0.02, 0.08))) +
  labs(y = NULL) +
  theme_science +
  theme(axis.text.y = element_text(size = 15),
        legend.position = "inside", legend.position.inside = c(0.98, 0.98),
        legend.justification = c(1, 1), legend.text = element_text(size = 14),
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        legend.key = element_blank())

## ============================================================
## e - ONE panel carrying the whole argument: the two series plus the model that links them
##     grey, observed occupancy; red band, occupancy predicted from the proportion of old fish
##     (AR(1) GAMM, 95% band); red symbols, the selectivity index itself (right axis, gears apart)
## ============================================================
pr_e <- predict(m_old$gam, se.fit = TRUE)
e_fit <- dm %>% transmute(Year,
                          fit = plogis(pr_e$fit),
                          lwr = plogis(pr_e$fit - 1.96 * pr_e$se.fit),
                          upr = plogis(pr_e$fit + 1.96 * pr_e$se.fit)) %>%
  arrange(Year) %>% mutate(run = cumsum(c(1, diff(Year) > 1)))

LAB_OCC <- "Spawning occupancy"
LAB_OLD <- "Old fish (age \u2265 5) in catch"
LAB_FIT <- "Model prediction with 95% CI"
LEV <- c(LAB_OCC, LAB_OLD, LAB_FIT)

p3e <- ggplot() +
  bands +
  geom_ribbon(data = e_fit, aes(Year, ymin = lwr, ymax = upr, group = run),
              fill = "grey75", alpha = 0.55) +
  geom_line(data = annual, aes(Year, occupancy, colour = LAB_OCC), linewidth = 0.9, lineend = "round") +
  geom_point(data = annual, aes(Year, occupancy, colour = LAB_OCC), size = 1.8, alpha = 0.85, shape = 16) +
  geom_line(data = sel, aes(Year, prop_old_wt, group = gear, colour = LAB_OLD),
            linewidth = 0.7, alpha = 0.55, lineend = "round") +
  geom_point(data = sel, aes(Year, prop_old_wt, colour = LAB_OLD, shape = gear), size = 2.2, alpha = 0.75) +
  geom_line(data = e_fit, aes(Year, fit, group = run, colour = LAB_FIT), linewidth = 1.3) +
  scale_colour_manual(values = setNames(c(col_occ, col_old, col_fit), LEV), breaks = LEV, name = NULL) +
  scale_shape_manual(values = c("Reduction seine" = 17, "Roe gillnet" = 16), name = NULL,
                     guide = guide_legend(order = 2,
                       override.aes = list(colour = col_old, size = 2.8, linetype = "blank"))) +
  guides(colour = guide_legend(order = 1, override.aes = list(
           shape     = c(16, NA, NA),
           linetype  = c("solid", "solid", "solid"),
           linewidth = c(0.9, 0.7, 1.3),
           size      = c(1.8, 0, 0)))) +
  scale_y_continuous("Percent", labels = percent_format(accuracy = 1),
                     limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.02, 0.34))) +
  labs(x = "Year") +
  common_panel +
  theme(legend.position = "inside", legend.position.inside = c(0.97, 0.965),
        legend.justification = c(1, 1), legend.text = element_text(size = 18),
        axis.text = element_text(size = 24, colour = "black"),
        axis.title.x = element_text(size = 26, margin = margin(t = 14)),
        axis.title.y = element_text(size = 26, margin = margin(r = 14)),
        legend.key.height = unit(1.1, "lines"), legend.key.width = unit(1.8, "lines"),
        legend.box = "vertical", legend.box.just = "right", legend.box.spacing = unit(0, "pt"),
        legend.justification.inside = c(1, 1), legend.location = "plot",
        legend.margin = margin(1, 3, 1, 3),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key = element_blank(), legend.spacing.y = unit(1, "pt"))

## ============================================================
## b (final) - every candidate driver against its own randomised null
##     grey density = adjusted R2 from 1000 randomisations in which that driver is shifted in
##                    time within each collapse phase (temporal structure kept, alignment with
##                    occupancy destroyed); mark = observed adjusted R2, coloured by driver class
##     Draws: Outputs/derived/fig3b_null_draws.csv (Fig3b_null_draws.R + _extra.R)
## ============================================================
CLASS_COLS <- c("Fishery" = col_old, "Climate" = col_sst, "SSB" = "#1F3B73")
nd <- read_csv(file.path(ALIGNED, "Outputs/derived/fig3b_null_draws.csv"), show_col_types = FALSE) %>%
  mutate(driver = recode(driver,
           "Old fish removed"            = "Proportion of old fish in catch",
           "Total catch"                 = "Total catch\n(positive association)",
           "Oceanic Nino Index"          = "ONI",
           "Pacific Decadal Oscillation" = "PDO",
           "Spawning stock biomass"      = "SSB"),
         class = case_when(
           driver %in% c("Mean SST", "ONI", "PDO") ~ "Climate",
           driver == "SSB"                         ~ "SSB",
           TRUE                                    ~ "Fishery"))
ord <- nd %>% distinct(driver, observed) %>% arrange(observed) %>% pull(driver)
nd  <- nd %>% mutate(driver = factor(driver, levels = ord))
X_MIN <- -0.07; X_MAX <- 0.74

dens <- nd %>% group_by(driver, class) %>%
  group_modify(function(x, k) { d <- density(x$null_r2, adjust = 1.1); tibble(x = d$x, h = d$y/max(d$y)) }) %>%
  ungroup() %>%
  mutate(yi = as.numeric(driver), ymin = yi - 0.04, ymax = yi - 0.04 + 0.60 * h) %>%
  filter(x >= X_MIN, x <= X_MAX)

pts <- nd %>% distinct(driver, class, observed, p_emp) %>%
  mutate(yi   = as.numeric(driver),
         off  = observed < X_MIN,
         xpos = pmax(observed, X_MIN + 0.004),
         lab  = if_else(p_emp < 0.001, "P < 0.001", sprintf("P = %.3f", p_emp)))

p3r2 <- ggplot() +
  geom_ribbon(data = dens, aes(x = x, ymin = ymin, ymax = ymax, group = driver),
              fill = "grey84", colour = "grey58", linewidth = 0.3) +
  geom_segment(data = pts, aes(x = xpos, xend = xpos, y = yi - 0.10, yend = yi + 0.46, colour = class),
               linewidth = 1.6, lineend = "round") +
  geom_point(data = filter(pts, off), aes(x = xpos, y = yi + 0.18, colour = class),
             shape = 60, size = 4.5, show.legend = FALSE) +
  geom_text(data = filter(pts, off), aes(x = xpos + 0.022, y = yi + 0.30,
                                         label = sprintf("%.2f", observed)),
            hjust = 0, size = 6.2, colour = "grey25") +
  geom_text(data = pts, aes(x = X_MAX - 0.01, y = yi + 0.18, label = lab),
            hjust = 1, size = 6.8, colour = "grey25") +
  scale_colour_manual(values = CLASS_COLS, breaks = names(CLASS_COLS), name = NULL) +
  guides(colour = guide_legend(override.aes = list(linewidth = 2.2))) +
  scale_y_continuous(NULL, breaks = seq_along(ord), labels = ifelse(grepl("\n", ord), ord, str_wrap(ord, width = 20)),
                     expand = expansion(add = c(0.18, 0.75))) +
  scale_x_continuous("Variance explained (adjusted R\u00b2)",
                     limits = c(X_MIN, X_MAX), expand = c(0, 0)) +
  theme_science +
  theme(axis.text.y = element_text(size = 21, lineheight = 0.9), axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 26, margin = margin(t = 14)),
        panel.grid = element_blank(),
        legend.position = "inside", legend.position.inside = c(0.99, 1.0),
        legend.justification = c(1, 1), legend.direction = "horizontal",
        legend.text = element_text(size = 20), legend.key.width = unit(2.2, "lines"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key = element_blank())

## ---- save ----
## only the two final panels are saved, as separate files for assembly in Inkscape
ggsave(file.path(OUT, "Fig3a_occupancy_prediction.pdf"), p3e,  width = 10, height = 6)
ggsave(file.path(OUT, "Fig3b_variance_vs_chance.pdf"),   p3r2, width = 10, height = 8.5)
## final Figure 3: both panels with bold a/b tags on a compact canvas (text stays large once placed)
fig3 <- (p3e | p3r2) + plot_layout(widths = c(1.6, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 44, family = "tagbold"),
        legend.background = element_blank(), legend.box.background = element_blank())
ggsave(file.path(OUT, "Figure3_final.pdf"), fig3, width = 24, height = 9.5)
showtext_opts(dpi = 300)
ggsave(file.path(OUT, "Figure3_final.png"), fig3, width = 24, height = 9.5, dpi = 300)
showtext_opts(dpi = 96)
cat("Saved to", OUT, "\n")

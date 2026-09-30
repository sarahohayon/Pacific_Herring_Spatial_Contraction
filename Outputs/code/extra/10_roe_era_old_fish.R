## ---------------------------------------------------------------------------
## The roe era on its own terms (1975-2024). Three panels, one argument.
##
##  a  The fishery takes a far larger share of the old fish than of the stock.
##     "Only 15-20% of SSB is taken" is true and beside the point: the unfished
##     remainder is mostly ages 3 and 4.
##  b  So the population never ages. Under a TAC since 1986, biomass drifts up
##     while the share of fish at ages 5, 6 and 7+ drifts flat or down.
##  c  And spawning occupancy, which fell once in 1984, never returns.
## ---------------------------------------------------------------------------
source("00_setup.R")
suppressMessages({library(tidyr); library(patchwork); library(scales)})

a  <- read_csv(file.path(DERIVED, "annual_series.csv"), show_col_types = FALSE)
bs <- read_csv(file.path(DERIVED, "biosample_annual_by_source.csv"), show_col_types = FALSE)
tf <- read_csv(file.path(DERIVED, "test_fishery_age.csv"), show_col_types = FALSE)

d <- a %>%
  left_join(bs %>% filter(src == "Test fishery") %>% transmute(Year, share = prop_old_wt), by = "Year") %>%
  mutate(pop = catch_total_t + SSB, old_pool = share * pop,
         U     = 100 * catch_total_t / pop,
         U_old = 100 * old_removed_t / old_pool) %>%
  filter(Year >= 1975, Year <= 2024)

X <- list(scale_x_continuous(limits = c(1975, 2024), breaks = seq(1980, 2020, 10),
                             expand = expansion(mult = c(0.02, 0.02))),
          annotate("rect", xmin = 1983.5, xmax = 1986.5, ymin = -Inf, ymax = Inf,
                   fill = "grey55", alpha = 0.30),
          annotate("rect", xmin = 1977.5, xmax = 1979.5, ymin = -Inf, ymax = Inf,
                   fill = "#E0A93C", alpha = 0.25),
          geom_vline(xintercept = 1986, linetype = "dashed", colour = "grey40", linewidth = 0.6),
          theme_science,
          theme(axis.title = element_text(size = 17), axis.text = element_text(size = 15),
                plot.title = element_text(size = 17, face = "bold"),
                panel.grid.minor = element_blank()))

col_stock <- "grey55"; col_old <- "#B23A34"

## a  removal rates
pa <- ggplot(d, aes(Year)) + X +
  geom_line(aes(y = U_old), colour = col_old, linewidth = 1.6) +
  geom_line(aes(y = U),     colour = col_stock, linewidth = 1.4) +
  annotate("text", x = 2005, y = 57, hjust = 0, size = 5.6, colour = col_old, fontface = "bold",
           label = "share of the OLD FISH removed") +
  annotate("text", x = 2005, y = 8,  hjust = 0, size = 5.6, colour = "grey35", fontface = "bold",
           label = "share of the whole stock removed") +
  annotate("text", x = 1986.4, y = 66, hjust = 0, size = 5, colour = "grey35", label = "TAC era") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 70)) +
  labs(y = "Removed per year", x = NULL,
       title = "a   The fishery takes a third of the old fish each year, not a sixth of the stock")

## b  age structure of the spawning population
ages <- tf %>% filter(Year >= 1975) %>%
  transmute(Year, `age 5` = 100*prop_age5, `age 6` = 100*prop_age6, `age 7 and older` = 100*prop_age7) %>%
  pivot_longer(-Year, names_to = "age", values_to = "pct")
age_cols <- c(`age 5` = "#E8A598", `age 6` = "#C4635A", `age 7 and older` = "#7B2119")
pb <- ggplot(ages, aes(Year, pct, colour = age)) + X +
  geom_line(linewidth = 1.4) +
  scale_colour_manual(values = age_cols, name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(y = "Of the spawning population", x = NULL,
       title = "b   Biomass recovers, the fish do not get older") +
  theme(legend.position = c(0.99, 0.97), legend.justification = c(1, 1),
        legend.text = element_text(size = 15), legend.direction = "horizontal")

## c  the response
pc <- ggplot(d, aes(Year, 100*occupancy)) + X +
  geom_line(colour = "grey20", linewidth = 1.5) +
  geom_point(size = 2, colour = "grey20") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 80)) +
  labs(y = "Sections with spawn", x = "Year",
       title = "c   Spawning occupancy falls once, in 1984, and stays down")

fig <- pa / pb / pc
save_fig(fig, "Roe_era_old_fish_and_occupancy", width = 12, height = 13)

cat("\n=== means by period ===\n")
d2 <- d %>% mutate(ph = ifelse(Year <= 1983, "1977-1983", ifelse(Year <= 1986, "1984-1986", "1987-2024")))
print(as.data.frame(d2 %>% group_by(ph) %>%
  summarise(SSB_t = round(mean(SSB)), stock_removed_pct = round(mean(U, na.rm=TRUE)),
            old_removed_pct = round(mean(U_old, na.rm=TRUE)),
            occupancy_pct = round(100*mean(occupancy)), .groups = "drop")), row.names = FALSE)
cat("\nfigure: Outputs/figures/Roe_era_old_fish_and_occupancy\n")

## ---------------------------------------------------------------------------
## Figure: everything in tonnes. Whole population, the old fish inside it, and the
## old fish removed each year. Labels sit in the top right; the lines carry nothing.
## ---------------------------------------------------------------------------
col_pop <- "grey50"; col_pool <- "#B23A34"; col_rem <- "#7B2119"; col_era <- "grey45"
key <- function(y, colour, text)
  annotate("text", x = 2024, y = y, hjust = 1, size = 5.8, colour = colour,
           fontface = "bold", label = text)

p_t <- ggplot(d, aes(Year)) +
  annotate("rect", xmin = 1983.5, xmax = 1986.5, ymin = -Inf, ymax = Inf,
           fill = "grey55", alpha = 0.30) +
  geom_area(aes(y = pop), fill = "grey70", alpha = 0.28) +
  geom_line(aes(y = pop), colour = col_pop, linewidth = 1.3) +
  geom_area(aes(y = old_pool), fill = col_pool, alpha = 0.22) +
  geom_line(aes(y = old_pool), colour = col_pool, linewidth = 1.6) +
  geom_col(aes(y = old_removed_t), fill = col_rem, width = 0.6) +
  key(212000, col_pop,  "Whole population (catch + spawn index)") +
  key(198000, col_pool, "Old fish in the population (age 5+)") +
  key(184000, col_rem,  "Old fish removed by the fishery") +
  key(170000, "grey40", "Occupancy collapse window (1984-1986)") +
  scale_y_continuous(labels = label_number(big.mark = ","), limits = c(0, 222000),
                     expand = expansion(mult = c(0, 0.01))) +
  scale_x_continuous(limits = c(1975, 2024), breaks = seq(1980, 2020, 10),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Year", y = "Tonnes") +
  theme_science +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16),
        panel.grid.minor = element_blank())
save_fig(p_t, "Old_fish_tonnes", width = 12, height = 7.5)

## ---------------------------------------------------------------------------
## Figure: two exploitation rates, same catch, different denominators.
##   overall   catch / (catch + spawn index)
##   old fish  old fish removed / old fish in the population
## Reference lines are DFO's harvest rate: 20% from 1983, and the 14% target now
## used for the Strait under the MSE management procedure.
## ---------------------------------------------------------------------------
dp <- d %>% mutate(pct_catch      = 100 * catch_total_t / pop,
                   pct_old_of_old = 100 * old_removed_t / old_pool)
col_catch <- "#2C3E70"; col_ofo <- "#B23A34"
keyp <- function(y, colour, text)
  annotate("text", x = 2024, y = y, hjust = 1, size = 5.8, colour = colour,
           fontface = "bold", label = text)

p_p <- ggplot(dp, aes(Year)) +
  annotate("rect", xmin = 1983.5, xmax = 1986.5, ymin = -Inf, ymax = Inf,
           fill = "grey55", alpha = 0.30) +
  ## DFO harvest rate. The TAC record we hold begins in 1986, so the rule is drawn
  ## from there; 14% is the target now applied to the Strait under the MSE procedure.
  annotate("segment", x = 1986, xend = 2022, y = 20, yend = 20,
           linetype = "dashed", colour = "grey35", linewidth = 0.8) +
  annotate("segment", x = 2022, xend = 2024, y = 14, yend = 14,
           linetype = "dashed", colour = "grey35", linewidth = 0.8) +
  annotate("text", x = 2019, y = 21.9, hjust = 1, size = 4.9, colour = "grey35", label = "20%") +
  annotate("text", x = 2024, y = 12.1, hjust = 1, size = 4.9, colour = "grey35", label = "14%") +
  geom_line(aes(y = pct_old_of_old), colour = col_ofo,   linewidth = 1.9) +
  geom_line(aes(y = pct_catch),      colour = col_catch, linewidth = 1.5) +
  keyp(72, col_ofo,   "Old fish exploitation rate (old fish removed / old fish in population)") +
  keyp(66, col_catch, "Overall exploitation rate (catch / catch + SSB)") +
  keyp(60, "grey40",  "Occupancy collapse window (1984-1986)") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 78),
                     breaks = seq(0, 70, 10), expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(limits = c(1975, 2024), breaks = seq(1980, 2020, 10),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Year", y = "Exploitation rate") +
  theme_science +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16),
        panel.grid.minor = element_blank())
save_fig(p_p, "Old_fish_percent_of_population", width = 12, height = 7.5)
cat("figures: Old_fish_tonnes | Old_fish_percent_of_population\n")

## ---------------------------------------------------------------------------
## The same age classes, but as tonnes in the spawning population rather than as
## shares of the old-fish pool. Within-pool percentages hide half the loss: the
## pool shrank and the fish inside it got younger.
## Tonnes at age = (share of sample weight at that age) x (catch + spawn index).
## ---------------------------------------------------------------------------
w_age <- read_csv(F_BIO, show_col_types = FALSE, name_repair = "minimal", guess_max = 1e5)
names(w_age) <- BIO_COLS[seq_len(ncol(w_age))]
w_age <- w_age %>%
  mutate(age = suppressWarnings(as.numeric(age)),
         weight_g = suppressWarnings(as.numeric(weight_g)),
         Year = suppressWarnings(as.numeric(year))) %>%
  mutate(Year = ifelse(Year < 30, Year + 2000, ifelse(Year < 100, Year + 1900, Year))) %>%
  filter(!is.na(age), age >= 1, !is.na(weight_g), weight_g > 0,
         stat_area %in% BIO_AREAS, grepl("test", source, ignore.case = TRUE))

## Age 3 and older only: the spawn index is spawning biomass, and herring mature at
## about age 3, so ages 1 and 2 are not part of it. Shares are renormalised within
## ages 3+ so the bands sum to the spawning population.
age_share <- w_age %>%
  filter(age >= 3) %>%
  mutate(age_c = pmin(age, 9)) %>%
  group_by(Year) %>% mutate(tot_w = sum(weight_g)) %>%
  group_by(Year, age_c) %>% summarise(share_w = sum(weight_g) / first(tot_w), .groups = "drop")

tons_at_age <- age_share %>%
  left_join(d %>% select(Year, pop), by = "Year") %>%
  filter(!is.na(pop)) %>%
  mutate(tonnes = share_w * pop,
         age_lab = factor(ifelse(age_c == 9, "9+", as.character(age_c)),
                          levels = c("5", "6", "7", "8", "9+")),
         phase = case_when(Year <= 1983 ~ "Pre-collapse\n(1977-1983)",
                           Year <= 1986 ~ NA_character_,
                           TRUE         ~ "Post-collapse\n(1987-2024)")) %>%
  filter(!is.na(phase)) %>%
  mutate(phase = factor(phase, levels = c("Pre-collapse\n(1977-1983)", "Post-collapse\n(1987-2024)")))

## annual series of tonnes at age, coloured as in Figure 2c of the manuscript
age_palette <- c("3" = "#40916C", "4" = "#74C69D", "5" = "#B7E4C7", "6" = "#FFD166",
                 "7" = "#F4A261", "8" = "#E76F51", "9+" = "#C0392B")

series <- age_share %>%
  left_join(d %>% select(Year, pop), by = "Year") %>%
  filter(!is.na(pop)) %>%
  mutate(tonnes = share_w * pop / 1000,
         age_lab = factor(ifelse(age_c == 9, "9+", as.character(age_c)),
                          levels = names(age_palette)))
save_tab(series %>% transmute(Year, age = age_lab, kt = round(tonnes, 2)), "tonnes_at_age_series")

## stacked: the height of the band is the whole old-fish pool, the bands within it
## are the age classes that make it up
p_ton_age <- ggplot(series, aes(Year, tonnes, fill = age_lab)) +
  geom_area(colour = "white", linewidth = 0.25) +
  annotate("rect", xmin = 1983.5, xmax = 1986.5, ymin = -Inf, ymax = Inf,
           fill = "grey25", alpha = 0.30) +
  annotate("text", x = 1987.5, y = 158, hjust = 0, size = 4.8, colour = "grey30",
           lineheight = 0.95, label = "occupancy collapse\nwindow (1984-1986)") +
  scale_fill_manual(values = age_palette, name = "Age",
                    guide = guide_legend(reverse = TRUE, ncol = 2)) +
  scale_y_continuous(labels = label_number(big.mark = ","),
                     expand = expansion(mult = c(0.01, 0.05))) +
  scale_x_continuous(breaks = seq(1980, 2020, 10),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Year", y = "Spawning population (\u00d7 1,000 tonnes)") +
  theme_science +
  theme(legend.position = c(0.005, 0.995), legend.justification = c(0, 1),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12),
        legend.key.size = unit(0.42, "cm"), legend.spacing.y = unit(1, "pt"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        axis.title = element_text(size = 18), axis.text = element_text(size = 16),
        panel.grid.minor = element_blank())
save_fig(p_ton_age, "Old_fish_tonnes_at_age", width = 12, height = 7.5)
cat("figure: Outputs/figures/Old_fish_tonnes_at_age\n")

## ---------------------------------------------------------------------------
## Panel d: the same stack rescaled so each year fills 100%. Absolute tonnage is
## gone, composition is all that is left, and the warm band's collapse is plain.
## ---------------------------------------------------------------------------
p_prop <- ggplot(series, aes(Year, tonnes, fill = age_lab)) +
  geom_area(position = "fill", colour = "white", linewidth = 0.25) +
  annotate("rect", xmin = 1983.5, xmax = 1986.5, ymin = -Inf, ymax = Inf,
           fill = "grey25", alpha = 0.30) +
  scale_fill_manual(values = age_palette, name = "Age",
                    guide = guide_legend(reverse = TRUE, ncol = 2)) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(limits = c(1975, 2024), breaks = seq(1980, 2020, 10),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Year", y = "Share of the spawning population") +
  theme_science +
  theme(legend.position = c(0.005, 0.995), legend.justification = c(0, 1),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12),
        legend.key.size = unit(0.42, "cm"), legend.spacing.y = unit(1, "pt"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        axis.title = element_text(size = 18), axis.text = element_text(size = 16),
        panel.grid.minor = element_blank())
save_fig(p_prop, "Old_fish_proportional", width = 12, height = 6)

## ---------------------------------------------------------------------------
## The three panels as one figure, shared x axis, no tags (added in Inkscape).
##   a  tonnes: the population, the old fish in it, and what was removed
##   b  the two exploitation rates the same catch implies
##   c  the spawning population by age class
## ---------------------------------------------------------------------------
strip_x <- theme(axis.title.x = element_blank())
fig_all <- (p_t + strip_x) / (p_p + strip_x) / (p_ton_age + strip_x) / p_prop +
  plot_layout(heights = c(1, 1, 1, 0.85))
save_fig(fig_all, "Old_fish_summary", width = 12, height = 23)
cat("figure: Outputs/figures/Old_fish_summary\n")

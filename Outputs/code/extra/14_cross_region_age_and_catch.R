## ---------------------------------------------------------------------------
## Coast-wide version of the Strait figures: the spawning biomass of each stock
## broken into age classes, and the catch that was taken from it.
##
##   age composition  test-fishery samples, by weight, ages 3 to 9+, 1976 onward
##                    (1975 excluded, see 13_cross_region_occupancy.R)
##   biomass          spawn index, surveyed records only
##   catch            DFO catch by section, 1950-2016, mapped to region by area
## ---------------------------------------------------------------------------
source("00_setup.R")
suppressMessages({library(data.table); library(tidyr); library(patchwork); library(scales)})

BIO_DB <- "/Users/sarah/Documents/Postdoc/Canada Research/Pacific herring/CLAUDE CODE/Data/Herring Biosample Database"
CATCH_F <- "/Users/sarah/Documents/Postdoc/Canada Research/Pacific herring/CLAUDE CODE/Data/DFO_Herring_Catch_Data.csv"
REG <- c(SoG = "Strait_of_Georgia.csv", WCVI = "West_Coast_Vancouver_Island.csv",
         CC = "Central_Coast.csv", PRD = "Prince_Rupert_District.csv", HG = "Haida_Gwaii.csv")
REG_NAME <- c(SoG = "Strait of Georgia", WCVI = "West Coast Vancouver Island",
              CC = "Central Coast", PRD = "Prince Rupert District", HG = "Haida Gwaii")

## ---- biomass and the area-to-region key -----------------------------------
sp <- read_csv(F_SPAWN, show_col_types = FALSE, guess_max = 1e5) %>% keep_surveyed() %>%
  mutate(bio = rowSums(cbind(Surface, Understory, Macrocystis), na.rm = TRUE))
key <- sp %>% filter(!is.na(Region)) %>% count(StatisticalArea, Region) %>%
  group_by(StatisticalArea) %>% slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(StatisticalArea, reg = Region)
ssb <- sp %>% left_join(key, by = "StatisticalArea") %>%
  filter(reg %in% names(REG), Year %in% 1951:2024) %>%
  group_by(reg, Year) %>% summarise(SSB = sum(bio, na.rm = TRUE), .groups = "drop")

## ---- age composition by weight --------------------------------------------
age_w <- bind_rows(lapply(names(REG), function(rg) {
  dt <- fread(file.path(BIO_DB, REG[rg]), select = c(2, 3, 15, 22, 26),
              col.names = c("sample", "year", "source", "weight_g", "age"), showProgress = FALSE)
  as_tibble(dt) %>%
    mutate(age = suppressWarnings(as.numeric(age)), weight_g = suppressWarnings(as.numeric(weight_g)),
           Year = suppressWarnings(as.numeric(year)),
           Year = ifelse(Year < 30, Year + 2000, ifelse(Year < 100, Year + 1900, Year))) %>%
    filter(source == "Test Fishery", !is.na(age), age >= 3, !is.na(weight_g), weight_g > 0,
           Year >= 1976) %>%
    mutate(age_c = pmin(age, 9)) %>%
    group_by(Year) %>% mutate(tot = sum(weight_g)) %>%
    group_by(Year, age_c) %>% summarise(share = sum(weight_g)/first(tot), .groups = "drop") %>%
    mutate(reg = rg)
}))

age_pal <- c("3" = "#40916C", "4" = "#74C69D", "5" = "#B7E4C7", "6" = "#FFD166",
             "7" = "#F4A261", "8" = "#E76F51", "9+" = "#C0392B")
ser <- age_w %>% left_join(ssb, by = c("reg", "Year")) %>% filter(!is.na(SSB)) %>%
  mutate(kt = share * SSB / 1000,
         age_lab = factor(ifelse(age_c == 9, "9+", as.character(age_c)), levels = names(age_pal)),
         region = factor(REG_NAME[reg], levels = unname(REG_NAME)))
write_csv(ser %>% transmute(region, Year, age = age_lab, kt = round(kt, 2)),
          file.path(DERIVED, "cross_region_ssb_at_age.csv"))

p_age <- ggplot(ser, aes(Year, kt, fill = age_lab)) +
  geom_area(colour = "white", linewidth = 0.2) +
  facet_wrap(~ region, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = age_pal, name = "Age", guide = guide_legend(reverse = TRUE, ncol = 2)) +
  scale_x_continuous(limits = c(1976, 2024), breaks = seq(1980, 2020, 10)) +
  scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Year", y = "Spawning biomass (× 1,000 tonnes)") +
  theme_science +
  theme(strip.background = element_blank(), strip.text = element_text(size = 15, face = "bold"),
        legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(size = 13), legend.title = element_text(size = 13),
        legend.key.size = unit(0.45, "cm"),
        axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        panel.grid.minor = element_blank())
save_fig(p_age, "CrossRegion_SSB_at_age", width = 11, height = 15)

## ---- catch -----------------------------------------------------------------
ct <- read_csv(CATCH_F, show_col_types = FALSE, guess_max = 1e5) %>%
  mutate(Year = suppressWarnings(as.numeric(Year)), Section = suppressWarnings(as.numeric(Section)),
         across(c(Total_Catch_tonnes, Gillnet, Seine, Trawl, SOK), ~suppressWarnings(as.numeric(.x)))) %>%
  filter(!is.na(Year), !is.na(Section), Year >= 1951) %>%
  ## the spawn index stores the area as a zero-padded string, so match that
  mutate(StatisticalArea = sprintf("%02d", Section %/% 10)) %>%
  left_join(key, by = "StatisticalArea") %>% filter(reg %in% names(REG))

catch_gear <- ct %>%
  select(reg, Year, Gillnet, Seine, Trawl, SOK) %>%
  pivot_longer(c(Gillnet, Seine, Trawl, SOK), names_to = "gear", values_to = "t") %>%
  filter(!is.na(t)) %>% group_by(reg, Year, gear) %>% summarise(kt = sum(t)/1000, .groups = "drop") %>%
  mutate(region = factor(REG_NAME[reg], levels = unname(REG_NAME)),
         gear = factor(gear, levels = c("Seine", "Gillnet", "Trawl", "SOK")))
write_csv(catch_gear %>% transmute(region, Year, gear, kt = round(kt, 2)),
          file.path(DERIVED, "cross_region_catch.csv"))

gear_pal <- c(Seine = "#2C3E70", Gillnet = "#B23A34", Trawl = "#8C8C8C", SOK = "#E0A93C")
p_catch <- ggplot(catch_gear, aes(Year, kt, fill = gear)) +
  geom_col(width = 0.9) +
  facet_wrap(~ region, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = gear_pal, name = NULL) +
  scale_x_continuous(limits = c(1951, 2017), breaks = seq(1960, 2010, 10)) +
  scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Year", y = "Catch (× 1,000 tonnes)") +
  theme_science +
  theme(strip.background = element_blank(), strip.text = element_text(size = 15, face = "bold"),
        legend.position = "top", legend.direction = "horizontal",
        legend.text = element_text(size = 13), legend.key.size = unit(0.45, "cm"),
        axis.title = element_text(size = 16), axis.text = element_text(size = 14),
        panel.grid.minor = element_blank())
save_fig(p_catch, "CrossRegion_catch", width = 11, height = 15)

cat("\n=== proportion aged 5+ in the spawning biomass, by region and period ===\n")
print(as.data.frame(ser %>% mutate(old = age_c >= 5, per = ifelse(Year <= 1983, "1976-1983",
        ifelse(Year <= 1990, "1984-1990", "1991-2024"))) %>%
  group_by(region, per) %>% summarise(pct_old = round(100*sum(share[old])/sum(share)), .groups = "drop") %>%
  pivot_wider(names_from = per, values_from = pct_old)), row.names = FALSE)
cat("\n=== catch by era (kt/yr) ===\n")
print(as.data.frame(catch_gear %>% mutate(era = ifelse(Year <= 1967, "1951-1967",
        ifelse(Year <= 1971, "closure", ifelse(Year <= 1983, "1972-1983", "1984-2016")))) %>%
  group_by(region, era) %>% summarise(kt_yr = round(sum(kt)/n_distinct(Year), 1), .groups = "drop") %>%
  pivot_wider(names_from = era, values_from = kt_yr)), row.names = FALSE)
cat("\nfigures: CrossRegion_SSB_at_age | CrossRegion_catch\n")

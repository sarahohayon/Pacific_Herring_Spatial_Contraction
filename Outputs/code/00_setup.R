## =============================================================================
##  00_setup.R - shared paths, constants and helpers
##  Pacific herring, Strait of Georgia: spatial collapse & collective memory
##
##  Sourced by every 0X_*.Rmd in this folder. Nothing here fits a model or
##  writes a result; it only defines where things live and the few constants
##  that must be identical across all analyses.
## =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

## ---- paths ------------------------------------------------------------------
## All paths are relative to this R_code/ folder, so the whole Aligned_2026_08_06
## directory can be moved or shared without editing anything.

if (!exists("PROJ_ROOT")) {
  PROJ_ROOT <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."),
                             mustWork = FALSE)
}
## When knitted, knitr sets the working directory to the Rmd's folder (R_code/),
## so the simple relative path below is the one that actually applies.
if (!dir.exists(file.path(PROJ_ROOT, "Raw_data"))) PROJ_ROOT <- ".."
if (!dir.exists(file.path(PROJ_ROOT, "Raw_data"))) PROJ_ROOT <- file.path("..", "..")  # code lives in Outputs/code

RAW     <- file.path(PROJ_ROOT, "Raw_data")
DERIVED <- file.path(PROJ_ROOT, "Outputs", "derived")
FIGS    <- file.path(PROJ_ROOT, "Outputs", "figures")
TABS    <- file.path(PROJ_ROOT, "Outputs", "tables")
for (p in c(DERIVED, FIGS, TABS)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

## ---- raw data files ---------------------------------------------------------
F_SPAWN   <- file.path(RAW, "Pacific_herring_spawn_index_data_2025_EN.csv")
F_BIO     <- file.path(RAW, "Biosample_Strait_of_Georgia.csv")
F_CATCH   <- file.path(RAW, "SOG_herring_catch_all_fishing_types.csv")
F_CATCH2  <- file.path(RAW, "SOG_herring_catch_all_fishing_types_2026_03_11.csv")
F_TAC     <- file.path(RAW, "SOG_TAC.csv")
F_SST     <- file.path(RAW, "Lighthouse_Stations_SST_Combined.csv")
F_ONI     <- file.path(RAW, "ONI_index_1950_2025.csv")
F_PDO     <- file.path(RAW, "PDO_1854_2025.csv")
F_VALUE   <- file.path(RAW, "Pacific herring catch and landed value dfo_INFLATION.xlsx")
F_SECT    <- file.path(RAW, "sections", "SectionsIntegrated.shp")

## ---- analysis constants -----------------------------------------------------
## These four decisions propagate through every analysis and are stated in
## Methods. Change them here and nowhere else.

OCC_AREAS  <- c(13, 14, 15, 16, 17, 18, 19, 28, 29) # SoG statistical areas: occupancy domain
SOG_REGION <- "SoG"  # DFO stock assessment region in the spawn index. Statistical areas 13 and 29
                     # extend beyond the SoG stock (Area 13 sections 131, 133, 134, 136 = Johnstone
                     # Strait/Discovery Islands; Area 29 section 293), which DFO leaves unassigned.
                     # Spawn records are kept only if Region == SOG_REGION.
BIO_AREAS  <- c(13, 14, 15, 17, 18)                 # areas with sufficient biosample coverage
MIN_ACTIVE <- 6      # a section is "historically used" if spawn recorded in >= 6 distinct years
MIN_FISH   <- 30     # a biosample is used only if it contains >= 30 aged fish
## Test fishery: all months, seine sets, from 1977. 1975 holds February samples
## only and no samples exist for 1976, so the series starts in 1977.
TF_FIRST_YEAR <- 1977
TF_GEARS      <- c("Seine", "Other seine")
keep_tf <- function(source, year, gear)
  source != "Test Fishery" | (year >= TF_FIRST_YEAR & gear %in% TF_GEARS)
OLD_AGE    <- 5      # "old fish" = age >= 5 (repeat spawners; age-3 are first-time spawners)
YRS        <- 1951:2024
BREAK_YEAR <- 1984   # collapse boundary: pre <= 1983, post >= 1984 (largest single-year drop)

## Spawn presence rule. A spawn record counts as observed spawning, including records for
## which DFO gives no biomass estimate (Surface, Understory and Macrocystis all NA, e.g. every
## Area 28 / Howe Sound record), EXCEPT records whose Method is "Incomplete": that category
## first appears in the 1990s and grows through the 2020s, so counting it would add presence
## only after the collapse (a reporting change, not a survey). All Incomplete records lack biomass.
keep_surveyed <- function(df) {
  dplyr::filter(df, !(is.na(Surface) & is.na(Understory) & is.na(Macrocystis) &
                        Method %in% "Incomplete"))
}

## Three-phase split used for the demographic contrasts (Results 2.2)
PHASE3 <- function(y) dplyr::case_when(y <= 1983 ~ "Pre-collapse",
                                       y <= 1986 ~ "Collapse",
                                       TRUE      ~ "Post-collapse")

## ---- helpers ----------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

## DFO biosamples store a 2-digit year; century break at 29/30.
recode_two_digit_year <- function(y) {
  y <- suppressWarnings(as.numeric(y))
  ifelse(is.na(y), NA_real_, ifelse(y < 100, ifelse(y <= 29, 2000 + y, 1900 + y), y))
}

## Column names of the raw DFO biosample extract (the file ships without a header
## row that matches these, so they are assigned positionally).
BIO_COLS <- c("season","sample_number","year","month","day","representative_set",
  "stock_assessment_region","stat_area","section","location_code","location",
  "latitude","longitude","source_code","source","gear_code","gear",
  "preservation_method_code","preservation_method","fish_number","length_mm","weight_g",
  "sex","maturity_code","maturity","age","dual_aged","gonad_length","gonad_weight")

## Catch gear groups, matched to the three gear groups in the biosample database.
gear_group <- function(g) dplyr::case_when(
  g %in% c("Seine","Salmon Seine","Other seine") ~ "Seine",
  g %in% c("Gillnet")                            ~ "Gillnet",
  g %in% c("Trawl","Other Trawl","Herring trawl","Other trawl") ~ "Trawl",
  TRUE ~ "Other")

## p-value formatting used in every figure annotation and table
fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p))

## ---- shared plotting theme --------------------------------------------------
## Two themes live here on purpose:
##   theme_herring() - compact, for multi-panel supplementary figures
##   theme_science   - the main-figure style (large type, Helvetica).

col_occ <- "#555555"; col_old <- "#D73027"; col_sst <- "#3461A8"

## Management-area palette. Identical to cluster_cols_9 in the exploratory Rmd.
area_cols <- c("13"="#5B8C5A","14"="#E07B3F","15"="#7EB5D6","16"="#C03B3B",
               "17"="#3A9E8F","18"="#3461A8","19"="#8B6BAE","28"="#2D6A4F","29"="#1A2F6B")

## Fishery x gear categories (Fig. 2a) and biosample sources (Fig. 2b).
fishery_cols <- c("Reduction"     = "#7B5B9D",
                  "Roe - Seine"   = "#E8894A",
                  "Roe - Gillnet" = "#B23A34",
                  "Food and bait" = "#5B8C5A")

src_cols <- c("Reduction"     = "#7B5B9D",
              "Roe seine"     = "#E8894A",
              "Roe gillnet"   = "#B23A34",
              "Food and bait" = "#5B8C5A",
              "Test fishery"  = "#3461A8")

## Era colours for the hysteresis panel: light grey pre-collapse, dark post.
PHASE_LEVELS <- c("Pre-collapse (1951-1983)", "Post-collapse (1984-2024)")
phase_fill   <- setNames(c("grey85", "grey35"), PHASE_LEVELS)
phase_edge   <- setNames(c("grey45", "grey15"), PHASE_LEVELS)
phase_hull   <- setNames(c("grey70", "grey30"), PHASE_LEVELS)

theme_herring <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(colour = "grey35", size = base_size - 1),
      axis.text     = ggplot2::element_text(colour = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text    = ggplot2::element_text(face = "bold", hjust = 0))
}

## Helvetica if it can be registered, otherwise the device default. Wrapped so
## the pipeline still knits on a machine without showtext or without the font.
HERRING_FONT <- ""
if (requireNamespace("showtext", quietly = TRUE) &&
    requireNamespace("sysfonts", quietly = TRUE)) {
  .ok <- tryCatch({
    if (!"helvetica" %in% sysfonts::font_families()) {
      .hv <- "/System/Library/Fonts/Helvetica.ttc"
      if (file.exists(.hv)) sysfonts::font_add(family = "helvetica", regular = .hv)
    }
    "helvetica" %in% sysfonts::font_families()
  }, error = function(e) FALSE)
  if (isTRUE(.ok)) { showtext::showtext_auto(); HERRING_FONT <- "helvetica" }
}

theme_science <- ggplot2::theme_classic(base_size = 20, base_family = HERRING_FONT) +
  ggplot2::theme(
    axis.title.x      = ggplot2::element_text(size = 20, margin = ggplot2::margin(t = 14)),
    axis.title.y      = ggplot2::element_text(size = 20, margin = ggplot2::margin(r = 14)),
    axis.text         = ggplot2::element_text(size = 20, colour = "black"),
    axis.ticks        = ggplot2::element_line(linewidth = 0.6),
    axis.line         = ggplot2::element_line(linewidth = 0.5),
    axis.ticks.length = ggplot2::unit(5, "pt"),
    plot.title        = ggplot2::element_text(size = 24, face = "plain"),
    legend.text       = ggplot2::element_text(size = 18),
    legend.title      = ggplot2::element_text(size = 18),
    plot.margin       = ggplot2::margin(15, 15, 15, 15, "pt")
  )

scale_x_decades <- ggplot2::scale_x_continuous(
  limits = c(1950, 2025),
  breaks = seq(1950, 2020, by = 10),
  expand = ggplot2::expansion(mult = c(0.01, 0.02)))

## common_panel suppresses the legend (most main-text panels carry none).
## Use common_panel_legend where a legend is wanted; it enlarges the keys so
## they stay readable against the small plotting symbols.
common_panel <- list(scale_x_decades, theme_science,
                     ggplot2::theme(legend.position = "none"))

common_panel_legend <- function(pos = c(0.01, 0.01), just = c(0, 0),
                                dir = "horizontal") {
  list(scale_x_decades, theme_science,
       ggplot2::theme(
         legend.position      = pos,
         legend.justification = just,
         legend.direction     = dir,
         legend.key.size      = ggplot2::unit(1.4, "lines"),
         legend.background    = ggplot2::element_rect(
           fill = scales::alpha("white", 0.65), colour = NA)))
}

## Shading annotations shared by every time-series panel:
## red band = 1968-1971 fishery closure; grey = post-collapse era.
SHADE_START <- BREAK_YEAR
SHADE_END   <- 2025
shade_df    <- tibble::tibble(xmin = SHADE_START, xmax = SHADE_END,
                              ymin = -Inf, ymax = Inf)

band_closure  <- function() ggplot2::annotate("rect", xmin = 1968, xmax = 1971,
                                              ymin = -Inf, ymax = Inf,
                                              fill = "#f4a3a3", alpha = 0.45)
band_postcoll <- function(xmax = SHADE_END) ggplot2::annotate("rect", xmin = BREAK_YEAR, xmax = xmax,
                                                         ymin = -Inf, ymax = Inf,
                                                         fill = "grey85", alpha = 0.45)

save_fig <- function(plot, name, width = 9, height = 6, dpi = 300) {
  ## showtext renders text at 96 dpi by default; match the PNG resolution so PNG text
  ## is the same size as in the PDF, then restore the default for the PDF device.
  if (nzchar(HERRING_FONT)) showtext::showtext_opts(dpi = dpi)
  ggplot2::ggsave(file.path(FIGS, paste0(name, ".png")), plot,
                  width = width, height = height, dpi = dpi, bg = "white")
  if (nzchar(HERRING_FONT)) showtext::showtext_opts(dpi = 96)
  ggplot2::ggsave(file.path(FIGS, paste0(name, ".pdf")), plot,
                  width = width, height = height, bg = "white")
  invisible(plot)
}

save_tab <- function(x, name) {
  readr::write_csv(x, file.path(TABS, paste0(name, ".csv")))
  invisible(x)
}

message("00_setup.R loaded | RAW = ", normalizePath(RAW, mustWork = FALSE))

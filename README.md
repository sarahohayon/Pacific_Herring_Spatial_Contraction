# Biomass recovery masks collective memory loss and spatial collapse in Pacific herring

Analysis code for the Strait of Georgia Pacific herring study (Ohayon, Dingwall & Bates).
Every statistic, table and figure in the paper and its supplement is produced by the files
here, starting from the raw DFO extracts. 

## Layout

```
Outputs/code/          the analysis: 00_setup.R, 01-09, and the figure scripts
Outputs/derived/       series built by 01 and read by everything else
Outputs/tables/        one CSV per supplementary table, plus supporting tables
Outputs/figures/       figure panels drawn by 01-08 (Figs. 1a,e,f, 2a, 4a-d) and Figs. S1-S5
Figures/Main/          Figures 1-4 as published, and the panels drawn by the Fig* scripts
Raw_data/              input data (see Data below)
run_all.R              runs the whole analysis in order
```

## How to run

```bash
cd Outputs/code && Rscript ../../run_all.R
```

`01` runs first: it builds the derived series that everything else reads. The commercial catch
and total allowable catch (TAC) data are not in this repository (see Data). Without them, `01`,
the catch panels of Figures 1 and 2 (`03`, `Fig2bcd_manuscript_style.R`) and the catch-based parts
of `04` (Figure 2a and removals at age) cannot be rerun; everything else runs from the series in
`Outputs/derived/`.

A full run takes about 15 minutes, mostly in `01` (a 100 MB biosample file) and `05` and the
randomization scripts (thousands of model fits). If `rmarkdown::render` cannot find pandoc
outside RStudio:

```bash
export RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64"
```

## What each file produces

| File | Produces |
|---|---|
| `00_setup.R` | paths, analysis constants, the presence rule, the test-fishery rule, themes and palettes. Read by every other file. |
| `01_data_preparation.Rmd` | all derived series in `Outputs/derived/` |
| `02_contraction_metrics.Rmd` | eleven contraction metrics and the PCA behind the choice of response (**Table S6**) |
| `03_biomass_occupancy_decoupling.Rmd` | biomass-occupancy decoupling and spawning phenology (**Fig. 1a, e, f**, **Fig. S1**, **Tables S1, S7**) |
| `04_demography_selectivity.Rmd` | fishery selectivity and age structure (**Fig. 2a**, **Fig. S3**); loss of older fish versus recruitment dilution, mortality and removals at age (**Table S10**) |
| `05_occupancy_drivers_GAMM.Rmd` | models of spawning occupancy, robustness checks and subset refits (**Figs. S4, S5**, **Tables S2-S5, S7**) |
| `06_spatial_asymmetry_refuge.Rmd` | spatial collapse, historical exploitation and the Area 14 refuge (**Fig. 4**) |
| `08_roe_gillnet_age_composition.Rmd` | age composition of the roe-gillnet catch against the test fishery (**Fig. S2**) |
| `09_prop_old_by_area.R` | older fish in Area 14 versus the peripheral areas (Results) |
| `Fig1_bcd_SoG_pipeline.R` | Figure 1 b, c, d |
| `Fig2bcd_manuscript_style.R` | Figure 2 b, c, d |
| `Fig2c_length_at_age_tests.R` | length-at-age tests quoted in the Results |
| `Fig3b_null_draws.R`, `Fig3b_null_draws_extra.R` | randomization nulls for Figure 3b and **Table S8b** |
| `Fig3b_sst_increment_test.R` | whether SST adds to the older-fish model (**Table S8c**) |
| `Fig3_randomisation_test.R` | slope randomizations and the gear-era model (**Table S8a**) |
| `Fig3_manuscript_style.R` | Figure 3 a, b |

Main-text figures are drawn as single panels and assembled into the final layout by hand;
the assembled Figures 1-4 are in `Figures/Main/`. The HTML report next to each `.Rmd` holds
its printed output.

## Decisions used throughout

They are set at the top of `00_setup.R` and described in the Methods.

* **Presence rule** (`keep_surveyed`): a spawn record counts as presence if it was surveyed,
  including records without a biomass estimate; records DFO classes as *Incomplete* are
  excluded, because that category first appears in the 1990s and would add presence only after
  the collapse. Table S7 repeats the results under the two alternative rules.
* **Historically used sections** (`MIN_ACTIVE = 6`): spawn recorded in at least six years,
  giving 19 sections.
* **Biosamples** (`MIN_FISH = 30`): a sample is used only if it holds at least 30 aged fish.
* **Older fish** (`OLD_AGE = 5`): age 5 and older.
* **Test fishery** (`keep_tf`): all months, seine sets, from 1977. The 1975 samples were
  collected in February only, and there are no samples for 1976.
* **Collapse boundary** (`BREAK_YEAR = 1984`): the year of the largest single-year decline in
  occupancy.
* **Map rule** (in `06`): the two-period maps show sites used in at least 2 years of any
  5-year window.

`01` writes two site-by-year files. `site_year_spawn.csv` keeps only records with coordinates
and is used for anything placed in space (contraction metrics, maps, per-area change).
`site_year_spawn_all.csv` keeps every record and is used for site occupancy histories, so a
site is never dropped from its own history for lacking a position.

## Data

| File | Source | Included |
|---|---|---|
| `Pacific_herring_spawn_index_data_2025_EN.csv` | DFO Pacific Herring Spawn Index, 1951-2025 | yes |
| `Lighthouse_Stations_SST_Combined.csv` | DFO BC lighthouse sea-surface temperature | yes |
| `ONI_index_1950_2025.csv`, `PDO_1854_2025.csv` | NOAA climate indices | yes |
| `Pacific herring catch and landed value dfo_INFLATION.xlsx` | DFO landed value, deflated with the CPI | yes |
| `SOG_TAC.csv` | total allowable catch, provided by DFO | no, available from DFO on request |
| `sections/SectionsIntegrated.shp` | DFO herring section boundaries | yes |
| `Biosample_Strait_of_Georgia.csv.gz` | DFO herring biosample database, Strait of Georgia, 1946-2024 (gzip-compressed; read directly) | yes |
| `SOG_herring_catch_all_fishing_types*.csv` | DFO commercial catch (reduction, food and bait, roe) by season, area, fishery and gear, provided by DFO | no, available from DFO on request |

Commercial catch and TAC data, and the detailed catch series derived from them
(`catch_decomposition.csv`), are available from DFO on request. Roe-fishery catch is also
published by DFO on the Open Government Portal. Annual catch totals and the older-fish removal
series used in the models are included in `Outputs/derived/annual_series.csv`.

The biosample extract is published by DFO as a CSV and stored here gzip-compressed; it has no
usable header row, so column names are assigned from `BIO_COLS` in
`00_setup.R`, and two-digit years are converted to calendar years.

## Requirements

R 4.5 or later, with tidyverse, mgcv, nlme, lme4, segmented, strucchange, sf, adehabitatHR,
spdep, geosphere, ineq, rnaturalearth (and rnaturalearthhires), gratia, patchwork, ggrepel,
showtext, scales, knitr and rmarkdown.

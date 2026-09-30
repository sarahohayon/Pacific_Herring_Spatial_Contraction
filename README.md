# Biomass recovery masks collective memory loss and spatial collapse in Pacific herring

Analysis code for the Strait of Georgia Pacific herring study (Ohayon, Dingwall & Bates).
Every statistic, table and figure in the paper and its supplement is produced by the files
here, starting from the raw DFO extracts. Nothing reads a cached result.

## Layout

```
Outputs/code/          the analysis: 00_setup.R, 01-09, and the figure scripts
Outputs/code/extra/    analyses explored during the study but not reported in the paper
Outputs/derived/       series built by 01 and read by everything else
Outputs/tables/        one CSV per supplementary table, plus supporting tables
Outputs/figures/       every generated figure, PDF and PNG
Figures/               final figures as used in the paper
Raw_data/              input data (see Data below)
run_all.R              runs the whole analysis in order
```

## How to run

```bash
cd Outputs/code && Rscript ../../run_all.R
```

`01` runs first: it builds the derived series that everything else reads. It is also the only
step that needs the two DFO extracts that are not in this repository; without them, start at
`02`, which runs from the series already in `Outputs/derived/`.

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
| `03_biomass_occupancy_decoupling.Rmd` | biomass-occupancy decoupling and spawning phenology (**Fig. 1**, **Fig. S1**, **Tables S1, S7**) |
| `04_demography_selectivity.Rmd` | fishery selectivity and age structure (**Fig. 2**, **Fig. S3**); loss of older fish versus recruitment dilution, mortality and removals at age (**Table S10**) |
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

Main-text figures are saved both as single panels (for assembling the final layout) and as a
combined figure. The HTML report next to each `.Rmd` holds its printed output.

The scripts in `extra/` (a roe-era summary, the comparison with other BC stocks, age-specific
exploitation and supplementary catch plots) were used to check the results but are not
reported in the paper; some of them need DFO data that are not in this repository.

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
| `SOG_TAC.csv` | total allowable catch | yes |
| `sections/SectionsIntegrated.shp` | DFO herring section boundaries | yes |
| `Biosample_Strait_of_Georgia.csv` | DFO herring biosample database, Strait of Georgia, 1946-2024 | no, available from DFO on request |
| `SOG_herring_catch_all_fishing_types*.csv` | DFO commercial catch by season, area, fishery and gear | no, available from DFO on request |

Every series derived from the two restricted extracts is included in `Outputs/derived/`.

The biosample extract has no usable header row; column names are assigned from `BIO_COLS` in
`00_setup.R`, and two-digit years are converted to calendar years.

## Requirements

R 4.5 or later, with tidyverse, mgcv, nlme, lme4, segmented, strucchange, sf, adehabitatHR,
spdep, geosphere, ineq, rnaturalearth (and rnaturalearthhires), gratia, patchwork, ggrepel,
showtext, scales, knitr and rmarkdown.

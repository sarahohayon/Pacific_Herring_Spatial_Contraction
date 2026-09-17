# Biomass recovery masks collective memory loss and spatial collapse in Pacific herring

Analysis code for the Strait of Georgia Pacific herring study. Every statistic, table and
figure in the manuscript, main text and supplement alike, is produced by the files here,
starting from the raw DFO extracts. Nothing reads a cached result, so a figure and the
statistic it illustrates cannot drift apart.

## Layout

```
Outputs/code/      the pipeline: 00_setup.R, 01-09, and the figure scripts
Outputs/derived/   derived series written by 01 and read by everything else
Outputs/tables/    one CSV per supplementary table, plus supporting tables
Outputs/figures/   every generated figure, PDF and PNG
Figures/           final figures as used in the manuscript, assembled from the panels
Raw_data/          DFO and NOAA inputs (see Data below)
run_all.R          renders the whole pipeline in order
```

## How to run

```bash
cd Outputs/code && Rscript ../../run_all.R
```

`01` must run first: it writes the derived series that everything else reads. It is also the
only file that needs the restricted DFO extracts, so without them start at `02`, which runs
from the published `Outputs/derived/` series.

Run time is dominated by `01` (a 100 MB biosample file) and `05` (76 GAMM fits); budget about
15 minutes. If `rmarkdown::render` cannot find pandoc outside RStudio:

```bash
export RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64"
```

## What each file produces

| File | Produces |
|---|---|
| `00_setup.R` | paths, analysis constants, presence rule, shared theme and palettes. Sourced by every file; the only place to change the constants below. |
| `01_data_preparation.Rmd` | all derived series into `Outputs/derived/` |
| `02_contraction_metrics.Rmd` | eleven contraction metrics and the PCA behind the choice of response: **Table S6** |
| `03_biomass_occupancy_decoupling.Rmd` | biomass-occupancy decoupling, phenology: **Fig. 1a-f**, **Fig. S1**, **Tables S1, S7** |
| `04_demography_selectivity.Rmd` | fishery selectivity and age structure: **Fig. 2a-d**, **Fig. S3** |
| `05_occupancy_drivers_GAMM.Rmd` | driver models, robustness and subset refits: **Fig. 3a-c**, **Figs. S4-S5**, **Tables S2-S5, S7** |
| `06_spatial_asymmetry_refuge.Rmd` | spatial collapse and the refuge: **Fig. 4a-d** |
| `07_supplementary_catch_figures.Rmd` | **Figs. S6-S7** |
| `08_roe_gillnet_age_composition.Rmd` | **Fig. S2** and the age-5 threshold diagnostics |
| `09_prop_old_by_area.R` | per-area age structure behind the refuge result (Area 14 odds ratios) |
| `Fig1_bcd_SoG_pipeline.R` | Figure 1 b, c, d in manuscript style, fed by the pipeline series |
| `Fig2bcd_manuscript_style.R` | Figure 2 b, c, d |
| `Fig2c_length_at_age_tests.R` | length-at-age mixed models quoted in the Results |
| `Fig3b_null_draws.R` | the null draws behind Fig. 3b, **Table S8** and the randomisation P values in the text |
| `Fig3b_sst_increment_test.R` | what SST adds beyond selectivity (the increment test in Table S8) |
| `Fig3_randomisation_test.R`, `Fig3_randomisation_R2_by_driver.R` | slope and variance randomisations (**Table S8**) |
| `Fig3_manuscript_style.R` | Figure 3 panels a and b, drawn from those null draws |

Main-text figures are written twice: one file per panel (`Fig1a_SSB`, `Fig4c_site_loss_vs_exploitation`,
and so on) for assembling the final layout, and one pre-assembled composite
(`Figure1_biomass_occupancy_decoupling`, `Figure4_four_panels_untagged`) for checking that the
figure reads as a whole. The rendered HTML report beside each Rmd carries its printed output.

## The decisions that propagate everywhere

They live at the top of `00_setup.R` and are stated in the Methods:

* **Presence rule** (`keep_surveyed`): a spawn record counts as presence if it was surveyed,
  including records DFO leaves without a biomass estimate; records DFO classes as *Incomplete*
  are excluded, because that category first appears in the 1990s and would add presence only
  after the collapse. Of 5,751 Strait of Georgia records, 257 carry no biomass estimate: 115 are
  surveyed and kept, 142 are Incomplete and dropped. Table S7 shows the results under all three
  possible rules.
* `MIN_ACTIVE = 6` - a section is historically used if spawn was recorded there in at least six
  distinct years, giving the occupancy denominator of **19 sections**.
* `MIN_FISH = 30` - a biosample is used only if it holds at least 30 aged fish.
* `OLD_AGE = 5` - old fish are age 5 and older, separating repeat spawners from first-time
  spawners at age 3.
* `BREAK_YEAR = 1984` - the collapse boundary, set at the largest single-year decline in
  occupancy; change-point analysis places the structural break after 1978 (95% CI 1977-1980).
* **Map rule** (in `06`): the two-period maps show sites in regular use, meaning at least 2 spawn
  years in any 5-year moving window within the period.

## Two data objects that look alike and are not

`01` writes both `site_year_spawn.csv` and `site_year_spawn_all.csv`. The first requires a
location to have coordinates and feeds anything that places spawning in space (contraction
metrics, maps, per-area change). The second keeps every record and feeds the site-level
occupancy histories, so a site is never dropped from its own history for lacking a position.
Ten records in the Strait have no coordinates; using the wrong object shifts the site-loss
counts by one site.

## Data

| File | Source | In this repository |
|---|---|---|
| `Pacific_herring_spawn_index_data_2025_EN.csv` | DFO Pacific Herring Spawn Index, 1951-2025 | yes |
| `Lighthouse_Stations_SST_Combined.csv` | DFO BC lighthouse SST programme | yes |
| `ONI_index_1950_2025.csv`, `PDO_1854_2025.csv` | NOAA climate indices | yes |
| `Pacific herring catch and landed value dfo_INFLATION.xlsx` | DFO landed value with CPI deflator | yes |
| `SOG_TAC.csv` | total allowable catch | yes |
| `sections/SectionsIntegrated.shp` | DFO herring section boundaries | yes |
| `Biosample_Strait_of_Georgia.csv` | DFO herring biosample database, SoG extract, 1946-2024 | no, available from DFO on request |
| `SOG_herring_catch_all_fishing_types*.csv` | DFO commercial catch by season, area, fishery, gear | no, available from DFO on request |

The two restricted extracts are the only inputs not included. Every series derived from them is
published in `Outputs/derived/`, so `02` onward reproduce without them.

The biosample file ships from DFO without a usable header row; column names are assigned
positionally from `BIO_COLS` in `00_setup.R`, and the two-digit sampling year is recoded to a
calendar year with a century break at 29/30.

## Requirements

R 4.5 or later. Packages: tidyverse, mgcv, nlme, lme4, segmented, strucchange, sf, adehabitatHR,
spdep, geosphere, ineq, rnaturalearth (with rnaturalearthhires), gratia, patchwork, ggrepel,
showtext, scales, knitr, rmarkdown.

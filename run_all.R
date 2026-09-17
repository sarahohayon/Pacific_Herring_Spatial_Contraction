## ---------------------------------------------------------------------------
## Reproduce every statistic and figure in the manuscript, in the order the
## paper presents them.
##
##   cd Outputs/code && Rscript ../../run_all.R
##
## or open Herring_analysis.Rproj in RStudio and source this file.
##
## 01 must run first: it writes the derived series that everything else reads.
## 01 is also the only file that needs the two restricted DFO extracts (see
## README). Without them, start at 02: the published Outputs/derived/ series
## carry every input the later steps use.
## ---------------------------------------------------------------------------

if (!file.exists("00_setup.R"))
  stop("Run this from Outputs/code/ (the folder holding 00_setup.R).")

reports <- sort(list.files(pattern = "^0[1-8].*\\.Rmd$"))

scripts <- c(
  "09_prop_old_by_area.R",           # per-area age structure behind the refuge result
  "Fig1_bcd_SoG_pipeline.R",         # Figure 1 b, c, d
  "Fig2bcd_manuscript_style.R",      # Figure 2 b, c, d
  "Fig2c_length_at_age_tests.R",     # length-at-age mixed models (Results)
  "Fig3b_null_draws.R",              # null draws: Fig 3b, Table S8 and the text P values
  "Fig3b_sst_increment_test.R",      # what SST adds beyond selectivity
  "Fig3_randomisation_test.R",       # slope randomisations (Table S8)
  "Fig3_randomisation_R2_by_driver.R",
  "Fig3_manuscript_style.R"          # Figure 3 panels, drawn from the null draws above
)

for (f in reports) { message("== render ", f); rmarkdown::render(f, quiet = TRUE) }
for (f in scripts) { message("== run    ", f); source(f, echo = FALSE) }

message("\nDone. Derived series in ../derived, tables in ../tables, figures in ../figures.")

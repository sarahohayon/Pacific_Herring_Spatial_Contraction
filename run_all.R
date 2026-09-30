## ---------------------------------------------------------------------------
## Reproduce every statistic, table and figure in the paper, in order.
##
##   cd Outputs/code && Rscript ../../run_all.R
##
## 01 runs first: it builds the derived series that everything else reads, and it
## needs the DFO catch and TAC data, which are available on request and not in this
## repository (see README); so do the catch panels of Figures 1 and 2.
##
## Each step runs in its own R session. Scripts in extra/ are not part of the paper.
## ---------------------------------------------------------------------------

if (!file.exists("00_setup.R"))
  stop("Run this from Outputs/code/ (the folder that holds 00_setup.R).")

reports <- sort(list.files(pattern = "^0[1-8].*\\.Rmd$"))   # 01-06 and 08

scripts <- c(
  "09_prop_old_by_area.R",          # Area 14 versus the peripheral areas (Results)
  "Fig1_bcd_SoG_pipeline.R",        # Figure 1 b, c, d
  "Fig2bcd_manuscript_style.R",     # Figure 2 b, c, d
  "Fig2c_length_at_age_tests.R",    # length-at-age tests (Results)
  "Fig3b_null_draws.R",             # randomization nulls for Figure 3b and Table S8b
  "Fig3b_null_draws_extra.R",       # ... the remaining drivers (run after the line above)
  "Fig3b_sst_increment_test.R",     # does SST add to older fish? (Table S8c)
  "Fig3_randomisation_test.R",      # slope randomizations and gear-era model (Table S8a)
  "Fig3_manuscript_style.R"         # Figure 3 a, b
)

for (f in reports) {
  message("== render ", f)
  rmarkdown::render(f, quiet = TRUE, envir = new.env())
}
for (f in scripts) {
  message("== run    ", f)
  status <- system2("Rscript", f)
  if (status != 0) stop(f, " failed")
}

message("\nDone. Derived series in Outputs/derived, tables in Outputs/tables, figures in Outputs/figures.")

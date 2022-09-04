# libraries ----
library(here)

# processing scripts ----
source(here("scripts", "processing_scripts", "01-summarize_trapnest_data.R"))
source(here("scripts", "processing_scripts", "02-calculate_land_cover_metrics.R"))
source(here("scripts", "processing_scripts", "03-clean_trait_data.R"))

# analysis scripts -----
source(here("scripts", "analysis_scripts", "01-crit1_func_alpha.R"))
source(here("scripts", "analysis_scripts", "02-crit2_func_alpha.R"))
source(here("scripts", "analysis_scripts", "03-crit3_RLQ.R"))
source(here("scripts", "analysis_scripts", "04_create_figure_S1_and_S2.R"))
################################################################################
# Accompanying code for the paper: 
#   No environmental filtering of wasps and bees during urbanization
#
# Paper authors: 
#   Garland Xie (1), 
#   Nicholas Sookhan (1), 
#   Kelly Carscadden (2), 
#   and Scott MacIvor (1)
#
# Corresponding author: 
#   Garland Xie (1)
#
# Affiliations: 
#   (1) Department of Biological Sciences, 
#       University of Toronto Scarborough,
#       1265 Military Trail, Toronto, ON, M1C 1A4, Canada
#       email: garland.xie@mail.utoronto.ca, 
#              nicholas.sookhan@mail.utoronto.ca
#              scott.macivor@mail.utoronto.ca
#   (2) Department of Ecology and Evolutionary Biology,
#       University of Colorado Boulder,
#
# Purpose of this R script: to reproduce Figure S3 (Appendix)

# libraries --------------------------------------------------------------------
library(here)      # for creating relative file-paths
library(tidyverse) 
library(patchwork) # for creating multi-panel figures
library(readxl)    # for reading excel files
library(janitor)   # for cleaning column names

# import -----------------------------------------------------------------------

# abundance per site 
int_raw <- read_excel(
  here("data/input_data", "trap_nest_jsm_Aug10_2021.xlsx")
)

# site data
site <- read_excel(
  here("data/input_data", "site_jsm_edits_Aug10_2021.xlsx")
)

# data cleaning: sites ---------------------------------------------------------

# some prep
int_tidy <- int_raw %>%
  janitor::clean_names() %>%
  filter(!is.na(id)) %>%
  
  # keep only species-level analyses
  filter(lower_species != "Hylaeus_sp")

# get sites that outside the TO boundary
outside_TO <- c(
  "GAJVv", 
  "SQWq3",
  "N53op", 
  "auCMf", 
  "lWpWV",
  "Z42dv",
  "h6kO1",
  "sC5O0", 
  "wB2e4"
)

# get sites that were sampled across 2011-2013 within TO
all_years <- site %>%
  
  # keep sites sampled across all three years (2011-2013)
  filter(Year_2011 == "Y" & Year_2012 == "Y" & Year_2013 == "Y") %>%
  
  # keep sites within TO
  filter(!(ID %in% outside_TO)) %>%
  
  pull(ID)


# Figure S3: histograms for abundance -----------------------------------------------

(hist_bees_ab <- int_tidy %>%
   filter(taxa_ls == "Bee", id %in% all_years) %>%
   group_by(id, year) %>%
   summarize(sum_broods = sum(no_broodcells)) %>%
   ggplot(aes(x = sum_broods)) + 
   geom_histogram(binwidth = 1) + 
   ylim(0, 10) + 
   facet_wrap(~year) + 
   labs(
     title = "A)",
     x = "Number of completed brood cells (bees)",
     y = NULL
   ) + 
   theme_bw() 
)

(hist_wasp_ab <- int_tidy %>%
    filter(taxa_ls == "Wasp", id %in% all_years) %>%
    group_by(id, year) %>%
    summarize(sum_broods = sum(no_broodcells)) %>%
    ggplot(aes(x = sum_broods)) + 
    geom_histogram(binwidth = 1) + 
    facet_wrap(~year) + 
    ylim(0, 10) + 
    labs(
      title = "B)",
      x = "Number of completed brood cells (wasps)",
      y = NULL) + 
    theme_bw() 
)

(hist_total_ab <- int_tidy %>%
    filter(id %in% all_years) %>%
    group_by(id, year) %>%
    summarize(sum_broods = sum(no_broodcells)) %>%
    ggplot(aes(x = sum_broods)) + 
    geom_histogram(binwidth = 1) + 
    facet_wrap(~year) + 
    ylim(0, 10) + 
    labs(
      title = "C)",
      x = "Number of completed brood cells",
      y = NULL) + 
    theme_bw() 
)

# multi-panel histograms
(hist_ab <- hist_bees_ab / hist_wasp_ab / hist_total_ab)

# save to disk -----------------------------------------------------------------

ggsave(
  plot = hist_ab, 
  filename = here(
    "output/data_appendix_output", 
    "Xie_et_al-2021-FigureS3-JAE.png"
  ),
  device = "png", 
  height = 5, 
  width  = 5
)

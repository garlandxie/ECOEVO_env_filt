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
# Purpose of this R script: to clean data for the community matrix

# libraries --------------------------------------------------------------------
library(here)      # for creating relative file-paths
library(dplyr)     # for manipulating data
library(ggplot2)   # for visualizing data
library(tidyr)     # for making wide tables
library(patchwork) # for creating multi-panel figures
library(janitor)   # for cleaning column names

# import -----------------------------------------------------------------------

# abundance per site 
int_raw <- read.csv(
  here("data/input_data", "trap_nest.csv")
)

# site data
site <- read.csv(
    here("data/input_data", "site_data.csv")
)

# data cleaning: sites ---------------------------------------------------------

# some prep
int_tidy <- int_raw %>%
  janitor::clean_names() %>%
  filter(!is.na(id)) %>%
  
  # keep only species-level analyses
  filter(species != "Hylaeus_sp")

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

# data cleaning: bees ----------------------------------------------------------

# compare number of brood cells (including parasitized) and number of alive cells
int_raw %>%
  ggplot() + 
  geom_histogram(aes(x = No_broodcells), fill = 'red',alpha = 0.5) + 
  geom_histogram(aes(x = No_alive), fill = 'blue', alpha = 0.5) +
  labs(x = NULL) + 
  facet_wrap(~Year) + 
  theme_bw()

# community data matrix
# proxy of abundance: number of brood cells
# aggregrated across all years (2011-2013)
broods <- int_tidy %>%
  filter(id %in% all_years) %>%
  group_by(id, species) %>%
  summarize(total_alive = sum(no_broodcells)) %>%
  ungroup() %>%
  pivot_wider(names_from = species, values_from = total_alive) %>%
  mutate(across(everything(), ~replace_na(., 0))) %>%
  column_to_rownames(var = "id")

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

# Figure S4: histograms for species richness -----------------------------------

(hist_sr_bees <- int_tidy %>%
   filter(taxa_ls == "Bee", id %in% all_years) %>%
   group_by(year, id) %>%
   summarize(species_richness = length(unique(species))) %>%
   ggplot(aes(x = species_richness)) + 
   geom_histogram(bins = 5, binwidth = 1) + 
   ylim(0, 75) + 
   facet_wrap(~year) + 
   scale_x_continuous(
     breaks = c(2,4,6,8)
   ) +
   labs(
     title = "A)",
     x = "Species richness (bees)",
     y = NULL
   ) + 
   theme_bw()
)

(hist_sr_wasps <- int_tidy %>%
    filter(taxa_ls == "Wasp", id %in% all_years) %>%
    group_by(id, year) %>%
    summarize(species_richness = length(unique(species))) %>%
    ggplot(aes(x = species_richness)) + 
    geom_histogram(bins = 5, binwidth = 1) + 
    ylim(0, 75) + 
    facet_wrap(~year) + 
    scale_x_continuous(
      breaks = c(2,4,6,8)
    ) + 
    labs(
      title = "B)",
      x = "Species richness (wasps)",
      y = NULL
    ) + 
    theme_bw()
)

(hist_sr_total <- int_tidy %>%
    filter(id %in% all_years) %>%
    group_by(id, year) %>%
    summarize(species_richness = length(unique(species))) %>%
    ggplot(aes(x = species_richness)) + 
    geom_histogram(bins =5, binwidth = 1) + 
    ylim(0, 75) + 
    facet_wrap(~year) +
    scale_x_continuous(
      breaks = c(2,4,6,8)
    ) + 
    labs(
      title = "C)",
      x = "Species richness (bees and wasps)",
      y = NULL
    ) + 
    theme_bw()
)

# multi-panel figure
(hist_SR <- hist_sr_bees / hist_sr_wasps / hist_sr_total)

# save to disk -----------------------------------------------------------------

ggsave(
  plot = hist_SR, 
  filename = here(
    "output", "data_appendix_output", 
    "Xie_et_al-2021-FigureS4-JAE.png"
  ),
  device = "png", 
  height = 5, 
  width  = 5
)

ggsave(
  plot = hist_ab, 
  filename = here(
    "output", "data_appendix_output", 
    "Xie_et_al-2021-FigureS3-JAE.png"
  ),
  device = "png", 
  height = 5, 
  width  = 5
)

write.csv(
  broods, 
  file = here(
    "data", "analysis_data", 
    "comm_matrix_B.csv"
    )
  )

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
library(tidyverse) 
library(patchwork) # for creating multi-panel figures

# import -----------------------------------------------------------------------

# abundance per site 
int_raw <- readxl::read_excel(
  here(
    "data/original", 
    "JSM_Data_TrapnestSynthesis2020_FINAL.xlsx"
    ), 
  sheet = "B-interactiondata(per-nest)"
  )

# site data
site <- readxl::read_excel(
  here(
    "data/original", 
    "JSM_Data_TrapnestSynthesis2020_FINAL.xlsx"
  ), 
  sheet = "A-metadata(per-site)"
)

# data cleaning: sites ---------------------------------------------------------

# some prep
int_tidy <- int_raw %>%
  janitor::clean_names() %>%
  filter(!is.na(site_id)) %>%
  mutate(
    
    lower_species = case_when(
      lower_species == "Hylaeus_sp" ~ "Hylaeus_punctatus", 
      TRUE ~ lower_species)
  )

# get sites that outside the TO boundary
outside_TO <- c(
  "Dumesh", 
  "Kavanah",
  "Lynott", 
  "RangersGround", 
  "RangersRoof",
  "Chute"
)

# get sites that were sampled across 2011-2013 within TO

all_years <- site %>%
  filter(!Site_ID %in% outside_TO) %>%
  filter(Year_2011  == "Y" &
         Year_2012 ==  "Y" & 
         Year_2013 ==  "Y") %>%
  pull(Site_ID)

# plot: histograms -------------------------------------------------------------

hist_bees <- int_tidy %>%
  filter(taxa_ls == "Bee") %>%
  ggplot(aes(x = no_broodcells)) + 
  geom_histogram(binwidth = 1) + 
  ylim(0, 200) + 
  facet_wrap(~year) + 
  labs(
    title = "A)",
    x = "Number of occupied brood cells (for bees)",
    y = "Count"
  ) + 
  theme_bw()

hist_wasp <- int_tidy %>%
  filter(taxa_ls == "Wasp") %>%
  ggplot(aes(x = no_broodcells)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(~year) + 
  ylim(0, 200) + 
  labs(
    title = "B)",
    x = "Number of occupied brood cells (for wasps)",
    y = NULL) + 
  theme_bw()

# multi-panel histograms
hist_all <- hist_bees + hist_wasp

# data cleaning: bees ----------------------------------------------------------

# community data matrix
# proxy of abundance: number of brood cells
# aggregrated across all years (2011-2013)
broods <- int_tidy %>%
  group_by(site_id, lower_species) %>%
  summarize(total_alive = sum(no_broodcells)) %>%
  pivot_wider(names_from = lower_species, values_from = total_alive) %>%
  ungroup() %>%
  mutate(across(everything(), ~replace_na(., 0)))
   
# get sites that were sampled in all three years
b3 <- broods$site_id[broods$site_id %in% all_years]

# subset accordingly
broods_tidy <- broods %>% 
  filter(site_id %in% b3) %>%
  filter(!(site_id %in% outside_TO)) %>%
  column_to_rownames(var = "site_id")

# save to disk -----------------------------------------------------------------

write.csv(
  broods_tidy, 
  file = here(
    "data/final", 
    "comm_matrix_B.csv"
    )
  )

ggsave(
  plot = hist_all, 
  filename = here(
    "output/figures/supp", 
    "fig-supp-histogram.png"
  ),
  device = "png", 
  height = 3, 
  width  = 7
)





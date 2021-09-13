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
# Purpose of this R script: to clean data for the trait matrix

# libraries --------------------------------------------------------------------
library(here)     # for creating relative file-paths
library(dplyr)    # for manipulating data 
library(readxl)   # for reading excel files

# import -----------------------------------------------------------------------

traits <- read_excel(
  here(
    "data/original", 
    "bee_wasp_traits_aug10.xlsx")
  )

# data cleaning: trait matrix with voltinism -----------------------------------

traits_tidy <- traits %>%
  janitor::clean_names() %>%
  
  # insert missing values
  mutate(body_size = na_if(body_size, "NA")) %>% 
  
  # change Campunula to Campunulaecae 
  mutate(specialization = case_when(
    specialization == "Genus (Campanula)" ~ "Family (Campanulaceae)", 
    TRUE ~ specialization)
    ) %>%

  # coerce into factor variables
  mutate(
    
    origin = factor(origin),
    
    # uses most common material
    nesting_material = factor(nesting_material), 
    primary_diet = factor(primary_diet),
    voltinism = factor(voltinism), 
    specialization = factor(specialization),
    trophic_rank = factor(trophic_rank)
    
    ) %>%
  
  # change to numeric values
  mutate(body_size = as.numeric(body_size)) %>%
  
  # change species names so it is consistent 
  # with the community matrix dataset
  mutate(
    species = stringr::str_replace(
      species, 
      pattern = " ", 
      replacement = "_")
    ) %>%
  
  select(-bee_wasp, -authority, -family)

# data cleaning: trait matrix without voltinism -------------------------------

# remove volitinism (as requested by JSM)
# remember to re-run the data pipeline with this dataset!
traits_no_volt <- traits_tidy %>%
  select(-voltinism)

# save to disk -----------------------------------------------------------------

saveRDS(
  traits_tidy, 
  file = here(
    "data/final", 
    "traits_with_volt.rds")
  )

saveRDS(
  traits_no_volt, 
  file = here(
    "data/final",
    "traits_no_volt.RDS"
  )
)




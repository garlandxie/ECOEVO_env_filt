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
library(janitor)  # for cleaning column names

# import -----------------------------------------------------------------------

traits <- read.csv(
  here("data/input_data", "traits.csv"),
  stringsAsFactors = FALSE)

# data cleaning: edit existing traits ------------------------------------------

traits_tidy <- traits %>%
  
  janitor::clean_names() %>%

  # insert missing values
  mutate(its = na_if(its, "-")) %>% 
  
  mutate(nest_mat_type = case_when(
    nest_mat_type == "Resin + Mud" ~ "Resin",
    nest_mat_type == "Leaf pulp + Sand" ~ "Leaf pulp",
    nest_mat_type == "Leaf pulp + Pebbles" ~ "Leaf pulp",
    nest_mat_type == "Leaf cut + Pulp + Mud" ~ "Leaf pulp",
    nest_mat_type == "Mud + Leaf pulp" ~ "Mud", 
    TRUE ~ nest_mat_type)
    ) %>%
  
  mutate(specialization = case_when(
    
    # family
    specialization == "Family (Asteraceae)"          ~ "Family", 
    specialization == "Family (Chrysomelidae)"       ~ "Family",
    specialization == "Family (Aphididae)"           ~ "Family",
    specialization == "Family (Campanulaceae)"       ~ "Family",
    specialization == "Family (Aphididae)"           ~ "Family", 
    
    # order 
    
    specialization == "Order (Lepidoptera)" ~ "Order",
    specialization == "Order (Araneae)" ~ "Order",
    specialization == "Multi-Order (Coleoptera, Lepidoptera)" ~ "Multi-Order",
    
    TRUE ~ specialization)
    ) %>%
   
  # coerce into factor variables
  mutate(
    
    origin              = factor(origin),
    nest_mat_type       = factor(nest_mat_type), 
    num_nest_mat_types  = factor(num_nest_mat_types),
    diet                = factor(diet),
    specialization      = factor(specialization),
    rank                = factor(rank)
    
    ) %>%
  
  
  # change to numeric values
  mutate(its = as.numeric(its)) %>%
  
  # change species names so it is consistent 
  # with the community matrix dataset
  mutate(
    species = stringr::str_replace(
      species, 
      pattern = " ", 
      replacement = "_")
  ) %>%
  select(-bee_wasp, -authority, -family)

# check packaging --------------------------------------------------------------

glimpse(traits_tidy)

# save to disk -----------------------------------------------------------------

write.csv(
  x = traits_tidy,
  file = here("data", "analysis_data", "traits_tidy.csv"), 
  row.names = FALSE
)






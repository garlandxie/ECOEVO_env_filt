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
    "bee_wasp_traits.xlsx")
  )

# data cleaning: ---------------------------------------------------------------

traits_tidy <- traits %>%
  clean_names() %>%
  
  # insert missing values
  mutate(
    
    itd = na_if(itd, "na"),
    emer_per = na_if(emer_per, "na")
    
    ) %>% 

  # coerce into factor variables
  mutate(
    
    native_y_n = factor(native_y_n),
    nest       = factor(nest),
    diet       = factor(diet),
    volt       = factor(volt)
    
    ) %>%
  
  # change to numeric values
  mutate(
    
    itd      = as.numeric(itd),
    emer_per = as.numeric(emer_per) 
    
    ) %>%

  select(-1, -family)

# double-check
str(traits_tidy)

# table prep -------------------------------------------------------------------

table_S1 <- 
  traits %>%
  mutate(genus   = strsplit(spp, split = "_"),
         genus   = sapply(genus, "[", 1),
         species = strsplit(spp, split = "_"),
         species = sapply(species, "[", 2)
         ) %>%
  
  mutate(itd = na_if(itd, "na"),
         itd = as.numeric(itd),
         itd = round(itd, digits = 2)
         ) %>%
  
  select("Family"           = family,
         "Genus"            = genus, 
         "Species"          = species, 
         "Native status"    = native_y_n,
         "Nesting material" = nest, 
         "Diet"             = diet, 
         "Voltinism"        = volt, 
         "Body size (ITD)"  = itd
  ) 

# save to disk -----------------------------------------------------------------

saveRDS(
  traits_tidy, 
  file = here(
    "data/final", 
    "traits.rds")
  )


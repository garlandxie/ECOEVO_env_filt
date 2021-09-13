################################################################################
# Accompanying code for the paper: 
#   No environmental filtering of wasps and bees during urbanization
#
# Paper authors: 
#   Garland Xie          (1), 
#   Nicholas Sookhan     (1), 
#   Kelly Carscadden     (2), 
#   and J. Scott MacIvor (1)
#
# Corresponding authors for this script:  
#   Garland Xie      (1)
#   Nicholas Sookhan (2)
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
# Purpose of this R script: to encrypt confidential site data

# libraries --------------------------------------------------------------------
library(readxl)    # for reading excel files  
library(stringi)   # for creating random strings of characters
library(here)      # for creating relative file paths 
library(dplyr)     # for manipulating our data 

# import -----------------------------------------------------------------------
site <- readxl::read_excel(
  here::here(
    "data/original", 
    "JSM_Data_TrapnestSynthesis2020_FINAL.xlsx"
  ), 
  sheet = "A-metadata(per-site)"
)

int_raw <- readxl::read_excel(
  here(
    "data/original", 
    "JSM_Data_TrapnestSynthesis2020_FINAL.xlsx"
  ), 
  sheet = "B-interactiondata(per-nest)"
)

# generate random characters ---------------------------------------------------
set.seed(1)
random_ID <- stringi::stri_rand_strings(
  n = length(site$Site_ID), 
  length = 5, 
  pattern = "[A-Za-z0-9]"
  )

# clean data -------------------------------------------------------------------

# site
site_tidy <- site
site_tidy$ID <- random_ID

site2 <- site_tidy[,c("ID", "Site_ID")]

# trap nests
int_tidy <- int_raw %>%
  
  left_join(
    site_tidy, 
    by = c("SiteID" = "Site_ID")
    ) %>%
  
  select(ID, colnames(int_raw)[2:17])

# clean more for sites
site_tidy <- site_tidy %>%
  select(2:9) %>%
  relocate(ID)

# save to disk -----------------------------------------------------------------

write.csv(
  x = int_tidy, 
  file = here("data/original", "trap_nest.csv")
)

write.csv(
  x = site_tidy, 
  file = here("data/original", "site.csv")
)

write.csv(
  x = site2, 
  file = here("data/original", "site_IDs.csv")
)
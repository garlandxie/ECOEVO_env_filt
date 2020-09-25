library(readxl)
library(stringi)
library(here)
library(dplyr)

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

# trap nests
int_tidy <- int_raw %>%
  
  left_join(
    site_tidy, 
    by = c("SiteID" = "Site_ID")
    ) %>%
  
  select(ID, colnames(int_raw)[2:17])

# save to disk -----------------------------------------------------------------

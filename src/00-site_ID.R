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
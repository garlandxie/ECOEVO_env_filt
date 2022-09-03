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
# Purpose of this R script: to conduct statistical analyses for CT criteria II 

# libraries --------------------------------------------------------------------
library(here)        # for creating relative file-paths
library(dplyr)       # for manipulating data
library(ggplot2)     # for visualizing data
library(GGally)      # for visualizing pairwise distance matrices
library(car)         # for calculating variance inflation factors
library(sp)          # for geospatial analyses 
library(spdep)       # for analyzing spatial autocorrelation  
library(patchwork)   # for making multi-panel plots
library(gghighlight) # for highlighting ggplot features
library(tibble)      # for creating tibbles 

# import -----------------------------------------------------------------------

## |- site ----
site <- read.csv(
  here(
    "data", "input_data",
    "site_data.csv")
)

## |- ses mfd ----
ses_mfd <- readRDS(
  here("data", "intermediate_data",
       "ses_mfd.rds")
  )

## |- land cover -----
l_250 <- read.csv(
  here("data", "intermediate_data", 
       "land_use_250.csv")
  )

l_500 <- read.csv(
  here("data", "intermediate_data", 
       "land_use_500.csv")
  )

# data cleaning ----------------------------------------------------------------

## |- 250m scale ----
reg_250 <-  ses_mfd %>%
  full_join(l_250, by = c("site_id" = "site")) %>%
  full_join(site, by = c("site_id" = "ID")) %>%
  select(site  = site_id, 
         sr    = ntaxa, 
         Habitat_type,
         longs = Longitude, 
         lats  = Latitude, 
         ntaxa,
         mfd.obs,
         ses_mfd, 
         p_value,
         perc_tree_250, 
         perc_grass_250, 
         perc_urb_250) %>%
  
  filter(
    
    # remove sites that have only 1 species
    !is.na(ses_mfd),
    
    # remove sites that are outside TO boundary
    !is.na(perc_tree_250)  &
    !is.na(perc_grass_250) & 
    !is.na(perc_urb_250)  
  ) %>%
  
  mutate(Habitat_type = factor(Habitat_type))

## |- 500m scale ----
reg_500 <- ses_mfd %>%
  full_join(l_500, by = c("site_id" = "site")) %>%
  full_join(site, by = c("site_id" = "ID")) %>%
  select(site  = site_id, 
         sr    = ntaxa, 
         Habitat_type, 
         longs = Longitude, 
         lats  = Latitude, 
         ntaxa,
         mfd.obs,
         ses_mfd, 
         p_value,
         perc_tree_500, 
         perc_grass_500, 
         perc_urb_500) %>%
  
  filter(
  
    # remove sites with only 1 species
    !is.na(ses_mfd),
    
    # remove sites outside TO boundary
    !is.na(perc_tree_500)  &
    !is.na(perc_grass_500) & 
    !is.na(perc_urb_500)  
    ) %>%
  
    mutate(Habitat_type = factor(Habitat_type))

# Figure S5: correlation matrix ------------------------------------------------

# Pearson's correlation matrix for 250m
pairs_250 <- reg_250 %>%
  select(
    "% Closed Green"  = perc_tree_250, 
    "% Open Green" = perc_grass_250,
    "% Impervious" = perc_urb_250) %>%
  ggpairs() + 
  labs(title = "250m spatial scale")

# Figure S6: correlation matrix ------------------------------------------------

# Pearson's correlation matrix for 500m
pairs_500 <- reg_500 %>%
  select(
    "% Closed Green"  = perc_tree_500, 
    "% Open Green" = perc_grass_500,
    "% Impervious" = perc_urb_500) %>%
  ggpairs() + 
  labs(title = "500m spatial scale")

# hypothesis testing ----

## |- multiple regression (250m) ----

# first fit
lm_250_v1 <- lm(ses_mfd ~ 
                  perc_grass_250 + perc_tree_250 + perc_urb_250 + Habitat_type, 
             data = reg_250)
vif(lm_250_v1)

# remove tree cover
lm_250_v2 <- update(lm_250_v1, ~. -perc_tree_250)
vif(lm_250_v2)

# get summary
sum_250 <- summary(lm_250_v2)

# show model diagnostics in a non-interactive manner 
plot(lm_250_v2, which = c(1))
plot(lm_250_v2, which = c(2))
plot(lm_250_v2, which = c(3))
plot(lm_250_v2, which = c(5))

## |- multiple regression (500m) ---- 

# first fit
lm_500_v1 <- lm(ses_mfd ~ 
                  perc_grass_500 + perc_tree_500 + perc_urb_500 + Habitat_type, 
             data = reg_500)
vif(lm_500_v1)

# remove tree cover
lm_500_v2 <- update(lm_500_v1, ~. -perc_tree_500)
vif(lm_500_v2)

# get summary
sum_500 <- summary(lm_500_v2)

# show model diagnostics in a non-interactive manner 
plot(lm_500_v2, which = c(1))
plot(lm_500_v2, which = c(2))
plot(lm_500_v2, which = c(3))
plot(lm_500_v2, which = c(5))

# spatial autocorrelaton tests ----

## |- 250m----
coords_250 <- reg_250 %>%
  select(site, longs, lats) 

coordinates(coords_250) <- ~lats + longs
proj4string(coords_250) <- CRS("+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 
                            +units=m +no_defs")

# formal spatial autocorrelation on the regression residuals
# for Global Moran's I
lm.morantest(lm_250_v2, 
             listw = nb2listw(
               knn2nb(
                 knearneigh(x = coords_250,
                            k = 8)
                 ),
               style = "W"
               )
             )

## |- 500m ----
coords_500 <- reg_500 %>%
  select(site, longs, lats) 

coordinates(coords_500) <- ~lats + longs
proj4string(coords_500) <- CRS("+proj=utm +zone=17 +ellps=GRS80 +datum=NAD83 
                            +units=m +no_defs")

lm.morantest(lm_500_v2, 
             listw = nb2listw(
               knn2nb(
                 knearneigh(x = coords_500, 
                            k = 8)
                 ),
               style = "W")
             )

# prep for Figure 4 and Figure S7 ----------------------------------------------

# get R-squared values for labelling

rq_250 <- summary(lm_250_v2)$adj.r.squared
rq_500 <- summary(lm_500_v2)$adj.r.squared

rq_lab_250 <- bquote("Adj-R"^2: .(format(rq_250, digits = 3)))
rq_lab_500 <- bquote("Adj-R"^2: .(format(rq_500, digits = 2)))

# Figure 4: ses.MFD vs % impervious surfaces (250m scale) ----------------------

# scatterplots: ses.MFD vs % impervious surfaces
(part_urb_250 <- reg_250 %>%
  select(site, 
         ses_mfd, 
         p_value, 
         perc_urb_250) %>%
  mutate(pred_urb_250 = predict(lm_250_v2, terms = "perc_urb_250"),
         part_urb_250 = pred_urb_250 + resid(lm_250_v2)) %>%
  ggplot(aes(x = perc_urb_250, y = part_urb_250)) + 
  geom_point() + 
  gghighlight(p_value < 0.05, use_direct_label = FALSE) + 
  geom_hline(yintercept = 0) + 
  annotate(
    geom = "text",
    x = 75, 
    y = 2, 
    label = rq_lab_250
    ) + 
  labs(y = "ses.MFD (partial residuals)",
       x = "% Impervious surface (250m spatial scale)") + 
  theme_bw()
  )

# Figure S7: ses.MFD vs % impervious surfaces (500m scale) ---------------------

(part_urb_500 <- reg_500 %>%
  select(
         ses_mfd, 
         p_value, 
         perc_urb_500
         ) %>%
  mutate(pred_urb_500 = predict(lm_500_v2, terms = "perc_urb_500"),
         part_urb_500 = pred_urb_500 + resid(lm_500_v2)) %>%
  ggplot(
    aes(
      x = perc_urb_500, 
      y = part_urb_500
      )
    ) + 
  geom_point() + 
  gghighlight(p_value < 0.05, use_direct_label = FALSE) +
  geom_hline(yintercept = 0) + 
  annotate(
    "text",
    x = 75, 
    y = 2, 
    label = rq_lab_500
  ) + 
  labs(y = "ses.MFD (partial residuals)",
       x = "% Impervious surface (500m spatial scale)") + 
  theme_bw()
)

# additional analyses ----------------------------------------------------------

## |- 250 m scale -----

reg_250 %>%
  ggplot(aes(x = Habitat_type, y = ses_mfd)) + 
  geom_boxplot() + 
  geom_point(alpha = 0.1) + 
  labs(
    x = "Urban Green Space Type", 
    y = "SES.MFD") + 
  theme_bw()

reg_250 %>%
  ggplot(aes(x = perc_urb_250)) + 
  geom_histogram() + 
  labs(x = "Percent Impervious Surface (in 250m scale)") + 
  facet_wrap(~Habitat_type) +
  theme_bw()

reg_250 %>%
  ggplot(aes(x = perc_urb_250, y = ses_mfd)) + 
  geom_point() + 
  facet_wrap(~Habitat_type) + 
  labs(
    x = "Percent Impervious Surface (in 250m scale)", 
    y = "SES.MFD") + 
  theme_bw()

## |- 500 m scale -----

table(reg_500$Habitat_type)

reg_500 %>%
  ggplot(aes(x = perc_urb_500)) + 
  geom_histogram() + 
  labs(x = "Percent Impervious Surface (in 500m scale)") + 
  facet_wrap(~Habitat_type) +
  theme_bw()

reg_500 %>%
  ggplot(aes(x = perc_urb_500, y = ses_mfd)) + 
  geom_point() + 
  facet_wrap(~Habitat_type) + 
  labs(
    x = "Percent Impervious Surface (in 500m scale)", 
    y = "SES.MFD") + 
  theme_bw()

library(emmeans)

emm_250 <- emmeans(lm_250_v2, specs = "Habitat_type")
emm_500 <- emmeans(lm_500_v2, specs = "Habitat_type")

ref_grid(lm_250_v2)
ref_grid(lm_500_v2)

pairs(emm_250)
pairs(emm_500)

# save to disk -----------------------------------------------------------------

## |- figure 4 ----
ggsave(filename = 
         here("output", "results", 
              "Xie_et_al-2021-Figure4-JAE.png"
              ),
       plot = part_urb_250, 
       device = "png", 
       width = 5, 
       height = 4
       )

## |- figure S7 ----
ggsave(filename = 
         here("output/data_appendix_output", 
              "Xie_et_al-2021-FigureS7-JAE.png"
         ),
       plot = part_urb_500, 
       device = "png", 
       width = 5, 
       height = 4
)

## |- figure S4 ----
ggsave(filename = 
         here("output/figures/supp", 
              "Xie_et_al-2021-FigureS5-JAE.png"
         ),
       plot = pairs_250, 
       device = "png", 
       width = 5, 
       height = 4
)

## |- figure S6 ----
ggsave(filename = 
         here("output/figures/supp", 
              "Xie_et_al-2021-FigureS6-JAE.png"
         ),
       plot = pairs_500, 
       device = "png", 
       width = 5, 
       height = 4
)

## |- reg mfd 250 ----
write.csv(x = reg_250,
          file = here(
            "data/final",
            "reg_mfd_250.csv")
          )

## |- reg mfd 500 ----
write.csv(x = reg_500, 
          file = here(
            "data/final",
            "reg_mfd_500.csv")
          )



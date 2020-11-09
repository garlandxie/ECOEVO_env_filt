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
library(here)    # for creating relative file-paths
library(dplyr)
library(ggplot2)
library(GGally)
library(car)
library(sp)
library(spdep)
library(patchwork)
library(gghighlight)
library(tibble)

# import -----------------------------------------------------------------------

# site info
site <- read.csv(
  here("data/original", 
       "site.csv")
  )

# ses mfd
ses_mfd <- readRDS(
  here("data/working",
       "ses_mfd.rds")
  )

# land cover
l_250 <- read.csv(
  here("data/working", 
       "land_use_250.csv")
  )

l_500 <- read.csv(
  here("data/working", 
       "land_use_500.csv")
  )

# data cleaning ----------------------------------------------------------------

# 250m 
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
         prop_tree_250, 
         prop_grass_250, 
         prop_urb_250) %>%
  
  filter(
    
    # remove sites that have only 1 species
    !is.na(ses_mfd),
    
    # remove sites that are outside TO boundary
    !is.na(prop_tree_250)  &
    !is.na(prop_grass_250) & 
    !is.na(prop_urb_250)  
  )

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
         prop_tree_500, 
         prop_grass_500, 
         prop_urb_500) %>%
  
  filter(
  
    # remove sites with only 1 species
    !is.na(ses_mfd),
    
    # remove sites outside TO boundary
    !is.na(prop_tree_500)  &
    !is.na(prop_grass_500) & 
    !is.na(prop_urb_500)  
    )

# exploratory data analysis: relationships between X and Y variables -----------

# check for clear patterns between response and explanatory relationships
# use a LOESS smoother to aid visual interpretation
# remove 95% CI's since we're still exploring data

# 250 m
reg_250 %>%
  ggplot(aes(x = prop_urb_250, y = ses_mfd)) + 
  geom_smooth(method = "loess", se = FALSE) + 
  geom_point() + 
  labs(x = "Percent impervious surface (250m)",
       y = "ses MFD")

reg_250 %>%
  ggplot(aes(x = prop_tree_250, y = ses_mfd)) + 
  geom_smooth(method = "loess", se = FALSE) + 
  geom_point() + 
  labs(x = "Percent tree cover (250m)",
       y = "ses MFD")

reg_250 %>%
  ggplot(aes(x = prop_grass_250, y = ses_mfd)) + 
  geom_smooth(method = "loess", se = FALSE) + 
  geom_point() + 
  labs(x = "Percent grass cover (250m)",
       y = "ses MFD")

# 500m
reg_500 %>%
  ggplot(aes(x = prop_urb_500, y = ses_mfd)) + 
  geom_smooth(method = "loess", se = FALSE) + 
  geom_point() + 
  labs(x = "Percent impervious surface (500m)",
       y = "ses MFD")

reg_500 %>%
  ggplot(aes(x = prop_tree_500, y = ses_mfd)) + 
  geom_smooth(method = "loess", se = FALSE) + 
  geom_point() + 
  labs(x = "Percent tree cover (500m)",
       y = "ses MFD")

reg_500 %>%
  ggplot(aes(x = prop_grass_500, y = ses_mfd)) + 
  geom_smooth(method = "loess", se = FALSE) + 
  geom_point() + 
  labs(x = "Percent grass cover (500m)",
       y = "ses MFD")

# exploratory data analysis: independent observations for Y variable -----------

reg_250 %>%
  ggplot(aes(x = longs, y = lats, col = ses_mfd)) + 
  geom_point() + 
  labs(x = "Longitude", 
       y = "Latitude") + 
  scale_colour_gradientn(colours = terrain.colors(10)) + 
  theme_minimal()

# correlation matrix -----------------------------------------------------------

# Pearson's correlation matrix for 250m
pairs_250 <- reg_250 %>%
  select(
    "% Closed Green"  = prop_tree_250, 
    "% Open Green" = prop_grass_250,
    "% Impervious" = prop_urb_250) %>%
  ggpairs() + 
  labs(title = "250m spatial scale")

# Pearson's correlation matrix for 500m
pairs_500 <- reg_500 %>%
  select(
    "% Closed Green"  = prop_tree_500, 
    "% Open Green" = prop_grass_500,
    "% Impervious" = prop_urb_500) %>%
  ggpairs() + 
  labs(title = "500m spatial scale")

# hypothesis testing: multiple regression (250m) -------------------------------

# first fit
lm_250_v1 <- lm(ses_mfd ~ prop_grass_250 + prop_tree_250 + prop_urb_250, 
             data = reg_250)
vif(lm_250_v1)

# remove tree cover
lm_250_v2 <- update(lm_250_v1, ~. -prop_tree_250)
vif(lm_250_v2)

# get summary
summary(lm_250_v2)

# show model diagnostics in a non-interactive manner 
plot(lm_250_v2, which = c(1))
plot(lm_250_v2, which = c(2))
plot(lm_250_v2, which = c(3))
plot(lm_250_v2, which = c(4))
plot(lm_250_v2, which = c(5))
plot(lm_250_v2, which = c(6))

# hypothesis testing: multiple regression (500m) -------------------------------

# first fit
lm_500_v1 <- lm(ses_mfd ~ prop_grass_500 + prop_tree_500 + prop_urb_500, 
             data = reg_500)
vif(lm_500_v1)

# remove tree cover
lm_500_v2 <- update(lm_500_v1, ~. -prop_tree_500)
vif(lm_500_v2)

# get summary
summary(lm_500_v2)

# show model diagnostics in a non-interactive manner 
plot(lm_500_v2, which = c(1))
plot(lm_500_v2, which = c(2))
plot(lm_500_v2, which = c(3))
plot(lm_500_v2, which = c(4))
plot(lm_500_v2, which = c(5))
plot(lm_500_v2, which = c(6))

# spatial autocorrelaton tests -------------------------------------------------

# 250m 
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

# 500m 
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

# plots ------------------------------------------------------------------------

# get R-squared values for labelling

rq_250 <- summary(lm_250_v2)$adj.r.squared
rq_500 <- summary(lm_500_v2)$adj.r.squared

rq_lab_250 <- bquote("Adj-R"^2: .(format(rq_250, digits = 3)))
rq_lab_500 <- bquote("Adj-R"^2: .(format(rq_500, digits = 2)))

# scatterplots: ses.MFD vs % impervious surfaces
(part_urb_250 <- reg_250 %>%
  select(site, 
         ses_mfd, 
         p_value, 
         prop_urb_250) %>%
  mutate(pred_urb_250 = predict(lm_250_v2, terms = "prop_urb_250"),
         part_urb_250 = pred_urb_250 + resid(lm_250_v2)) %>%
  ggplot(aes(x = prop_urb_250, y = part_urb_250)) + 
  geom_point() + 
  gghighlight(p_value < 0.05, use_direct_label = FALSE) + 
  geom_hline(yintercept = 0) + 
  annotate(
    "text",
    x = 75, 
    y = 2, 
    label = rq_lab_250
    ) + 
  labs(y = "ses.MFD (partial residuals)",
       x = "% Impervious surface (250m spatial scale)") + 
  theme_bw()
  )

(part_urb_500 <- reg_500 %>%
  select(
         ses_mfd, 
         p_value, 
         prop_urb_500
         ) %>%
  mutate(pred_urb_500 = predict(lm_500_v2, terms = "prop_urb_500"),
         part_urb_500 = pred_urb_500 + resid(lm_500_v2)) %>%
  ggplot(
    aes(
      x = prop_urb_500, 
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


# save to disk -----------------------------------------------------------------

ggsave(filename = 
         here("output/figures/main", 
              "fig-crit2_mfd_250.png"
              ),
       plot = part_urb_250, 
       device = "png", 
       width = 5, 
       height = 4
       )


ggsave(filename = 
         here("output/figures/supp", 
              "fig-crit2_mfd_500.png"
         ),
       plot = part_urb_500, 
       device = "png", 
       width = 5, 
       height = 4
)

ggsave(filename = 
         here("output/figures/supp", 
              "fig-corr_matrix_250.png"
         ),
       plot = pairs_250, 
       device = "png", 
       width = 5, 
       height = 4
)

ggsave(filename = 
         here("output/figures/supp", 
              "fig-corr_matrix_500.png"
         ),
       plot = pairs_500, 
       device = "png", 
       width = 5, 
       height = 4
)

write.csv(x = reg_250,
          file = here(
            "data/final",
            "reg_mfd_250.csv")
          )

write.csv(x = reg_500, 
          file = here(
            "data/final",
            "reg_mfd_500.csv")
          )



################################################################################
# Accompanying code for the paper: 
#   No environmental filtering of wasps and bees during urbanization
#
# Paper authors: 
#   Garland Xie      (1), 
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
# Purpose of this R script: to create a figure of the land cover map

# library ----------------------------------------------------------------------
library(here)      # for creating relative file-paths
library(vegan)     # for analyzing ecological community data
library(dplyr)     # for manipulating data
library(tibble)    # for manipulating data frames
library(ggplot2)   # for visualizing data
library(sf)        # for manipulating geospatial data
library(ggsn)      # for adding cartographical elements
library(patchwork) # for creating multi-panel figures
library(opendatatoronto)

# site -------------------------------------------------------------------------
site <- read.csv(
  here("data/original", 
       "site.csv")
  )

# ses mfd
comm <- read.csv(
  here("data/final", 
       "comm_matrix_B.csv"),
  row.names = 1
  )

# land use 
l_250 <- read.csv(
  here("data/working", 
       "land_use_250.csv")
  )

l_500 <- read.csv(
  here("data/working", 
       "land_use_500.csv")
  )

# TO boundary 
if(!file.exists(
  here(
    "data/original",
    "citygcs_regional_mun_wgs84.shp"
    )
  )
) {
  reg_bound_ID <- "841fb820-46d0-46ac-8dcb-d20f27e57bcc"
  packages <- show_package(reg_bound_ID)
  resources <- list_package_resources(packages)
  bound <- get_resource(resources)
}


# clean ------------------------------------------------------------------------

# species richness
SR <- comm %>%
  decostand(method = "pa") %>%
  rowSums() %>% 
  data.frame()

SR$site <- rownames(SR)
colnames(SR) <- c("ntaxa", "site_id")

# all relevant info for 250m
tidy_250 <- site %>%
  inner_join(l_250, by = c("ID" = "site")) %>%
  inner_join(SR, by = c("ID" = "site_id")) %>%
  select(ID, 
         Latitude, 
         Longitude, 
         prop_urb_250, 
         ntaxa,
         Habitat_type) 

# all relevant info for 500m
tidy_500 <- site %>%
  inner_join(l_500, by = c("ID" = "site")) %>%
  inner_join(SR, by = c("ID" = "site_id")) %>%
  select(ID, 
         Latitude, 
         Longitude, 
         prop_urb_500, 
         ntaxa,
         Habitat_type) 

# land cover map: 250m -------------------------------------------------------------------

(lc_map_250 <- ggplot(data = bound) + 
  geom_sf(fill = NA) + 
  geom_point(
    data = tidy_250, 
    aes(
      x = Longitude, 
      y = Latitude, 
      colour = prop_urb_250,
      size = ntaxa
      )
    ) +
  
  # environmental gradient
  scale_colour_gradientn(
    colours = terrain.colors(10), 
    name = "% Impervious Surface") + 
  labs(title = NULL,
       x = "Longitude", 
       y = "Latitude"
       ) + 
  
  # species richness
  scale_size_continuous(
    name   = "Species Richness",
    breaks = c(2, 5, 8, 10, 13)
  ) +
  
  # north arrow
  north(bound, symbol = 9) +
  
  # scale-bar
  scalebar(
    data = bound, 
    dist = 10 ,
    transform = TRUE, 
    dist_unit = "km",
    st.size = 2) +
  
  theme_bw() 
)

# land cover map: 500m ---------------------------------------------------------

(lc_map_500 <- ggplot(data = bound) + 
  geom_sf(fill = NA) + 
  geom_point(
    data = tidy_500, 
    aes(
      x = Longitude, 
      y = Latitude, 
      colour = prop_urb_500,
      size = ntaxa
    )
  ) +
  
   
   # environmental gradient
   scale_colour_gradientn(
     colours = terrain.colors(10), 
     name = "% Impervious Surface") + 
   labs(
        x = "Longitude", 
        y = "Latitude") + 
   
  #
  scale_size_continuous(
    name   = "Species Richness",
    breaks = c(2, 5, 8, 10, 13)
  ) +
  
  # north arrow
   
  north(bound, symbol = 9) +
  
  # scale-bar
  scalebar(
    data = bound, 
    dist = 10 ,
    transform = TRUE, 
    dist_unit = "km",
    st.size = 2) + 
  
  # theme 
  theme_bw()
)

# UGS map ----------------------------------------------------------------------

tidy_250 <- tidy_250 %>%
  mutate(Habitat_type = case_when(
    Habitat_type == "Communty" ~ "Community Garden",
    Habitat_type == "Roof"     ~ "Green Roof",
    Habitat_type == "Garden"   ~ "Home Garden",
    TRUE ~ Habitat_type 
  ))

tidy_500 <- tidy_500 %>%
  mutate(Habitat_type = case_when(
    Habitat_type == "Communty" ~ "Community Garden",
    Habitat_type == "Roof"     ~ "Green Roof",
    Habitat_type == "Garden"   ~ "Home Garden",
    TRUE ~ Habitat_type 
  ))

(ugs_map_250 <- ggplot(data = bound) + 
   
   # geometry
   geom_sf(fill = NA) + 
   geom_point(
     data = tidy_250, 
     aes(
       x = Longitude, 
       y = Latitude, 
       color = Habitat_type
     ),
     size = 2
   ) + 
   
   # scale-bar
   scalebar(
     data = bound, 
     dist = 10 ,
     transform = TRUE, 
     dist_unit = "km",
     st.size = 4)  +
   
   # labels
   labs(
     title = NULL,
     x = "Longitude",
     y = "Latitude"
   ) + 
   
   # legend
   scale_color_discrete(name = "UGS type") + 
   
   # theme
   theme_bw() +
   theme(
     legend.position = c(0.86, 0.25),
     legend.background = element_rect(fill = "white", color = "black"),
     axis.text.x = element_text(size = 14),
     axis.text.y = element_text(size = 14),
     axis.title.x = element_text(size = 16),
     axis.title.y = element_text(size = 16)
   )
)

(ugs_map_500 <- ggplot(data = bound) + 
    
    # geometry
    geom_sf(fill = NA) + 
    geom_point(
      data = tidy_500, 
      aes(
        x = Longitude, 
        y = Latitude, 
        color = Habitat_type
      ),
      size = 2,
      width = 0.01
    ) + 
    
    # north arrow
    north(bound, symbol = 9) +
    
    # scale-bar
    scalebar(
      data = bound, 
      dist = 10 ,
      transform = TRUE, 
      dist_unit = "km",
      st.size = 2)  +
    
    # labels
    labs(
      title = NULL,
      x = "Longitude",
      y = "Latitude"
    ) + 
    
    # legend
    scale_color_discrete(name = "UGS type") + 
    
    # theme
    theme_bw() + 
    theme(
      legend.position = c(0.86, 0.25),
      legend.background = element_rect(fill = "white", color = "black"),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    )
)

# save to disk -----------------------------------------------------------------

ggsave(
  plot = lc_map_250, 
  here(
  "output/figures/main", 
  "fig-lc-map_250.png"
  ),
  device = "png",
  width = 10, 
  height = 10
  )

ggsave(
  plot = lc_map_500, 
  here(
    "output/figures/supp", 
    "fig-lc-map_500.png"
  ),
  device = "png",
  width = 10, 
  height = 10
)

# note sure if UGS 250 and 500 differ so just plot both

ggsave(
  plot = ugs_map_250, 
  here(
    "output/figures/supp", 
    "fig-ugs-map_250.png"
  ),
  device = "png",
  width = 10, 
  height = 10
)

ggsave(
  plot = ugs_map_500, 
  here(
    "output/figures/supp", 
    "fig-ugs-map_500.png"
  ),
  device = "png",
  width = 10, 
  height = 10
)
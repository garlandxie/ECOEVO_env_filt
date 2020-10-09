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
# Purpose of this R script: to conduct statistical analyses for CT criteria III 

# Load libraries ---------------------------------------------------------------
library(ade4)
library(here)
library(tidyverse)
library(patchwork)
library(grid)

# Import files -----------------------------------------------------------------

# R tables - environmental data
env_250 <- read.csv(
  here("data/working", 
       "land_use_250.csv")
)

env_500 <- read.csv(
  here("data/working", 
       "land_use_500.csv")
)

# L table - species abundance
comm <- read.csv(
  here("data/final",
       "comm_matrix_B.csv"
       ),
  row.names = 1
  )


# Q table - traits
traits <- readRDS(
  here("data/final", 
       "traits.rds")
)

# RLQ: 250 ---------------------------------------------------------------------

# some prep
keep_spp <- traits %>%
  filter(!is.na(itd)) %>%
  pull(spp)

# Procedure of RLQ follows closely to:
# Dray et al. 2014. Ecology. 
# https://doi.org/10.1890/13-0196.1
# Supplementary Info 1

# correspondance analysis 
# community data matrix

L_250 <- comm %>%
  rownames_to_column(var = "site_id") %>%
  select(site_id, all_of(keep_spp)) %>%
  filter(site_id %in% env_250$site) %>%
  column_to_rownames(var = "site_id") %>%
  as.data.frame() 

dudiL_250 <- dudi.coa(df = L_250, scannf = F)

# principal component analysis
# environmental variables
# weighted by site weights from previous correspondance analysis

R_250 <- env_250 %>%
  filter(site %in% rownames(comm)) %>%
  column_to_rownames(var = "site") %>%
  select(prop_urb_250, prop_tree_250, prop_grass_250) %>%
  as.data.frame() 

dudiR_250 <- 
  dudi.pca(
    df = R_250, 
    row.w = dudiL_250$lw, 
    scannf = F
    )

# hill-smith method: PCA of discrete/continuous data 
# reference: Hill and Smith. 1976. Taxon
# weighted by sites weight from previous correspondance analysis

Q_250 <- traits %>%
  column_to_rownames(var = "spp") %>%
  select(-emer_per, -taxa) %>%
  filter(!is.na(itd)) %>%
  as.data.frame()

dudiQ_250 <- 
  dudi.hillsmith(
    df = Q_250,
    row.w = dudiL_250$cw, 
    scannf = F
    )

# RLQ: 250m spatial scale 
RLQ_250 <- rlq(dudiR = dudiR_250, 
               dudiQ = dudiQ_250, 
               dudiL = dudiL_250,
               scannf = FALSE, 
               nf = 3)

glob_RLQ_250 <- randtest(xtest = RLQ_250,
         nrepet = 49999, 
         modeltype = 6,
         alter = "greater")

# cumulative projected inertia (%)
summary(RLQ_250)

# RLQ: 500 ---------------------------------------------------------------------

# Procedure of RLQ follows closely to:
# Dray et al. 2014. Ecology. 
# https://doi.org/10.1890/13-0196.1
# Supplementary Info 1

# correspondance analysis 
# community data matrix
L_500 <- comm %>%
  rownames_to_column(var = "site_id") %>%
  select(site_id, all_of(keep_spp)) %>%
  filter(site_id %in% env_500$site) %>%
  column_to_rownames(var = "site_id") %>%
  as.data.frame() 

dudiL_500 <- dudi.coa(
  df = L_500, 
  scannf = F
  )

# principal component analysis
# environmental variables
# weighted by site weights from previous correspondance analysis

R_500 <- env_500 %>%
  filter(site %in% rownames(comm)) %>%
  column_to_rownames(var = "site") %>%
  select(prop_urb_500, prop_tree_500, prop_grass_500) %>%
  as.data.frame()

dudiR_500 <- 
  dudi.pca(
    df = R_500, 
    row.w = dudiL_500$lw, 
    scannf = F
  )

# hill-smith method: PCA of discrete/continuous data 
# reference: Hill and Smith. 1976. Taxon
# weighted by sites weight from previous correspondance analysis

Q_500 <- traits %>%
  column_to_rownames(var = "spp") %>%
  select(-emer_per, -taxa) %>%
  filter(!is.na(itd)) %>%
  as.data.frame()

dudiQ_500 <- 
  dudi.hillsmith(
    df = Q_500,
    row.w = dudiL_500$cw, 
    scannf = F
  )

# RLQ: 250m spatial scale 
# obtain first two RLQ axes in a non-iteractive manner 
RLQ_500 <- rlq(dudiR = dudiR_500, 
               dudiQ = dudiQ_500, 
               dudiL = dudiL_500,
               scannf = FALSE, 
               nf = 2)

glob_RLQ_500 <- 
  randtest(
    xtest = RLQ_500,
    nrepet = 49999, 
    modeltype = 6,
    alter = "greater"
    )

# cumulative projected inertia (%)
summary(RLQ_500)

# plots: scree plot ------------------------------------------------------------

eigs_250 <- data.frame(
  Axes = paste("Axis", 1:3, sep = " "),
  Eigenvalues =  RLQ_250$eig,
  Projected.Inertia = round(RLQ_250$eig/sum(RLQ_250$eig)*100, digits = 2)
)

eigs_500 <- data.frame(
  Axes = paste("Axis", 1:3, sep = " "),
  Eigenvalues =  RLQ_500$eig,
  Projected.Inertia = round(RLQ_500$eig/sum(RLQ_500$eig)*100, digits = 2)
)

scree_250 <- 
  ggplot(data = eigs_250, aes(x = Axes, y = Projected.Inertia)) +
  geom_point() + 
  geom_segment(
    aes(
      x = Axes, 
      xend = Axes,
      y = 0, 
      yend = Projected.Inertia) 
  ) + 
  ylim(0, 100) + 
  labs(title = "A) 250m spatial scale",
       x = "RLQ Axes", 
       y = "Projected Inertia (%)") + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
  ) 

scree_500 <- 
  ggplot(data = eigs_500, aes(x = Axes, y = Projected.Inertia)) +
  geom_point() + 
  geom_segment(
    aes(
      x = Axes, 
      xend = Axes,
      y = 0, 
      yend = Projected.Inertia) 
  ) + 
  ylim(0, 100) + 
  labs(title = "B) 500m spatial scale",
       x = "RLQ Axes", 
       y = "Projected Inertia (%)") + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
  ) 

scree <- scree_250 + scree_500

# plots: env loadings  ---------------------------------------------------------

# set up grob objects to place text underneath the plots

more_urb <- textGrob(
  "More Urban", 
  gp = gpar(
    fontsize = 13
  )
)

less_urb <- textGrob(
  "Less Urban", 
  gp = gpar(
    fontsize = 13
  )
)

less_urb_arrow <- linesGrob(
  arrow = arrow(
    type   = "open", 
    ends   = "first", 
    length = unit(3, "mm")
    ),
  gp = gpar(
    col = "black", 
    lwd = 1)
  )
more_urb_arrow <- linesGrob(
  arrow = arrow(
    type   = "open", 
    ends   = "last", 
    length = unit(3, "mm")
  ),
  gp = gpar(
    col = "black", 
    lwd = 1)
)
  
R_250_load <- RLQ_250$l1 %>%
  rownames_to_column(var = "class") %>%
  select(class, RS1) %>%
  mutate(
    
    class = case_when(
      class == "prop_urb_250"   ~ "% Impervious surface", 
      class == "prop_tree_250"  ~ "% Tree cover",
      class == "prop_grass_250" ~ "% Grass cover"),
    
    class = factor(class),
    class = fct_reorder(class, RS1)) %>%
  ggplot(aes(x = class, y = RS1, fill = class)) + 
  geom_point() + 
  geom_segment(
    aes(
      x = class, 
      xend = class,
      y = 0, 
      yend = RS1)
    )+ 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(title = "A) 250m spatial scale",
       x = NULL,
       y = "Relative importance in environmnental scores") + 
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_text(size = 10))

R_500_load <- RLQ_500$l1 %>%
  rownames_to_column(var = "class") %>%
  select(class, RS1) %>%
  ggplot(aes(x = class, y = RS1, fill = class)) + 
  geom_point() + 
  geom_segment(
    aes(
      x = class, 
      xend = class,
      y = 0, 
      yend = RS1)
  )+   
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(title = "B) 500m spatial scale",
       x = NULL,
       y = "Relative importance in environmnental scores") + 
  coord_flip() +
  theme_bw() + 
  theme(axis.text.y = element_blank(), 
        legend.position = "none",
        axis.title.x = element_text(size = 10))     

R_load <- R_250_load + R_500_load

# plots: species scores ---------------------------------------------------------------

(L_250_load <- RLQ_250$mQ %>%
  rownames_to_column(var = "species") %>%
  select(species, NorS1) %>%
  mutate(
    species = 
      str_replace(
        species, 
        pattern = "_",
        replacement = " "
        ),
    
    species = reorder(species, NorS1)
    )%>%
  ggplot(aes(y = species, x = NorS1)) + 
  geom_point() +
  geom_segment(
    aes(
      x = 0, 
      xend = NorS1,
      y = species, 
      yend = species
    ), 
    alpha = 0.5) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  labs(title = "250m spatial scale",
       x = "Relative importance in normed species scores",
       y = NULL) + 
  theme_bw() +
  theme(plot.margin  = unit(c(5, 5, 5, 5), "lines"),
        axis.text.y  = element_text(size = 13),
        axis.title.x = element_text(size = 13)) + 
   annotation_custom(
     more_urb,
     ymin = -5,
     ymax = -5,
     xmin = 2,
     xmax = 2
   ) + 
   annotation_custom(
     less_urb,
     ymin = -5,
     ymax = -5,
     xmin = -2,
     xmax = -2
   ) + 
   annotation_custom(
     less_urb_arrow,
     ymin = -5,
     ymax = -5,
     xmin = -1,
     xmax = 0
   ) +
   annotation_custom(
     more_urb_arrow,
     ymin = -5,
     ymax = -5,
     xmin = 1,
     xmax = 0
   ) + 
   coord_cartesian(clip = "off") 
)

L_500_load <- RLQ_500$mQ %>%
  rownames_to_column(var = "species") %>%
  select(species, NorS1) %>%
  mutate(
    species = 
      str_replace(
        species, 
        pattern = "_",
        replacement = " "
      ),
    
    species = reorder(species, NorS1)
    ) %>%
  ggplot(aes(x = NorS1, y = species)) + 
  geom_point() +
  geom_segment(
    aes(
      x    = 0, 
      xend = NorS1,
      y    = species, 
      yend = species
    ), 
    alpha = 0.5) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  labs(title = "500m spatial scale",
       x = "Relative importance in normed species scores",
       y = NULL) + 
  theme_bw() + 
  theme(plot.margin  = unit(c(5, 5, 5, 5), "lines"),
        axis.text.y  = element_text(size = 13),
        axis.title.x = element_text(size = 13)) + 
  annotation_custom(
    more_urb,
    ymin = -5,
    ymax = -5,
    xmin = 2,
    xmax = 2
  ) + 
  annotation_custom(
    less_urb,
    ymin = -5,
    ymax = -5,
    xmin = -2,
    xmax = -2
  ) + 
  annotation_custom(
    less_urb_arrow,
    ymin = -5,
    ymax = -5,
    xmin = -1,
    xmax = 0
  ) +
  annotation_custom(
    more_urb_arrow,
    ymin = -5,
    ymax = -5,
    xmin = 1,
    xmax = 0
  ) + 
  coord_cartesian(clip = "off") 

L <- L_250_load + L_500_load

# plot: trait scores -----------------------------------------------------------------

(RLQ_250_load <- RLQ_250$c1 %>%
  rownames_to_column(var = "traits")%>%
  select(traits, CS1) %>%
  mutate(traits = as.character(traits), 
         traits = case_when(
           traits == "volt.2"             ~ "Multivoltinism",
           traits == "volt.1"             ~ "Univoltinism",
           traits == "nest.secretions"   ~ "Nest (Secretions)", 
           traits == "nest.scrapings"    ~ "Nest (Scrapings)",
           traits == "nest.resin"        ~ "Nest (Resin)",
           traits == "nest.mud"          ~ "Nest (Mud)",
           traits == "nest.leaf.and.mud" ~ "Nest (Leaf + Mud)",
           traits == "nest.leaf"         ~ "Nest (Leaf)", 
           traits == "nest.grass"        ~ "Nest (Grass)",
           traits == "nativ.1"          ~ "Native Status", 
           traits == "nativ.0"          ~ "Exotic Status",
           traits == "itd"               ~ "Body Size (ITD)",
           traits == "diet.spider"       ~ "Diet (Spider)",
           traits == "diet.pollen"       ~ "Diet (Pollen)", 
           traits == "diet.katydid"      ~ "Diet (Katydid)", 
           traits == "diet.caterpillar"  ~ "Diet (Caterpillar)", 
           traits == "diet.beetle.larva" ~ "Diet (Beetle Larva)", 
           traits == "diet.aphid"        ~ "Diet (Aphid)", 
           TRUE ~ traits)
         ) %>%  
  ggplot(aes(
    x = CS1,
    y = traits %>% fct_reorder(CS1)) 
    ) + 
  geom_point(shape = 18, size = 2) +
  geom_segment(
    aes(
      y    = traits, 
      yend = traits,
      x    = 0, 
      xend = CS1
    ), 
    alpha = 0.5) + 
  xlim(-3, 3) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  labs(title = "250m spatial scale",
       y = NULL,
       x = "Relative importance in trait scores") + 
  theme_bw() + 
  theme(plot.margin  = unit(c(5, 5, 5, 5), "lines"),
        axis.text    = element_text(size = 13),
        axis.title.x = element_text(size = 13)) + 
  annotation_custom(
    more_urb,
    ymin = -3,
    ymax = -3,
    xmin = 3,
    xmax = 3
  ) + 
  annotation_custom(
    less_urb,
    ymin = -3,
    ymax = -3,
    xmin = -3,
    xmax = -3
  ) + 
   annotation_custom(
     less_urb_arrow,
     ymin = -3,
     ymax = -3,
     xmin = -2,
     xmax = 0
   ) +
   annotation_custom(
     more_urb_arrow,
     ymin = -3,
     ymax = -3,
     xmin = 2,
     xmax = 0
   ) + 
  coord_cartesian(clip = "off") 
)

(RLQ_500_load <- RLQ_500$c1 %>%
  rownames_to_column(var = "traits")%>%
  select(traits, CS1) %>%
  mutate(traits = as.character(traits), 
         traits = case_when(
           traits == "volt.2"             ~ "Multivoltinism",
           traits == "volt.1"             ~ "Univoltinism",
           traits == "nest.secretions"   ~ "Nest (Secretions)", 
           traits == "nest.scrapings"    ~ "Nest (Scrapings)",
           traits == "nest.resin"        ~ "Nest (Resin)",
           traits == "nest.mud"          ~ "Nest (Mud)",
           traits == "nest.leaf.and.mud" ~ "Nest (Leaf + Mud)",
           traits == "nest.leaf"         ~ "Nest (Leaf)", 
           traits == "nest.grass"        ~ "Nest (Grass)",
           traits == "nativ.1"          ~ "Native Status", 
           traits == "nativ.0"          ~ "Exotic Status",
           traits == "itd"               ~ "Body Size (ITD)",
           traits == "diet.spider"       ~ "Diet (Spider)",
           traits == "diet.pollen"       ~ "Diet (Pollen)", 
           traits == "diet.katydid"      ~ "Diet (Katydid)", 
           traits == "diet.caterpillar"  ~ "Diet (Caterpillar)", 
           traits == "diet.beetle.larva" ~ "Diet (Beetle Larva)", 
           traits == "diet.aphid"        ~ "Diet (Aphid)", 
           TRUE ~ traits)
  ) %>%  
  ggplot(aes(y = traits %>% fct_reorder(CS1), x = CS1)) + 
  geom_point(shape = 18, size = 2.5) +
  geom_segment(
    aes(
      y    = traits, 
      yend = traits,
      x    = 0, 
      xend = CS1
    ), 
    alpha = 0.5) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  xlim(-3, 3) +   
  labs(title = "500m spatial scale",
       y = NULL,
       x = "Relative importance in trait scores") + 
  theme_bw() + 
  theme(plot.margin = unit(c(5, 5, 5, 5), "lines"),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13)
        ) + 
  annotation_custom(
    more_urb,
    ymin = -3,
    ymax = -3,
    xmin = 3,
    xmax = 3
  ) + 
  annotation_custom(
      less_urb,
      ymin = -3,
      ymax = -3,
      xmin = -3,
      xmax = -3
    ) + 
  annotation_custom(
      less_urb_arrow,
      ymin = -3,
      ymax = -3,
      xmin = -2,
      xmax = 0
    ) +
  annotation_custom(
      more_urb_arrow,
      ymin = -3,
      ymax = -3,
      xmin = 2,
      xmax = 0
    ) +
  coord_cartesian(clip = "off") 
)  
 
# save to disk -----------------------------------------------------------------

# PCA environmental variable loadings
ggsave(
  plot = R_load, 
  filename = here(
    "output/figures/supp", 
    "fig_sup_env_loadings.png"
    ),
  device = "png",
  height = 5, 
  width = 8
)

# Scree plots
ggsave(
  plot = scree, 
  filename = here(
    "output/figures/supp", 
    "fig_sup_scree.png"),
  device = "png",
  height = 5, 
  width = 7
)

# species scores
ggsave(
  plot = L_250_load, 
  filename = here(
    "output/figures/supp", 
    "fig_sup_species_rlq_250.png"),
  device = "png",
  height = 10, 
  width = 9
)

ggsave(
  plot = L_500_load, 
  filename = here(
    "output/figures/supp", 
    "fig_sup_species_rlq_500.png"),
  device = "png",
  height = 10, 
  width = 9
)

# trait scores
ggsave(
  plot = RLQ_250_load, 
  filename = here(
    "output/figures/supp", 
    "fig_sup_traits_250.png"),
  device = "png",
  height = 7, 
  width = 9
)

ggsave(
  plot = RLQ_500_load, 
  filename = here(
    "output/figures/supp", 
    "fig_sup_traits_500.png"),
  device = "png",
  height = 7, 
  width = 9
)
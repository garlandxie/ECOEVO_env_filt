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
       "traits_no_volt.rds")
)

# RLQ: 250 ---------------------------------------------------------------------

# some prep
keep_spp <- traits %>%
  filter(!is.na(body_size)) %>%
  filter(species != "Hylaeus_punctatus") %>%
  pull(species)

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
  select(perc_urb_250, perc_tree_250, perc_grass_250) %>%
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
  filter(species != "Hylaeus_punctatus") %>%
  column_to_rownames(var = "species") %>%
  filter(!is.na(body_size)) %>%
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
  select(perc_urb_500, perc_tree_500, perc_grass_500) %>%
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
  filter(species != "Hylaeus_punctatus") %>%
  column_to_rownames(var = "species") %>%
  filter(!is.na(body_size)) %>%
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

# prep: grobs ------------------------------------------------------------------

# set up grob objects to place text underneath the plots

more_urb <- textGrob(
  "More Urban", 
  gp = gpar(
    fontsize = 13
  )
)

less_urb <- textGrob(
  "More Closed Green", 
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

# plots: env loadings  ---------------------------------------------------------
  
R_250_load <- RLQ_250$l1 %>%
  rownames_to_column(var = "class") %>%
  select(class, RS1) %>%
  mutate(
    
    class = case_when(
      class == "perc_urb_250"   ~ "% Impervious surface", 
      class == "perc_tree_250"  ~ "% Closed green cover",
      class == "perc_grass_250" ~ "% Open green cover"),
    
    class = factor(class),
    class = fct_reorder(class, RS1)) %>%
  ggplot(aes(x = RS1, y = class)) + 
  geom_point() + 
  geom_segment(
    aes(
      x = 0, 
      xend = RS1,
      y = class, 
      yend = class)
    ) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  labs(title = "250 m spatial scale",
       x = "Relative importance in environmnental scores",
       y = NULL) + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.margin  = unit(c(5, 5, 5, 5), "lines"),
        axis.title.x = element_text(size = 10)) + 
  annotation_custom(
    more_urb,
    ymin = -1,
    ymax = -1,
    xmin = 0.5,
    xmax = 0.5
  ) + 
  annotation_custom(
    less_urb,
    ymin = -1,
    ymax = -1,
    xmin = -0.5,
    xmax = -0.5
  ) + 
  annotation_custom(
    less_urb_arrow,
    ymin = -1,
    ymax = -1,
    xmin = -0.2,
    xmax = 0
  ) +
  annotation_custom(
    more_urb_arrow,
    ymin = -1,
    ymax = -1,
    xmin = 0.3,
    xmax = 0
  ) + 
  coord_cartesian(clip = "off") 

R_500_load <- RLQ_500$l1 %>%
  rownames_to_column(var = "class") %>%
  select(class, RS1) %>%
  mutate(
    
    class = case_when(
      class == "perc_urb_500"   ~ "% Impervious surface", 
      class == "perc_tree_500"  ~ "% Closed green cover",
      class == "perc_grass_500" ~ "% Open green cover"),
    
    class = factor(class),
    class = fct_reorder(class, RS1)) %>%
  ggplot(aes(x = RS1, y = class)) + 
  geom_point() + 
  geom_segment(
    aes(
      x = 0, 
      xend = RS1,
      y = class, 
      yend = class)
  )+   
  geom_vline(xintercept = 0, linetype = "dashed") + 
  labs(title = "500 m spatial scale ",
       x = "Relative importance in environmnental scores",
       y = NULL) + 
  theme_bw() + 
  theme( 
        legend.position = "none",
        plot.margin  = unit(c(5, 5, 5, 5), "lines"),
        axis.title.x = element_text(size = 10)) + 
  annotation_custom(
    more_urb,
    ymin = -1,
    ymax = -1,
    xmin = 0.5,
    xmax = 0.5
  ) + 
  annotation_custom(
    less_urb,
    ymin = -1,
    ymax = -1,
    xmin = -0.5,
    xmax = -0.5
  ) + 
  annotation_custom(
    less_urb_arrow,
    ymin = -1,
    ymax = -1,
    xmin = -0.2,
    xmax = 0
  ) +
  annotation_custom(
    more_urb_arrow,
    ymin = -1,
    ymax = -1,
    xmin = 0.3,
    xmax = 0
  ) + 
  coord_cartesian(clip = "off") 

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
  xlim(-3, 3) + 
  labs(title = "250 m spatial scale",
       x = "Relative importance in normed species scores",
       y = NULL) + 
  theme_bw() +
  theme(plot.margin  = unit(c(5, 5, 5, 5), "lines"),
        axis.text.y  = element_text(size = 13, face = "italic"),
        axis.title.x = element_text(size = 13)) + 
   annotation_custom(
     more_urb,
     ymin = -5,
     ymax = -5,
     xmin = 3,
     xmax = 3
   ) + 
   annotation_custom(
     less_urb,
     ymin = -5,
     ymax = -5,
     xmin = -3,
     xmax = -3
   ) + 
   annotation_custom(
     less_urb_arrow,
     ymin = -5,
     ymax = -5,
     xmin = -1.5,
     xmax = 0
   ) +
   annotation_custom(
     more_urb_arrow,
     ymin = -5,
     ymax = -5,
     xmin = 2,
     xmax = 0
   ) + 
   coord_cartesian(clip = "off") 
)

(L_500_load <- RLQ_500$mQ %>%
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
  xlim(-3, 3) + 
  labs(title = "500 m spatial scale",
       x = "Relative importance in normed species scores",
       y = NULL) + 
  theme_bw() + 
  theme(plot.margin  = unit(c(5, 5, 5, 5), "lines"),
        axis.text.y  = element_text(size = 13, face = "italic"),
        axis.title.x = element_text(size = 13)) + 
  annotation_custom(
    more_urb,
    ymin = -5,
    ymax = -5,
    xmin = 1,
    xmax = 1
  ) + 
  annotation_custom(
    less_urb,
    ymin = -5,
    ymax = -5,
    xmin = -3,
    xmax = -3
  ) + 
  annotation_custom(
    less_urb_arrow,
    ymin = -5,
    ymax = -5,
    xmin = -2,
    xmax = 0
  ) +
  annotation_custom(
    more_urb_arrow,
    ymin = -5,
    ymax = -5,
    xmin = 0.2,
    xmax = 0
  ) + 
  coord_cartesian(clip = "off") 
)


# plot: trait scores -----------------------------------------------------------------

(RLQ_250_load <- RLQ_250$c1 %>%
  rownames_to_column(var = "traits")%>%
  select(traits, CS1) %>%
  mutate(
    traits = as.character(traits), 
    taxon = case_when(
      
      # Nesting material preference
      traits == "nesti.Secretions"              ~ "Bees", 
      traits == "nesti.Nest.tube.scrapings"     ~ "Bees",
      traits == "nesti.Resin"                   ~ "Bees",
      traits == "nesti.Mud"                     ~ "Bees",
      traits == "nesti.Leaf.cut"                ~ "Bees",
      traits == "nesti.Leaf.hair"               ~ "Bees", 
      traits == "nesti.Grass"                   ~ "Bees",
      traits == "nesti.Leaf.pulp"               ~ "Bees",
      
      # Body size
      traits == "body_size"                     ~ "Both",
      
      # Origin
      traits == "origi.Holarctic"               ~ "Both",
      traits == "origi.Nearctic"                ~ "Both",
      traits == "origi.Palearctic"              ~ "Both",
      
      # Diet
      traits == "prima.Aphids"                  ~ "Wasps",
      traits == "prima.Beetle.larva"            ~ "Wasps",
      traits == "prima.Caterpillars"            ~ "Wasps",
      traits == "prima.Pollen"                  ~ "Wasps",
      traits == "prima.Spiders"                 ~ "Wasps",
      traits == "prima.Tree.crickets"           ~ "Wasps",
      
      # Specialization
      traits == "speci.Family..Aphididae."      ~ "Wasps", 
      traits == "speci.Family..Asteraceae."     ~ "Bees",
      traits == "speci.Family..Campanulaceae."  ~ "Bees",
      traits == "speci.Family..Chrysomelidae."  ~ "Wasps",
      traits == "speci.Family.Aphididae."       ~ "Wasps",
      traits == "speci.Genus..Campanula."       ~ "Bees",
      traits == "speci.Multi.Order"             ~ "Both",
      traits == "speci.Multi.Order..Coleoptera..Lepidoptera." ~ "Wasps",
      traits == "speci.Order..Araneae."         ~ "Wasps",
      traits == "speci.Order..Lepidoptera."     ~ "Wasps",
      traits == "speci.Order..Orthoptera."      ~ "Wasps",
      
      # Trophic level
      traits == "troph.First"  ~ "Bees",
      traits == "troph.Second" ~ "Wasps",
      traits == "troph.Third"  ~ "Wasps",
      
      TRUE ~ traits)
  ) %>%
    mutate(traits = case_when(
           
          # Nesting material preference
          traits == "nesti.Secretions"              ~ "Nest (Secretions)", 
          traits == "nesti.Nest.tube.scrapings"    ~ "Nest (Scrapings)",
          traits == "nesti.Resin"                   ~ "Nest (Resin)",
          traits == "nesti.Mud"                     ~ "Nest (Mud)",
          traits == "nesti.Leaf.cut"               ~ "Nest (Leaf cut)",
          traits == "nesti.Leaf.hair"               ~ "Nest (Leaf hair)", 
          traits == "nesti.Grass"                  ~ "Nest (Grass)",
          traits == "nesti.Leaf.pulp"              ~ "Nest (pulp)",
          
          # Body size
          traits == "body_size"                    ~ "Body Size (ITD)",
           
          # Origin
          traits == "origi.Holarctic"              ~ "Origin (Holarctic)",
          traits == "origi.Nearctic"               ~ "Origin (Nearctic)",
          traits == "origi.Palearctic"             ~ "Origin (Palearctic)",
           
          # Diet
          traits == "prima.Aphids"                 ~ "Diet (Alphids)",
          traits == "prima.Beetle.larva"           ~ "Diet (Beetle larvae)",
          traits == "prima.Caterpillars"           ~ "Diet (Caterpillars)",
          traits == "prima.Pollen"                 ~ "Diet (Pollen)",
          traits == "prima.Spiders"                ~ "Diet (Spiders)",
          traits == "prima.Tree.crickets"          ~ "Diet (Tree crickets)",
          
          # Specialization
          traits == "speci.Family..Aphididae."      ~ "Specialization (Aphididae)", 
          traits == "speci.Family..Asteraceae."     ~ "Specialization (Asteraceae)",
          traits == "speci.Family..Campanulaceae."  ~ "Specialization (Campanulaceae)",
          traits == "speci.Family..Chrysomelidae."  ~ "Specialization (Chrysomelidae)",
          traits == "speci.Family.Aphididae."       ~  "Specialization (Aphididae)",
          traits == "speci.Genus..Campanula."       ~ "Specialization (Campanula)",
          traits == "speci.Multi.Order"             ~ "Specialization (Multi-order)",
          traits == "speci.Multi.Order..Coleoptera..Lepidoptera." ~ "Specialization (Multi-order: Coleoptera + Lepidoptera)",
          traits == "speci.Order..Araneae."         ~ "Specialization (Aranae)",
          traits == "speci.Order..Lepidoptera."     ~ "Specialization (Lepidoptera)",
          traits == "speci.Order..Orthoptera."      ~ "Specialization (Orthoptera)",
          
          # Trophic level
          traits == "troph.First"  ~ "Trophic Rank (First)",
          traits == "troph.Second" ~ "Trophic Rank (Second)",
          traits == "troph.Third"  ~ "Trophic Rank (Third)",
       
          TRUE ~ traits)
          ) %>%  
  ggplot(aes(
    x = CS1,
    y = traits %>% fct_reorder(CS1)) 
    ) + 
  geom_point(aes(shape = taxon), size = 2) +
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
  labs(title = "250 m spatial scale",
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
    xmin = 3.5,
    xmax = 3.5
  ) + 
  annotation_custom(
    less_urb,
    ymin = -3,
    ymax = -3,
    xmin = -3.5,
    xmax = -3.5
  ) + 
   annotation_custom(
     less_urb_arrow,
     ymin = -3,
     ymax = -3,
     xmin = -1.5,
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
    mutate(
      traits = as.character(traits), 
      taxon = case_when(
        
        # Nesting material preference
        traits == "nesti.Secretions"              ~ "Bees", 
        traits == "nesti.Nest.tube.scrapings"     ~ "Bees",
        traits == "nesti.Resin"                   ~ "Bees",
        traits == "nesti.Mud"                     ~ "Bees",
        traits == "nesti.Leaf.cut"                ~ "Bees",
        traits == "nesti.Leaf.hair"               ~ "Bees", 
        traits == "nesti.Grass"                   ~ "Bees",
        traits == "nesti.Leaf.pulp"               ~ "Bees",
        
        # Body size
        traits == "body_size"                     ~ "Both",
        
        # Origin
        traits == "origi.Holarctic"               ~ "Both",
        traits == "origi.Nearctic"                ~ "Both",
        traits == "origi.Palearctic"              ~ "Both",
        
        # Diet
        traits == "prima.Aphids"                  ~ "Wasps",
        traits == "prima.Beetle.larva"            ~ "Wasps",
        traits == "prima.Caterpillars"            ~ "Wasps",
        traits == "prima.Pollen"                  ~ "Wasps",
        traits == "prima.Spiders"                 ~ "Wasps",
        traits == "prima.Tree.crickets"           ~ "Wasps",
        
        # Specialization
        traits == "speci.Family..Aphididae."      ~ "Wasps", 
        traits == "speci.Family..Asteraceae."     ~ "Bees",
        traits == "speci.Family..Campanulaceae."  ~ "Bees",
        traits == "speci.Family..Chrysomelidae."  ~ "Wasps",
        traits == "speci.Family.Aphididae."       ~ "Wasps",
        traits == "speci.Genus..Campanula."       ~ "Bees",
        traits == "speci.Multi.Order"             ~ "Both",
        traits == "speci.Multi.Order..Coleoptera..Lepidoptera." ~ "Wasps",
        traits == "speci.Order..Araneae."         ~ "Wasps",
        traits == "speci.Order..Lepidoptera."     ~ "Wasps",
        traits == "speci.Order..Orthoptera."      ~ "Wasps",
        
        # Trophic level
        traits == "troph.First"  ~ "Bees",
        traits == "troph.Second" ~ "Wasps",
        traits == "troph.Third"  ~ "Wasps",
        
        TRUE ~ traits)
    ) %>%
    mutate(traits = case_when(
      
      # Nesting material preference
      traits == "nesti.Secretions"              ~ "Nest (Secretions)", 
      traits == "nesti.Nest.tube.scrapings"    ~ "Nest (Scrapings)",
      traits == "nesti.Resin"                   ~ "Nest (Resin)",
      traits == "nesti.Mud"                     ~ "Nest (Mud)",
      traits == "nesti.Leaf.cut"               ~ "Nest (Leaf cut)",
      traits == "nesti.Leaf.hair"               ~ "Nest (Leaf hair)", 
      traits == "nesti.Grass"                  ~ "Nest (Grass)",
      traits == "nesti.Leaf.pulp"              ~ "Nest (pulp)",
      
      # Body size
      traits == "body_size"                    ~ "Body Size (ITD)",
      
      # Origin
      traits == "origi.Holarctic"              ~ "Origin (Holarctic)",
      traits == "origi.Nearctic"               ~ "Origin (Nearctic)",
      traits == "origi.Palearctic"             ~ "Origin (Palearctic)",
      
      # Diet
      traits == "prima.Aphids"                 ~ "Diet (Alphids)",
      traits == "prima.Beetle.larva"           ~ "Diet (Beetle larvae)",
      traits == "prima.Caterpillars"           ~ "Diet (Caterpillars)",
      traits == "prima.Pollen"                 ~ "Diet (Pollen)",
      traits == "prima.Spiders"                ~ "Diet (Spiders)",
      traits == "prima.Tree.crickets"          ~ "Diet (Tree crickets)",
      
      # Specialization
      traits == "speci.Family..Aphididae."      ~ "Specialization (Aphididae)", 
      traits == "speci.Family..Asteraceae."     ~ "Specialization (Asteraceae)",
      traits == "speci.Family..Campanulaceae."  ~ "Specialization (Campanulaceae)",
      traits == "speci.Family..Chrysomelidae."  ~ "Specialization (Chrysomelidae)",
      traits == "speci.Family.Aphididae."       ~  "Specialization (Aphididae)",
      traits == "speci.Genus..Campanula."       ~ "Specialization (Campanula)",
      traits == "speci.Multi.Order"             ~ "Specialization (Multi-order)",
      traits == "speci.Multi.Order..Coleoptera..Lepidoptera." ~ "Specialization (Multi-order: Coleoptera + Lepidoptera)",
      traits == "speci.Order..Araneae."         ~ "Specialization (Aranae)",
      traits == "speci.Order..Lepidoptera."     ~ "Specialization (Lepidoptera)",
      traits == "speci.Order..Orthoptera."      ~ "Specialization (Orthoptera)",
      
      # Trophic level
      traits == "troph.First"  ~ "Trophic Rank (First)",
      traits == "troph.Second" ~ "Trophic Rank (Second)",
      traits == "troph.Third"  ~ "Trophic Rank (Third)",
      
      TRUE ~ traits)
    ) %>%  
  ggplot(aes(y = traits %>% fct_reorder(CS1), x = CS1)) + 
  geom_point(aes(shape = taxon), size = 2.5) +
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
  labs(title = "500 m spatial scale",
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
    xmin = 3.5,
    xmax = 3.5
  ) + 
  annotation_custom(
      less_urb,
      ymin = -3,
      ymax = -3,
      xmin = -3.5,
      xmax = -3.5
    ) + 
  annotation_custom(
      less_urb_arrow,
      ymin = -3,
      ymax = -3,
      xmin = -1.5,
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
  plot = R_250_load, 
  filename = here(
    "output/figures/supp", 
    "fig_sup_env_loadings_250.png"
    ),
  device = "png",
  height = 5, 
  width = 8
)

ggsave(
  plot = R_500_load, 
  filename = here(
    "output/figures/supp", 
    "fig_sup_env_loadings_500.png"
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
  height = 8, 
  width = 11
)

ggsave(
  plot = RLQ_500_load, 
  filename = here(
    "output/figures/supp", 
    "fig_sup_traits_500.png"),
  device = "png",
  height = 8, 
  width = 11
)
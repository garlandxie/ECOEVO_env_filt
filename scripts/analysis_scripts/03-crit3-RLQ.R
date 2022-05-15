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
library(readr)

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
traits <- readr::read_csv(
  here("output", "tables", "traits_tidy.csv"), 
  col_types = cols(
    species          = col_character(),
    body_size        = col_double(), 
    nesting_material = col_factor(), 
    primary_diet     = col_factor(), 
    specialization   = col_factor(), 
    trophic_rank     = col_factor(), 
    native_status    = col_factor(), 
    num_nest_mat     = col_factor()
  ) 
  
)

# cheack packaging -------------------------------------------------------------

glimpse(comm)
glimpse(env_250)
glimpse(env_500)
glimpse(traits)

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

# Figure S8: scree plots (250m and 50mm scales) --------------------------------

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

# Prep: grobs!! ----------------------------------------------------------------

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

# Figure S9: Environmental loadings (250m) -------------------------------------
  
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

# Figure S10: Environmental loadings (500m scale) -------------------------------

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

# Figure 6: species scores (250m) ----------------------------------------------

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
  xlim(-4, 4) + 
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
     xmin = -3.2,
     xmax = -3.2
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

# Figure S12: species score (500 m) -------------------------------------------- 

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
  xlim(-7, 7) + 
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
    xmin = 4,
    xmax = 4
  ) + 
  annotation_custom(
    less_urb,
    ymin = -5,
    ymax = -5,
    xmin = -5,
    xmax = -5
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
    xmin = 2,
    xmax = 0
  ) + 
  coord_cartesian(clip = "off") 
)

# Figure 5: Trait scores (250 m) -----------------------------------------------

RLQ_250_load <- RLQ_250$c1 %>%
  rownames_to_column(var = "traits") %>%
  select(traits, CS1) %>%
   
  # obtain "trait" category 
  # i.e., nesting material, primary diet, specialization, trophic level
  # i.e., numbering of nesting materials, native status
  mutate(
    nest_cat = case_when(
      str_detect(traits, pattern = "body")   ~ "Body", 
      str_detect(traits, pattern = "nesti.") ~ "Nesti",
      str_detect(traits, pattern = "prima.") ~ "Prima", 
      str_detect(traits, pattern = "speci.") ~ "Speci",
      str_detect(traits, pattern = "troph.") ~ "Troph", 
      str_detect(traits, pattern = "nativ.") ~ "Native", 
      str_detect(traits, pattern = "num.n.") ~ "Num",
      TRUE ~ traits
    )
  ) %>%
   
  # group into sections by trait category 
  # then order by important score within each group
  # (instead of all results being globally ordered by important score)
  group_by(nest_cat) %>%
  arrange(nest_cat, CS1) %>%
   
  mutate(
    trait = as.character(traits), 
    taxon = case_when(
      
      # Nesting material preference
      traits == "nesti.Nest.tube.scrapings"     ~ "Bees",
      traits == "nesti.Leaf.hair"               ~ "Bees", 
      traits == "nesti.Mud"                     ~ "Bees",
      traits == "nesti.Resin"                   ~ "Both",
      traits == "nesti.Leaf.pulp"               ~ "Bees",
      traits == "nesti.Leaf.cut"                ~ "Bees",
      traits == "nesti.Secretions"              ~ "Both", 
      traits == "nesti.Grass"                   ~ "Wasps",

      # Body size
      traits == "body_size"                     ~ "Both",
      
      # Origin
      traits == "nativ.Native"                  ~ "Both",
      traits == "nativ.Non.Native"              ~ "Both",

      # Diet
      traits == "prima.Aphids"                  ~ "Wasps",
      traits == "prima.Beetle.larva"            ~ "Wasps",
      traits == "prima.Caterpillars"            ~ "Wasps",
      traits == "prima.Pollen"                  ~ "Bees",
      traits == "prima.Single.Spider"           ~ "Wasps",
      traits == "prima.Spiders"                 ~ "Wasps",
      traits == "prima.Tree.crickets"           ~ "Wasps",
      
      # Specialization
      traits == "speci.Family"                 ~ "Both",
      traits == "speci.Order"                  ~ "Both",
      traits == "speci.Multi.Order"            ~ "Both",
    
      # Trophic level
      traits == "troph.First"                  ~ "Bees",
      traits == "troph.Second"                 ~ "Wasps",
      traits == "troph.Third"                  ~ "Wasps",
      
      # Number of nesting materials
      traits == "num_n.Single"                 ~ "Both",
      traits == "num_n.Multi"                  ~ "Both", 
      traits == "num_n.None"                   ~ "Both",
      
      TRUE ~ traits)
    
  ) %>%
   
  mutate(traits = case_when(
     
     # Nesting material preference
     traits == "nesti.Nest.tube.scrapings"     ~ "Nest (Scrapings)",
     traits == "nesti.Leaf.hair"               ~ "Nest (Leaf hair)", 
     traits == "nesti.Mud"                     ~ "Nest (Mud)",
     traits == "nesti.Resin"                   ~ "Nest (Resin)",
     traits == "nesti.Leaf.pulp"               ~ "Nest (Pulp)",
     traits == "nesti.Leaf.cut"                ~ "Nest (Leaf cut)",
     traits == "nesti.Secretions"              ~ "Nest (Secretions)", 
     traits == "nesti.Grass"                   ~ "Nest (Grass)",

     # Body size
     traits == "body_size"                     ~ "Body Size (ITD)",
     
     # Origin
     traits == "nativ.Native"                  ~ "Native Status",
     traits == "nativ.Non.Native"              ~ "Exotic Status",
     
     # Diet
     traits == "prima.Pollen"                  ~ "Diet (Pollen)",
     traits == "prima.Caterpillars"            ~ "Diet (Caterpillars)",
     traits == "prima.Beetle.larva"            ~ "Diet (Beetle larvae)",
     traits == "prima.Single.Spider"           ~ "Diet (Single Spider)", 
     traits == "prima.Tree.crickets"           ~ "Diet (Tree crickets)",
     traits == "prima.Aphids"                  ~ "Diet (Aphids)",
     traits == "prima.Spiders"                 ~ "Diet (Spiders)",

     # Specialization
     traits == "speci.Multi.Order"             ~ "Specialization (Multi-Order)",
     traits == "speci.Family"                  ~ "Specialization (Family)",
     traits == "speci.Order"                   ~ "Specialization (Order)",

     # Trophic level
     traits == "troph.Second"                  ~ "Trophic Rank (Second)",
     traits == "troph.Third"                   ~ "Trophic Rank (Third)",
     traits == "troph.First"                   ~ "Trophic Rank (First)",

     # Number of nesting materials
     traits == "num_n.Single"                 ~ "Nesting Material Type Collected (Single)",
     traits == "num_n.Multi"                  ~ "Nesting Material Type Collected (Multiple)", 
     traits == "num_n.None"                   ~ "Nesting Material Type Collected (None)",
     
     TRUE ~ traits)
     ) 
 
level_order <- factor(RLQ_250_load$traits, levels = RLQ_250_load$traits)

RLQ_250_plot <- RLQ_250_load %>%
  ggplot(aes(x = CS1, y = factor(traits, levels = level_order))) + 
  #geom_point(aes(shape = taxon), size = 2) +
  geom_segment(
    aes(
      y    = traits, 
      yend = traits,
      x    = 0, 
      xend = CS1
    ), 
    alpha = 0.5) + 
  xlim(-4, 4) + 
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

# Figure S11: Trait scores (500m scale) ----------------------------------------

(RLQ_500_load <- RLQ_500$c1 %>%
  rownames_to_column(var = "traits")%>%
  select(traits, CS1) %>%
   mutate(
     traits = as.character(traits), 
     taxon = case_when(
       
       # Nesting material preference
       traits == "nesti.Secretions"              ~ "Both", 
       traits == "nesti.Nest.tube.scrapings"     ~ "Bees",
       traits == "nesti.Resin"                   ~ "Both",
       traits == "nesti.Mud"                     ~ "Bees",
       traits == "nesti.Leaf.cut"                ~ "Bees",
       traits == "nesti.Leaf.hair"               ~ "Bees", 
       traits == "nesti.Grass"                   ~ "Wasps",
       traits == "nesti.Leaf.pulp"               ~ "Bees",
       
       # Body size
       traits == "body_size"                     ~ "Both",
       
       # Origin
       traits == "nativ.Native"                  ~ "Both",
       traits == "nativ.Non.Native"              ~ "Both",
       
       # Diet
       traits == "prima.Aphids"                  ~ "Wasps",
       traits == "prima.Beetle.larva"            ~ "Wasps",
       traits == "prima.Caterpillars"            ~ "Wasps",
       traits == "prima.Pollen"                  ~ "Bees",
       traits == "prima.Single.Spider"           ~ "Wasps",
       traits == "prima.Spiders"                 ~ "Wasps",
       traits == "prima.Tree.crickets"           ~ "Wasps",
       
       # Specialization
       traits == "speci.Family"                 ~ "Both",
       traits == "speci.Order"                  ~ "Both",
       traits == "speci.Multi.Order"            ~ "Both",
       
       # Trophic level
       traits == "troph.First"                  ~ "Bees",
       traits == "troph.Second"                 ~ "Wasps",
       traits == "troph.Third"                  ~ "Wasps",
       
       # Number of nesting materials
       traits == "num_n.Single"                 ~ "Both",
       traits == "num_n.Multi"                  ~ "Both", 
       traits == "num_n.None"                   ~ "Both",
       
       TRUE ~ traits)
   ) %>%
   mutate(traits = case_when(
     
     # Nesting material preference
     traits == "nesti.Secretions"              ~ "Nest (Secretions)", 
     traits == "nesti.Nest.tube.scrapings"     ~ "Nest (Scrapings)",
     traits == "nesti.Resin"                   ~ "Nest (Resin)",
     traits == "nesti.Mud"                     ~ "Nest (Mud)",
     traits == "nesti.Leaf.cut"                ~ "Nest (Leaf cut)",
     traits == "nesti.Leaf.hair"               ~ "Nest (Leaf hair)", 
     traits == "nesti.Grass"                   ~ "Nest (Grass)",
     traits == "nesti.Leaf.pulp"               ~ "Nest (Pulp)",
     
     # Body size
     traits == "body_size"                     ~ "Body Size (ITD)",
     
     # Origin
     traits == "nativ.Native"                  ~ "Native Status",
     traits == "nativ.Non.Native"              ~ "Exotic Status",
     
     # Diet
     traits == "prima.Aphids"                  ~ "Diet (Aphids)",
     traits == "prima.Beetle.larva"            ~ "Diet (Beetle larvae)",
     traits == "prima.Caterpillars"            ~ "Diet (Caterpillars)",
     traits == "prima.Pollen"                  ~ "Diet (Pollen)",
     traits == "prima.Single.Spider"           ~ 'Diet (Single Spider)', 
     traits == "prima.Spiders"                 ~ "Diet (Spiders)",
     traits == "prima.Tree.crickets"           ~ "Diet (Tree crickets)",
     
     # Specialization
     traits == "speci.Family"                  ~ "Specialization (Family)",
     traits == "speci.Order"                   ~ "Specialization (Order)",
     traits == "speci.Multi.Order"             ~ "Specialization (Multi-Order)",
     
     # Trophic level
     traits == "troph.First"                   ~ "Trophic Rank (First)",
     traits == "troph.Second"                  ~ "Trophic Rank (Second)",
     traits == "troph.Third"                   ~ "Trophic Rank (Third)",
     
     # Number of nesting materials
     traits == "num_n.Single"                 ~ "Nesting Material Type Collected (Single)",
     traits == "num_n.Multi"                  ~ "Nesting Material Type Collected (Multiple)", 
     traits == "num_n.None"                   ~ "Nesting Material Type Collected (None)",
     
     
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
  xlim(-7, 7) +   
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
      xmin = -3.7,
      xmax = -3.7
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
 
# Other analyses ---------------------------------------------------------------

# 250 m spatial scale
# fourthcorner(
#  R_250, 
#  L_250, 
#  Q_250, 
#  modeltype = 6,
#  p.adjust.method.G = "fdr", 
#  p.adjust.method.D = "fdr",
#  nrepet = 49999
#  )

# 500 m spatial scale
# fourthcorner(
#  R_500, 
#  L_500, 
#  Q_500, 
#  modeltype = 6,
#  p.adjust.method.G = "fdr", 
#  p.adjust.method.D = "fdr",
#  nrepet = 49999
#  )

# Save To disk: main figures ---------------------------------------------------

# Figure 5
ggsave(
  plot = RLQ_250_load, 
  filename = here(
    "output/figures/main", 
    "Xie_et_al-2021-Figure5-JAE.png"),
  device = "png",
  height = 8, 
  width = 11
)

# Figure 6
ggsave(
  plot = L_250_load, 
  filename = here(
    "output/figures/main", 
    "Xie_et_al-2021-Figure6-JAE.png"),
  device = "png",
  height = 10, 
  width = 9
)

# Save to disk: supplementary figures ------------------------------------------

# Figure S8
ggsave(
  plot = scree, 
  filename = here(
    "output/figures/supp", 
    "Xie_et_al-2021-FigureS8-JAE.png"),
  device = "png",
  height = 5, 
  width = 7
)

# Figure S9
ggsave(
  plot = R_250_load, 
  filename = here(
    "output/figures/supp", 
    "Xie_et_al-2021-FigureS9-JAE.png"
    ),
  device = "png",
  height = 5, 
  width = 8
)

# Figure S10
ggsave(
  plot = R_500_load, 
  filename = here(
    "output/figures/supp", 
    "Xie_et_al-2021-FigureS10-JAE.png"
  ),
  device = "png",
  height = 5, 
  width = 8
)

# Figure S11
ggsave(
  plot = RLQ_500_load, 
  filename = here(
    "output/figures/supp", 
    "Xie_et_al-2021-FigureS11-JAE.png"),
  device = "png",
  height = 8, 
  width = 11
)

# Figure S12
ggsave(
  plot = L_500_load, 
  filename = here(
    "output/figures/supp", 
    "Xie_et_al-2021-FigureS12-JAE.png"),
  device = "png",
  height = 10, 
  width = 9
)
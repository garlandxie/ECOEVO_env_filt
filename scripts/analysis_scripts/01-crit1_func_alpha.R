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
# Purpose of this R script: to conduct statistical analyses for CT criteria I 

# libraries --------------------------------------------------------------------
library(StatMatch)   # for calculating gower's distance
library(picante)     # for analyzing community matrices
library(here)        # for creating relative file-paths
library(vegan)       # for calculating relative abundances
library(tidyverse)  
library(gghighlight) # for visualising data

# import -----------------------------------------------------------------------

comm <- read.csv(
  here(
    "data", "analysis_data", 
    "comm_matrix_B.csv"
    ),
  row.names = 1
  )


trait <- read.csv(
  here(
    "data", 
    "analysis_data", 
    "traits_tidy.csv"
    )
  )

site <- read.csv(
  here(
    "data", "input_data", 
    "site_data.csv"
    )
)

# functional distance matrix  --------------------------------------------------

keep_spp <- trait %>%
  filter(!is.na(its)) %>%
  filter(species != "Hylaeus_punctatus") %>%
  pull(species)

comm_rel <- comm %>%
  select(all_of(keep_spp)) %>%
  decostand(method = "total")
  
trait_tidy <- trait %>%
  as.data.frame() %>% 
  filter(species != "Hylaeus_punctatus") %>%
  column_to_rownames(var = "species") %>%
  filter(!is.na(body_size)) 

trait_dist <- gower.dist(trait_tidy)
colnames(trait_dist) <- rownames(trait_tidy)
rownames(trait_dist) <- rownames(trait_tidy)

# ses.mfd: ---------------------------------------------------------------------

set.seed(123)

ses_mfd  <- ses.mpd(
  samp = comm_rel,  
  dis  = trait_dist, 
  runs = 4999, 
  abundance.weighted = TRUE,
  null.model = "frequency"
  )

colnames(ses_mfd) <- 
  str_replace(
    colnames(ses_mfd), 
    "mpd",
    "mfd"
  )

# grobs: set up arrows ---------------------------------------------------------

cluster <- textGrob(
  "Clustering", 
  gp = gpar(
    fontsize = 12
  )
)

over_disp <- textGrob(
  "Overdispersion", 
  gp = gpar(
    fontsize = 12
  )
)

less_arrow <- linesGrob(
  arrow = arrow(
    type   = "open", 
    ends   = "first", 
    length = unit(3, "mm")
  ),
  gp = gpar(
    col = "black", 
    lwd = 1)
)

more_arrow <- linesGrob(
  arrow = arrow(
    type   = "open", 
    ends   = "last", 
    length = unit(3, "mm")
  ),
  gp = gpar(
    col = "black", 
    lwd = 1)
)

# Figure 3: CT Criteria I (evidence of clustering) -----------------------------

ses_mfd_tidy <- ses_mfd %>%
  rename(ses_mfd = mfd.obs.z, 
         p_value = mfd.obs.p) %>%
  filter(ntaxa != 1 & 
         ntaxa != 0 & 
         !is.na(p_value)) %>% 
  rownames_to_column(var = "site_id")  

ses_mfd_habitat <- ses_mfd_tidy %>%
  inner_join(site, by = c("site_id" = "ID")) %>%
  select(
    site_id, 
    ntaxa,
    mfd.obs, 
    mfd.rand.mean,
    mfd.rand.sd,
    mfd.obs.rank, 
    ses_mfd, 
    p_value, 
    Habitat_type)

(ses_mfd_ugs <- ses_mfd_habitat %>%
    mutate(Habitat_type = case_when(
      Habitat_type == "Community" ~ "Community Garden",
      Habitat_type == "Roof"      ~ "Green Roof",
      Habitat_type == "Garden"    ~ "Home Garden",
      Habitat_type == "Park"      ~ "Public Park",
      TRUE ~ Habitat_type )
      ) %>%
  ggplot(aes(y = Habitat_type, x = ses_mfd)) + 
  geom_point() + 
  gghighlight(
    p_value < 0.05, 
    use_direct_label = FALSE,
    use_group_by = FALSE
    ) +
  xlim(-4, 4) + 
  labs(
    x = "ses.MFD",
    y = "UGS Type"
  ) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_bw() +
  theme(plot.margin = unit(c(5, 5, 5, 5), "lines")) + 
  annotation_custom(
    cluster,
    xmin = -3,
    xmax = -3,
    ymin = -0.5,
    ymax = -0.5
    ) +
  annotation_custom(
    less_arrow,
    xmin = 0,
    xmax = -1.5,
    ymin = -0.5,
    ymax = -0.5
    ) + 
  annotation_custom(
    over_disp,
    xmin = 3,
    xmax = 3,
    ymin = -0.5,
    ymax = -0.5
  ) + 
  annotation_custom(
    more_arrow,
    xmin = 0,
    xmax = 1.3,
    ymin = -0.5,
    ymax = -0.5
    ) +
  coord_cartesian(clip = "off")
)
    
# save to disk -----------------------------------------------------------------

# mfd
saveRDS(
  ses_mfd_tidy, 
  here("data/working",
       "ses_mfd.rds"
  )
)

write.csv(
  ses_mfd_tidy, 
  here("output/tables/supp",
       "ses_mfd.csv"
       )
)

# plots
ggsave(filename = 
         here(
           "output/figures/main",
           "Xie_et_al-2021-Figure3-JAE.png"
         ),
       plot = ses_mfd_ugs,
       width  = 7, 
       height = 5,
       device = "png"
)



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
    "data/final", 
    "comm_matrix_B.csv"
    ),
  row.names = 1
  )

trait <- readRDS(
  here("data/final", 
       "traits.rds"
       )
  )

# functional distance matrix  --------------------------------------------------

keep_spp <- trait %>%
  filter(!is.na(itd)) %>%
  pull(spp)

comm_rel <- comm %>%
  select(all_of(keep_spp)) %>%
  decostand(method = "total")
  
trait_tidy <- trait %>%
  as.data.frame() %>% 
  column_to_rownames(var = "spp") %>%
  filter(!is.na(itd)) %>%
  select(native_y_n, 
         nest, 
         diet, 
         volt, 
         itd) 

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

# how many sites have 1 species?
sum(ses_mfd$ntaxa != 1)

# plots ------------------------------------------------------------------------

ses_mfd_tidy <- ses_mfd %>%
  rename(ses_mfd = mfd.obs.z, 
         p_value = mfd.obs.p) %>%
  filter(ntaxa != 1 & 
         ntaxa != 0 & 
         !is.na(p_value)) %>% 
  rownames_to_column(var = "site_id") %>%
  filter(!(site_id %in% c("Dumesh",
                        "Kavanah", 
                        "Lynott", 
                        "RangersGround", 
                        "RangersRoof",
                        "Chute")))

(crit1_mfd <- ses_mfd_tidy %>%
  ggplot(aes(x = p_value, y = ses_mfd)) +
  geom_point() + 
  gghighlight(p_value < 0.05) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  geom_rect(ymin = -4, ymax = 0,
            xmin = 0, xmax = 0.05,
            alpha = 0.01,
            fill = "red") + 
  ylim(-4, 2) + 
  labs(x = "p-value", 
       y= "ses.MFD") + 
  theme_bw()
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
           "fig-crit1.png"
         ),
       plot = crit1_mfd,
       width  = 4, 
       height = 4,
       device = "png"
       )





# Bees and wasps in the city  

The following documentation follows the [TIER 4.0 Protocol](https://www.projecttier.org/tier-protocol/protocol-4-0/root/).

## Manuscript 

Title: No evidence of environmental filtering of cavity-nesting solitary bees and wasps by urbanization using trap nests

Authors: Garland Xie, Nicholas Sookhan, Kelly A. Carscadden, and Scott MacIvor

Accepted in Ecology and Evolution (ECE-2022-02-00197.R1)

## Abstract

Spatial patterns in biodiversity are used to establish conservation priorities and ecosystem management plans. The environmental filtering of communities along urbanization gradients has been used to explain biodiversity patterns but demonstrating filtering requires precise statistical tests to link suboptimal environments at one end of a gradient to lower population sizes via ecological traits. Here we employ a three-part framework on observational community data to test: I) for trait clustering (i.e., phenotypic similarities among co-occurring species) by comparing trait diversity to null expectations, II) if trait clustering is correlated with an urbanization gradient, and III) if species’ traits relate to environmental conditions. If all criteria are met, then there is evidence that urbanization is filtering communities based on their traits. We use a community of 46 solitary cavity-nesting bee and wasp species sampled across Toronto, a large metropolitan city, over three years to test these hypotheses. None of the criteria were met, so we did not have evidence for environmental filtering. We do show that certain ecological traits influence which species perform well in urban environments. For example, Cellophane bees (Hylaeus: Colletidae) secrete their own nesting material and were overrepresented in urban areas, while native leafcutting bees (Megachile: Megachilidae) were most common in greener areas. For wasps, prey preference was important, with aphid-collecting (Psenulus and Passaloecus: Crabronidae) and generalist spider-collecting (Trypoxylon: Crabronidae) wasps overrepresented in urban areas and caterpillar- and beetle-collecting wasps (Euodynerus and Symmorphus: Vespidae, respectively) overrepresented in greener areas. We emphasize that changes in the prevalence of different traits across urban gradients without corresponding changes in trait diversity with urbanization, do not constitute environmental filtering. By applying this rigorous framework, future studies can test whether urbanization filters other nesting guilds (i.e., ground-nesting bees and wasps) or larger communities consisting of entire taxonomic groups. 


## Software and version

Code to reproducible the analysis and figures is written in the R programming language (version 4.1.1). 
The developer (Garland Xie) typically writes and runs the code using macOS Big Sur 11.6.

## Folder structure 

- Data
  - [Input data](data/input_data): all raw data 
    - [site.csv](data/input_data/site_data.csv): site information (i.e., survey years and urban green space type). Geographic coordinates were removed for this version of the dataset due to sensitivity regarding private households. 
    - [traits.csv](data/input_data/traits.csv): trait matrix of seven different ecological traits for a given wasp or bee taxa. 
    - [trap_nest.csv](data/input_data/trap_nest.csv): data on completed, dead, or alive brood cells for each species within a given tube across three sampled years (2011, 2012, 2013)
    - [metadata](data/input_data/metadata): containing information about the sources and contents of input data files
  - [Analysis data](data/analysis_data): data that is in a format suitable for the analysis   
  - [Intermediate data](data/intermediate_data): partially processed data to use in subsequent data analyses
  
- Scripts
  - [Processing scripts](scripts/processing_scripts): The commands in these scripts transform Input Data Files into Analysis Data Files
  - [Analysis scripts](scripts/analysis_scripts): The commands in these scripts generate main results in the manuscript
  - [Master script](scripts/master_script.R): script that reproduces the Results of your project by executing all the other scripts, in the correct order

- Output
  -  [Data appendix output](output/data_appendix_output): contains files with figures, tables and other output presented in the Data Appendix
  -  [Results](output/results): contains files with figures, tables, and other output presented in your report

## Instructions for reproducing the results

**Important note:** one of our datasets (i.e., site information) contains sensitive information regarding the latitude and longitude of private households. For this reason, we have not made this dataset available for public use. In addition, this data is of minimal importance for reproducing most of statistical analyses (with the exception of reproducing the spatial autocorrelation tests, maps, and calculating land cover metrics). If necessary, please contact the corresponding author, Garland Xie (garlandxie@gmail.com) to request access this particular dataset. 

To use the code in this repository to reproduce the manuscript's results,
please follow the following steps:
1. `git clone` this repository or download it as a zip folder
2. Open `Rstudio`, go to `file > Open project` and open the `D99_env_filt.Rproj`
Rproject associated with this repository
3. Run `renv::restore()` in your R console. Requires `renv` package (see [THIS](https://rstudio.github.io/renv/articles/renv.html) vignette)

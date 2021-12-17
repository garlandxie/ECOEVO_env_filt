# Bees and wasps in the city  

The following documentation of this README file follows the TIER 4.0 Protocol (https://www.projecttier.org/tier-protocol/protocol-4-0/root/)

## Manuscript 

Ecological traits, not level of urbanization, structure cavity-nesting bee and wasp diversity in the city.

Authors: Garland Xie, Nicholas Sookhan, Kelly A. Carscadden, and Scott MacIvor

Submitted to Journal of Animal Ecology (In Review)

## Abstract

1. Spatial patterns in biodiversity shed light on community assembly and can be used to prioritize conservation and ecosystem management. If urbanization leads to more functionally homogeneous solitary bee and wasp communities, then there should be impacts to critical ecosystem services, and in turn, urban food production and stability of remnant native vegetation. 

2. The ‘environmental filtering’ of communities along urbanization gradients has been used to explain biodiversity patterns. However, demonstrating environmental filtering requires precise statistical tests to link suboptimal environments at one end of a gradient to lower population sizes via ecological traits. 

3. Here, we employ a rigorous three-part framework to test I) for clustering by comparing standardized effect sizes of trait diversity to null expectations, II) if trait clustering is correlated to the urbanization gradient, and III) if species’ traits relate to environmental conditions. If all criteria are met, then there is evidence that urbanization is filtering communities based on their traits. We use a community of 51 solitary cavity-nesting bee and wasp species sampled across a large metropolitan city over three years to test for environmental filtering by urbanization. 

4. None of the criteria were met, so we did not have evidence for environmental filtering. We do however demonstrate that certain ecological traits influence which species perform well in the extreme ends of the urbanization gradient. Cellophane bees (Hylaeus: Colletidae) secrete their own nesting material and were overrepresented in urban areas, while native leafcutting bees (Megachile: Megachilidae), which line nests with cut or chewed plant material, were most common in more tree canopied (less urban) areas. For wasps, prey preference was important, with aphid-collecting (Psenulus and Passaloecus: Crabronidae) and generalist spider-collecting (Trypoxylon: Crabronidae) wasps overrepresented in urban areas and caterpillar- and beetle-collecting wasps (Euodynerus and Symmorphus: Vespidae, respectively) overrepresented towards greener areas. Non-native species were positively associated with urban areas. 

5. We emphasize that changes in the prevalence of different traits across urban gradients, without corresponding changes in trait diversity with urbanization, do not constitute environmental filtering. Nonetheless, ecological traits offer novel and practical insights into the relationship between urbanization and biodiversity.


## Software and version

Code to reproducible the analysis and figures is written in the R programming language (version 4.1.1). 
The developer (Garland Xie) typically writes and runs the code using macOS Big Sur 11.6

## Folder structure 

- Report
- Data
  - [InputData](data/input_data): all raw data 
    - [Metadata](data/input_data/metadata): containing information about the sources and contents of InputData files
  - [AnalysisData](data/analysis_data): data that is in a format suitable for the analysis   
  - [IntermediateData](data/intermediate_data): partially processed data to use in subsequent data analyses
- Scripts
  - [ProcessingScripts](scripts/processing_scripts): The commands in these scripts transform Input Data Files into Analysis Data Files
  - [DataAppendixScripts](scripts/data_appendix_scripts): The commands in theses scripts produce the figures, tables, and descriptive statistics in the Data Appendix
  - [AnalysisScripts](scripts/analysis_scripts): The commands in these scripts generate the results in the Rmarkdown report
- Document: contains the master script that reproduces the Results of your project by executing all the other scripts, in the correct order
- Output
  -  [DataAppendixOutput](output/data_appendix_output): contains files with figures, tables and other output presented in the Data Appendix
  -  [Results](output/results): contains files with figures, tables, and other output presented in your report

## Instructions for reproducing the results

To use the code in this repository to reproduce the manuscript's results,
please follow the following steps:
1. `git clone` this repository or download it as a zip folder
2. Open `Rstudio`, go to `file > Open project` and open the `Msc.Rproj`
Rproject associated with this repository
3. Run `renv::restore()` in your R console. Requires `renv` package (see [THIS](https://rstudio.github.io/renv/articles/renv.html) vignette)

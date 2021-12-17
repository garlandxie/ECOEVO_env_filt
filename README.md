# Bees and wasps in the city  

The following documentation of this README file follows the TIER 4.0 Protocol (https://www.projecttier.org/tier-protocol/protocol-4-0/root/)

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

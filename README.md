# Blue Whales and Lagrangian Features
 This repository contains the code and processed data to generate the analysis and figures in the manuscript *Blue whales increase feeding rates at fine-scale ocean features* (currently in review). Authors: James A. Fahlbush, Max F. Czapanskiy, John Calambokidis, David E. Cade, Briana Abrahms, Elliott L. Hazen and Jeremy A. Goldbogen

### How to cite

Following publication, please cite this compendium as:

> Fahlbusch, J. A., Czapanskiy, M. F., Calambokidis, J., Cade, D. E.,
> Abrahms, B., Hazen, E. L., and Goldbogen, J. A. (2022). *Compendium
> of R code and data for “Blue whales increase feeding rates at fine-scale ocean features”*. Accessed 28 Jan 2022.

A knitted R-Markdown with all of the Methods and Results (including the code to produce them) can be found at:
https://physalus.github.io/Blue-Whales-and-Lagrangian-Features/

All main text figures (1-5) and supplemental figures are generated using the scripts found in the main directory and data found in the folders in this repository, however most are then brought into Adobe Illustrator for final figure composition. 

The repository is organized as follows:
* [:file\_folder: Main directory](https://github.com/physalus/Blue-Whales-and-Lagrangian-Features): all scripts 
  * BlueWhaleLagrangianMethods.Rmd contains all of the analyses for the the manuscript and produces a PDF file when knitted
  * LagrangianMethods.Rmd is identical to BlueWhaleLagrangianMethods.Rmd but produces the HTML file found at the link above 
  * Appendix1_Blue_Whale_Data_Processing.R contains the processing steps for blue whale data including extracting FTLE (Requires the sub-directory "dataRaw")
  * Appendix2_Null_Model_Data_Processing.R contains the processing steps for creating null model data including extracting FTLE (Requires the sub-directory "dataRaw")
  * Figure1_BlueWhaleDives.R produces the panels for figure 1 in the manuscript
  * Figure2_MapTracersFTLE.R produces the panels for figure 2 in the manuscript     
* [:file\_folder: dataProcessed](https://github.com/physalus/Blue-Whales-and-Lagrangian-Features/tree/main/dataProcessed): all processed data used in the analyses of this manuscript
  * Files include the outputs from Appendix1_Blue_Whale_Data_Processing.R and Appendix2_Null_Model_Data_Processing.R, and are used as inputs to BlueWhaleLagrangianMethods.Rmd
* [:file\_folder: Figure_files](https://github.com/physalus/Blue-Whales-and-Lagrangian-Features/tree/main/Figure_files): all additional data used to produce figures 1 and 2 of this manuscript
* [:file\_folder: rmdFiles](https://github.com/physalus/Blue-Whales-and-Lagrangian-Features/tree/main/rmdFiles): Files associated with the creation of the knitted R-markdown PDF and HTML files 
* [:file\_folder: Output](https://github.com/physalus/Blue-Whales-and-Lagrangian-Features/tree/main/Output): Plots that are exported during the knitting process and data precessing steps
* [:file\_folder: functions](https://github.com/physalus/Blue-Whales-and-Lagrangian-Features/tree/main/functions): Custom functions used this analysis
* [:file\_folder: PackageReferences](https://github.com/physalus/Blue-Whales-and-Lagrangian-Features/tree/main/PackageReferences): bibtex citations for publicly-available software packages used in this analysis

*Note: Due to space limitations, the following directories are incomplete on github. The files to populate these folders can be found at: https://purl.stanford.edu/kn066mv7984*
* [:file\_folder: dataRaw](https://github.com/physalus/Blue-Whales-and-Lagrangian-Features/tree/main/dataRaw): raw tag-associated data and surface current data (including restored and processed FTLE)
  * *These files are only necessary if you want to recreate the processed data. Processed data can be found in the folder "dataProcessed"*
  * DeploymentMetadata.csv contains information about each of the blue whale deployments 
  * FTLEMetadata.xlsx contains the settings used to calculate FTLE for each deployment
* [:file\_folder: Models](https://github.com/physalus/Blue-Whales-and-Lagrangian-Features/tree/main/Models): contains the processed GLMM models in .rds format.  
  * *These files are the Generalized Linear Mixed Effects Models produced in the analysis and ake a considerable amount of time and processing power to produce.*

### Licenses

**Code, Text and Figures :** [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)

**Data :** [ODbL-1.0 Open Database License](https://opendatacommons.org/licenses/dbcl/1-0/)


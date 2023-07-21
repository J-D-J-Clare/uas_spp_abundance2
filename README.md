#### Authors:  

Andrii Zaiats, T. Trevor Caughlin, Jennyffer Cruz, David S. Pilliod, Megan E. Cattau, Rongsong Liu, Richard Rachman, Maisha Maliha, Donna Delparte, John D. J. Clare

---  

#### Overview:

This repository stores the scripts that were used for the article "Scalable monitoring of plant population abundance in rangelands with aerial imagery under high observation errors."  

___  

The repository contains data processing, tabular data, figures, modeling results, and simulations presented in the manuscript. 

Folders:\
    - **utils**: scripts for data pre-processing. Require external data.
    - **data**: contains .csv\ files and pre-processed spatial data that were used in the final statistical models.\
    - **figures**: includes .png/ figures from _figures.R_
    
R scripts:\
    - *00_helper_fn.R*: contains utility functions.  
    - *01_data_prep.R*: data preparations; require external data.  
    - *02_figures.R*: generates figures for the manuscript.  
    - *03_explore_nimblefit.R*: initial explorations of the model fit objects.  
    - *04_results.R*: derived quantities from the model fit objects presented in the paper.  
    - *05_power_results.R*: generates figures and results from the power analysis.

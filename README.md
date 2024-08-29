#### Authors:  

Andrii Zaiats, T. Trevor Caughlin, Jennyffer Cruz, David S. Pilliod, Megan E. Cattau, Rongsong Liu, Richard Rachman, Maisha Maliha, Donna Delparte, John D. J. Clare

---  

#### Overview:

This repository stores the scripts that were used for the article "Propagating observation errors to enable scalable and rigorous enumeration of plant population abundance with aerial imagery."  

___  

The repository contains data processing, tabular data, figures, modeling results, and simulations presented in the manuscript. 

Folders:\
    - **utils**: scripts for reproducibility of the data pre-processing steps. Require external data.  
    - **data**: contains .csv\ files and pre-processed spatial data that were used in the final statistical models.  
    - **figures**: includes .png/ figures from _figures.R_  
    - **model_scripts**: final R scripts for the NIMBLE models.  
    - **results**: summaries of the posterior distributions stored in .csv files.  
    - **Simulation scripts**: simulation R scripts for sapmpling scenarios 1 and 2 (see main text for details).  
    
    
R scripts:\
    - *00_helper_fn.R*: contains utility functions.  
    - *01_data_prep.R*: data preparations; require external data.  
    - *02_figures.R*: generates figures for the manuscript.  
    - *03_explore_nimblefit.R*: initial explorations of the model fit objects and extraction of posterior summaries into **results**.  
    - *04_results.R*: derived quantities from the model fit objects presented in the manuscript.  
    - *05_power_results.R*: generates figures and results from the power analysis.

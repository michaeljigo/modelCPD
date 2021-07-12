These directories contain the model, behavioral data and analyses associated with the paper:

            Michael Jigo, David J. Heeger & Marisa Carrasco; An image-computable model on how 
            endogenous and exogenous attention differentially alter visual perception.
            Proceedings of the National Academy of Sciences (in press)


-----------------------------------------
---------- Recreate figures -------------
-----------------------------------------

The main figures in the paper (Figures 4-9) can be recreated using the functions in the directory "create_figures".
Each function is titled with the figure it will recreate. 
Functions may take a few minutes to run because it will call and evaluate the model on several stimulus images.


-----------------------------------------
----------------- Model -----------------
-----------------------------------------

The model code is located in the "model" directory and is named "imAmodel" - IMage-computable Attention Model.
To generate population responses to input images, "imAmodel" calls "compute_model_response" which
computes the model's response to the image.


-----------------------------------------
---------- Directory structure ----------
-----------------------------------------

Directories are structured as follows:
Key: |--  = directory branch
     |-> = description of files contained in directory

|-- create_figures
|   |-> functions that recreate the main figures in the paper
|   |    
|-- data
|   |-- behavior
|   |   |-> .mat files contain behavioral data for each texture segmentation experiment used in the paper
|   |   |-> the README file provides more details
|   |    
|   |-- bootstrap_samples
|   |   |-> sub-directories contain bootstrap samples for each model variant fit to the data
|   |    
|-- generalize
|   |-> directory contains data for two separate experiments and model code to fit the data from each experiment
|   |    
|   |-- JigoCarrasco2020
|   |   |-- code
|   |   |   |-> functions that perform fit to behavioral data
|   |   | 
|   |   |-- data
|   |   |   |-- behavior
|   |   |   |-- fitted_parameters
|   |   | 
|   |   |-- figures
|   |   |   |-> preliminary figures for fits
|   |    
|   |-- MontagnaPestilliCarrasco2009
|   |   |-> organized similarly to JigoCarrasco2020
|   |    
|-- model
|   |-> contains functions that execute model and also fit model to texture segmentation experiments
|   |    
|   |-- helperfun
|   |   |-> contains functions that facilitate model execution
|   |    
|-- texture
|   |-> contains functions that will create texture, Landolt square and grating stimuli

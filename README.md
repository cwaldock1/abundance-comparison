# __abundance-comparison__

This README contains brief descriptions of the input data and the scripts used in the comparison of niche-abundance models. 


# __ecological survey data__

## Reef Life Survey

raw-data/RLS-spatial-fish-data-extract _provided by Rick Stuart-Smith on the 08.08.2019_

raw-data/RLS-site-subset-metadata _accessed from reeflifesurvey.org on 19.09.2019 for the purpose of providing and australian subset of data from which to assess preliminary models_

## Breeding Bird Survey

raw-data/breeding-bird-survey-usa/50-StopData/...
raw data downloaded via ftp server from the BBS website https://www.pwrc.usgs.gov/BBS/RawData/Choose-Method.cfm folder contains additional files descriptions and metadata.

## Environmental and spatial datasets

Documented in table XX in supporting materials.

# __scripts__

## __overall script organisation__: 

Scripts are organised into the following set of folders:  
  
__data-processing__: This set of scripts processes the raw survey data, subsets to our focal species and provides the inputs to the model scripts.   
__fitting-models__: This set of scripts provides a configuration script with calls a long script of model calls (model-functions are in seperate folder). These are run on the ethz landscape ecology group server.   
__model-functions__: Contains the set of functions to run given model (rf, glm, gam, brt).   
__evaluating-models__: Scripts to read and process the outputs of model calls. Here we want to produce one object that contains all the results and can be called into results scripts and figures scripts in standardized functions to produce a set of figures across the different ecological datasets.   
__figures-R2__: here a series of functions will act on the evaluating-models object to give a set of standard results of interest.   

## __data-processing__

### /functions

Function to extract absence records from a buffered size around a known occurrence record:  
_get_buffered_absences.R_  
  
Function to standardize rasters to a common grid scheme through inverse distance weighted interpolation:  
_standardize_raster.R_    

### /scripts  

Both of the RLS and BBS datasets have the following scripts: 

_data-processing_:  
• reads in full datasets and selects appropriate fields  
• aggregates abundances across years and local replicates  
• filters species by number of records  
• estimates species properties and filters to only 50 species in our different abundance classes  
• buffers absences for each species  
• saves a set of files containing species abundance information to which mdoels are later fitted and tests. 

_environmental-data_:  
• reads in and transforms all environmental rasters
• extracts environment values of sites from rasters
• performs PCA across relevant raster layers
• builds set of standarized rasters with a common grid system and transformations that match the covariate dataset
• produces the xy covariate dataset used in the models

_environmental-cross-validation_:  
• produces a set of data that have a in-the-bag sample to which models are fitted and a out-of-bag sample on which models are tests. 
• the out-the-bag is determined along one axis of a species environmental niche 
• produces a  set of oob cross-validation data to which models are later fitted and tested.



## __fitting-models__

all models are fitted on the server of the landscape ecology group at ethz. 
 
### /fit-models

_fit_all_models.R_:  
File that contains all model calls with varying parameter settings amongst functions. Loads all model functions and essentially runs all models. Is called by the configuration files. 

### /batchscripts-all  

Contains 4 .R scripts that produce 4 batch files which call the configuration files than run the fit_all_models.R script. 

### / rls or bbs 

Contains two scripts each, one which runs the fit_all_models call using in-the-bag cross validation, and one that calls the fit_all_models.R script but using out-the-bag cross validation. 



## __evaluating-models__

_00-bind-prediction-outputs-server.R_  
This step in the script binds together the prediction objects (observed vs. predictions for all the run models). This creates the input objects into the 01-evaluate-metrics-all-R2.R script. _bind-prediction-outputs.bat_ -> _00-bind-prediction-outputs-server.R_

_01-evaluate-metrics-all-R2.R_: 
For each species prediction and observations a set of evaluation metrics are calculated. Functions for this script are provided in _scripts/evaluating-models/functions/evaluation_functions.R_. Run with _batch_evaluate_metrics_all-R2_ -> _01-evaluate-metrics-all-R2.R_


## __figures__ 

_01-model-performance-figures.R_
Set of figures summarising the performance of different modelling frameworks. 
Functions for this script are provided in "scripts/figures/functions/model-performance-functions.R"

_02-model-predictions-figures.R_
Set of figures showing the observed and predicted values across different modelling frameworks. 
Functions for this script are provided in "scripts/figures/functions/model-prediction-functions.R"

_03-species-performance-figures.R_
to be written. aim is to provide figures summarising performance of models across species' attributes of abundance and occupancy. 

_04-spatial-projection-bbs-figures.R_
_04-spatial-projection-rls-figures.R_
to be writted. aim is to produce maps of difference in relative abundance and occurrence using 




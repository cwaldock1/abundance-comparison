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

## __organisation__: 

Scripts are organised into the following set of folders:  
  
__data-processing__: This set of scripts processes the raw survey data, subsets to our focal species and provides the inputs to the model scripts.   
__fitting-models__: This set of scripts provides a configuration script with calls a long script of model calls (model-functions are in seperate folder). These are run on the ethz landscape ecological group server.   
__model-functions__: Contains the set of functions to run focal models.   
__evaluating-models__: Scripts to read and process the outputs of model calls. Here we want to produce one object that contains all the results and can be called into results scripts and figures scripts in standardized functions to produce a set of figures across the different ecological datasets.   

Potentially necessary folders:   
I want to avoid producing different scripts for the BBS data and the RLS data, and instead have a set of functions that produce the results of interest across these inputs so that they are easily modifiable and I won't have to copy and past results and figures scripts across to sets of scripts each time I change something in one script.   
__results__: here a series of functions will act on the evaluating-models object to give a set of standard results of interest.   
__figures__: here a series of functions will act on the evaluating-models object to give a set of standard results of interest.  


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


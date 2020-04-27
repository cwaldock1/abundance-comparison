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

## organisation: 

Scripts are organised into the following set of folders:  
  
__data-processing__: This set of scripts processes the raw survey data, subsets to our focal species and provides the inputs to the model scripts.   
__fitting-models__: This set of scripts provides a configuration script with calls a long script of model calls (model-functions are in seperate folder). These are run on the ethz landscape ecological group server.   
__model-functions__: Contains the set of functions to run focal models.   
__evaluating-models__: Scripts to read and process the outputs of model calls. Here we want to produce one object that contains all the results and can be called into results scripts and figures scripts in standardized functions to produce a set of figures across the different ecological datasets.   

Potentially necessary folders:   
I want to avoid producing different scripts for the BBS data and the RLS data, and instead have a set of functions that produce the results of interest across these inputs so that they are easily modifiable and I won't have to copy and past results and figures scripts across to sets of scripts each time I change something in one script.   
__results__: here a series of functions will act on the evaluating-models object to give a set of standard results of interest.   
__figures__: here a series of functions will act on the evaluating-models object to give a set of standard results of interest  


## data-processing - RLS


scripts/01_rls-data-processing _aim of script is to produce a subset of the RLS data and species which provide an analysis of abundance. The output of this script is a matrix of sites by species that is sparse with all surveys performed in the autralian subset_

scripts/02_rls-environmental-data... 


## data-processing - BBS

### ../bbs/01-bbs-data-processing.R

Script to read in the raw BBS datasets and process following the same procedures as above. 


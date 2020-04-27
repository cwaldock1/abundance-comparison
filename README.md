# abundance-comparison

This README contains brief descriptions of the input data and the scripts used in the comparison of niche-abundance models. 

# ecological survey data

## Reef Life Survey

raw-data/RLS-spatial-fish-data-extract _provided by Rick Stuart-Smith on the 08.08.2019_

raw-data/RLS-site-subset-metadata _accessed from reeflifesurvey.org on 19.09.2019 for the purpose of providing and australian subset of data from which to assess preliminary models_

## Breeding Bird Survey

raw-data/breeding-bird-survey-usa/50-StopData/...
raw data downloaded via ftp server from the BBS website https://www.pwrc.usgs.gov/BBS/RawData/Choose-Method.cfm folder contains additional files descriptions and metadata.

## Environmental and spatial datasets

Documented in table XX in supporting materials.

# scripts

## data-processing - RLS


scripts/01_rls-data-processing _aim of script is to produce a subset of the RLS data and species which provide an analysis of abundance. The output of this script is a matrix of sites by species that is sparse with all surveys performed in the autralian subset_

scripts/02_rls-environmental-data... 


## data-processing - BBS

### ../bbs/01-bbs-data-processing.R

Script to read in the raw BBS datasets and process following the same procedures as above. 


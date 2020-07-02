# RUN  RANDOM FOREST MODEL FOR TWO STAGE APPROACH ----

# run occurrence models - two stage ---- 

source('scripts/model-functions/occ_rf.R')

occ_rf(abundance = abundance_input,
       validation = validation_input,
       covariates = covariates,
       species_name = unique(abundance_input$TAXONOMIC_NAME),
       n_bootstrap = n_boots,
       base_dir = base_dir,
       n.cores=1)

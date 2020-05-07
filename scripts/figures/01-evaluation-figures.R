# script to produce figures


# load in evaluation data ----

all_assessments <- lapply(list.files('results/model_assessment', full.names = T), readRDS)

plot_data <- do.call(rbind, all_assessments)

# full plots across all species and model scenarios for supporting materials 

plot_data %>% group_by(dataset, 
                       abundance_response, 
                       cross_validation)

all_model_plots(plot_data, )

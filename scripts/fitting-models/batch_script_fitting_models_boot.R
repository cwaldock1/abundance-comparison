# script to create batch_file to run RLS data models

setwd('/Volumes/Simulation/conor/abundance-comparison')

load("data/rls_abun_modelling_data_v2.RData")

length_i <- length(rls_abun_list)

Nodes <- c('@RunAsMultiple,@NodeSet_16')

file_out <- c()
for(i in 1:length_i){
  file_out[i] <- paste0('%_shared%R-latest%r.bat %conor%abundance-comparison%scripts%fitting-models%fitting_models_boot.R',
                        ' ', 
                        i)}

file_out <- gsub("%","\\", file_out, fixed=TRUE)

fileConn<-file("batch_fit_models_boot.bat")
writeLines(c(Nodes,file_out), fileConn)
close(fileConn)


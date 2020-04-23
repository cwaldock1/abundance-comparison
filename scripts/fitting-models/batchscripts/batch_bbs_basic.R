# script to create batch_file to run bbs data models

#setwd('/Volumes/Simulation/conor/abundance-comparison')

load("data/bbs_abun_modelling_data.RData")

length_i <- length(bbs_abun_list)

Nodes <- c('@RunAsMultiple,@NodeSet_16')

file_out <- c()
for(i in 1:length_i){
  file_out[i] <- paste0('%_shared%R-latest%r.bat %conor%abundance-comparison%scripts%fitting-models%bbs%bbs_config_basic.R',
                        ' ', 
                        i)}

file_out <- gsub("%","\\", file_out, fixed=TRUE)

fileConn<-file("scripts/fitting-models/batchscripts/batch_bbs_basic.bat")
writeLines(c(Nodes,file_out), fileConn)
close(fileConn)


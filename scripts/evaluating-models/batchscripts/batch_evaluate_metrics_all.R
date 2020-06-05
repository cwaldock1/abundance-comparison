# script to create batch_file to run bbs data models

#setwd('/Volumes/Simulation/conor/abundance-comparison')

Nodes <- c('@Node21')

file_out <- paste0('%_shared%R-latest%r.bat %conor%abundance-comparison%scripts%evaluating-models%01-evaluate-metrics-all.R')
file_out <- gsub("%","\\", file_out, fixed=TRUE)

fileConn<-file("scripts/evaluating-models/batchscripts/batch_evaluate_metrics_all.bat")
writeLines(c(Nodes,file_out), fileConn)
close(fileConn)


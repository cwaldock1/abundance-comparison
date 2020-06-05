# script to create batch_file to run bbs data models

#setwd('/Volumes/Simulation/conor/abundance-comparison')

Nodes <- c('@RunAsMultiple, @Node_NODE18')

spatial_scale <- c(0.1, 1, 5, 10, 20, 35, 50)
bind_files <- list.files('results/predictions_all/bind', full.names = T)



file_out <- lapply(seq_along(spatial_scale), function(x){
  file_out <- c()
  for(i in 1:length(bind_files)){
    file_out[i] <- paste0('%_shared%R-latest%r.bat %conor%abundance-comparison%scripts%evaluating-models%02-scale-evaluate-metrics.R', ' ', x, ' ', i)
  }
  return(file_out)
  })

file_out <- do.call(c, file_out)
file_out <- gsub("%","\\", file_out, fixed=TRUE)

fileConn<-file("scripts/evaluating-models/batchscripts/batch_scale_evaluate_metrics.bat")
writeLines(c(Nodes,file_out), fileConn)
close(fileConn)


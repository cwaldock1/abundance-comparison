# script to create batch_file to run RLS data models

rls_abun_list <- list.files('data/rls_all_basic')

length_i <- length(rls_abun_list)

Nodes <- c('@RunAsMultiple, @NodeSet_16, @noLog')

file_out <- c()
for(i in 1:length_i){
  file_out[i] <- paste0('%_shared%R-latest%r.bat %conor%abundance-comparison%scripts%fitting-rfocc-models%rls-all%rls_config_occ_rf.R',
                        ' ', 
                        i)}

file_out <- gsub("%","\\", file_out, fixed=TRUE)

fileConn<-file("scripts/fitting-rfocc-models/batchscripts-all/batch_rls_occ_rf.bat")
writeLines(c(Nodes,file_out), fileConn)
close(fileConn)


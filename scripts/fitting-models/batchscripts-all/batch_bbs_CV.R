# script to create batch_file to run RLS data models

bbs_abun_cv_list <- list.files('data/bbs_all_CV')

length_i <- length(bbs_abun_cv_list)

Nodes <- c('@RunAsMultiple, @NodeSet_16, @noLog')

file_out <- c()
for(i in 1:length_i){
  file_out[i] <- paste0('%_shared%R-latest%r.bat %conor%abundance-comparison%scripts%fitting-models%bbs-all%bbs_config_CV.R',
                        ' ', 
                        i)}

file_out <- gsub("%","\\", file_out, fixed=TRUE)

fileConn<-file("scripts/fitting-models/batchscripts-all/batch_bbs_CV.bat")
writeLines(c(Nodes,file_out), fileConn)
close(fileConn)



# load reqired libararies -------------------------------------------------

library(parallel)
library(stringr)
# One-hot encoding --------------------------------------------------------

AA <- c("A","C","T","G");
pat_list <- lapply(1:4, function(n){
  v <- rep(0,4);
  v[n]<-1;
  s <- v;
  return(s);
})
names(pat_list) <- AA;


# reading ASA and dynaseq static feature files ----------------------------

fltr_static <- read.table("/storage/obj2-29March2019/RCSB/dynaseq_feature_static.tsv")[,1:12]

fltr_ASA <- read.table("/storage/obj2-29March2019/RCSB/all_dna.rsa")[,c(1,2,3,5,7,9,11,13)]


# Prepare ASA for DNA pairs -----------------------------------------------

dna_name <- unique(sort(fltr_ASA[,1]))

system("rm -rf /storage/obj2-29March2019/RCSB/pat_pair")
system("mkdir -p /storage/obj2-29March2019/RCSB/pat_pair")

lapply(1:length(dna_name), function(n){
  mat <- fltr_ASA[fltr_ASA$V1 == dna_name[n],]
  mat_A <- mat[mat[,3] == "A",]
  mat_B <- mat[mat[,3] == "B",]
  
  mat_start <- cbind(mat_A[,c(1,2)], mat_B[seq(nrow(mat_B),1),2])
  mat_end <- mat_A[,c(4,5,6,7,8)] + mat_B[seq(nrow(mat_B),1),c(4,5,6,7,8)]
  
  chainA <- sapply(mat_start[,2], function(x) substr(x, 2,2))
  chainB <- sapply(mat_start[,3], function(x) substr(x, 2,2))
  
  CA <- paste("pat_list[[chainA[",seq(1,nrow(mat_A)),"]]]", sep = "")
  C_A <- lapply(1:nrow(mat_A),function(n3){eval(parse(text=CA[n3]))})
  patA <- t(do.call(cbind, C_A))
  
  CB <- paste("pat_list[[chainB[",seq(1,nrow(mat_B)),"]]]", sep = "")
  C_B <- lapply(1:nrow(mat_B),function(n3){eval(parse(text=CB[n3]))})
  patB <- t(do.call(cbind, C_B))
  
  
  my_dynamicA <- mydyna(paste(chainA, collapse = ""), dna_name[n])
  my_staticA <- dyna2static(my_dynamicA)[,1:12]
  
  my_dynamicB <- mydyna(paste(chainB, collapse = ""), dna_name[n])
  my_staticB <- dyna2static(my_dynamicB)[,1:12]
 
  my_static <- cbind(my_staticA, my_staticB)
  Mat <- rbind(NA,NA,my_static,NA,NA)
  
  myMat <- cbind(mat_start, patA, patB, Mat, mat_end)
  
  chName1 <- paste(colnames(my_staticA),"-1", sep = "")
  chName2 <- paste(colnames(my_staticB),"-2", sep = "")
  
  colnames(myMat) <- c("name", "chain_A", "chain_B", paste("pat-",seq(1,8), sep = ""), chName1, chName2, "All-atoms", "Total-Side", "Main-Chain", "Non-polar", "All polar")
  
  
  write.table(myMat , file = paste("/storage/obj2-29March2019/RCSB/pat_pair/", dna_name[n],".tsv", sep = ""), sep = "\t", row.names = F)
  
})
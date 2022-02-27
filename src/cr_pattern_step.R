
# load reqired libararies -------------------------------------------------

library(parallel)
library(stringr)
# One-hot encoding --------------------------------------------------------

AA <- c("AA/TT","AT/AT","AG/CT","AC/GT","CC/GG","CG/CG","CA/TG","CT/AG","GA/TC","GT/AC");
pat_list <- lapply(1:10, function(n){
  v <- rep(0,10);
  v[n]<-1;
  s <- v;
  return(s);
})
names(pat_list) <- AA;
pat_list[["TA/TA"]] <- pat_list$`AT/AT`
pat_list[["TT/AA"]] <- pat_list$`AA/TT`
pat_list[["TG/CA"]] <- pat_list$`AC/GT`
pat_list[["TC/GA"]] <- pat_list$`AG/CT`
pat_list[["GG/CC"]] <- pat_list$`CC/GG`
pat_list[["GC/GC"]] <- pat_list$`CG/CG`


# reading ASA and dynaseq static feature files ----------------------------

fltr_static <- read.table("/storage/obj2-29March2019/RCSB/dynaseq_feature_static.tsv")[,1:12]

fltr_ASA <- read.table("/storage/obj2-29March2019/RCSB/all_dna.rsa")[,c(1,2,3,5,7,9,11,13)]


# Prepare ASA for DNA pairs -----------------------------------------------

dna_name <- unique(sort(fltr_ASA[,1]))

system("rm -rf /storage/obj2-29March2019/RCSB/pat_step")
system("mkdir -p /storage/obj2-29March2019/RCSB/pat_step")

lapply(1:length(dna_name), function(n){
  mat <- fltr_ASA[fltr_ASA$V1 == dna_name[n],]
  mat_A <- mat[mat[,3] == "A",]
  mat_B <- mat[mat[,3] == "B",]
  
  mat_start <- cbind(mat_A[,c(1,2)], mat_B[seq(nrow(mat_B),1),2])
  mat_end <- mat_A[,c(4,5,6,7,8)] + mat_B[seq(nrow(mat_B),1),c(4,5,6,7,8)]
  
  chainA <- sapply(mat_start[,2], function(x) substr(x, 2,2))
  chainB <- sapply(mat_start[,3], function(x) substr(x, 2,2))
  
  my_dynamicA <- mydyna(paste(chainA, collapse = ""), dna_name[n])
  my_staticA <- dyna2static(my_dynamicA)[,1:12]
  
  my_dynamicB <- mydyna(paste(chainB, collapse = ""), dna_name[n])
  my_staticB <- dyna2static(my_dynamicB)[,1:12]
  
  my_static <- cbind(my_staticA, my_staticB)
  Mat <- rbind(NA,NA,my_static,NA,NA)
  
  enkod <- all_ftr <- all_trgt <- patt <- c()
  
  lapply(1:(length(chainA)-1), function(m){
    pat <- paste(paste(chainA[c(m,(m+1))], collapse = ""), "/", paste(chainB[c((m+1),m)], collapse = ""), sep = "")
    patt <<- c(patt,pat)
    enkod <<- rbind(enkod, pat_list[[pat]])
    
    ftr <- (Mat[m,] + Mat[(m+1),])/2
    all_ftr <<- rbind(all_ftr, ftr)
    
    trgt <- mat_end[m,] + mat_end[(m+1),]
    all_trgt <<- rbind(all_trgt, trgt)
  })
  
 mat_stp <- cbind(patt, enkod, all_ftr, all_trgt)
 
 chName1 <- paste(colnames(my_staticA),"-1", sep = "")
 chName2 <- paste(colnames(my_staticB),"-2", sep = "")
 
 colnames(mat_stp) <- c("base_step", paste("pat-",seq(1,10), sep = ""), chName1, chName2, "All-atoms", "Total-Side", "Main-Chain", "Non-polar", "All polar")
 
 write.table(mat_stp, file = paste("/storage/obj2-29March2019/RCSB/pat_step/", dna_name[n],".tsv", sep = ""), sep = "\t", row.names = F)
  
})
  
  
  
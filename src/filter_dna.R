

library(seqinr)



# To filter DNA sequences from RCBS-PDB database --------------------------


fasta_file <- dir("/storage/obj2-29March2019/RCSB/FASTA")

lapply(1:length(fasta_file), function(n){
  fasta <- read.fasta(paste("/storage/obj2-29March2019/RCSB/FASTA/",fasta_file[n], sep = ""))
  num_chain <- length(fasta)
  system("mkdir -p /storage/obj2-29March2019/RCSB/fltr_fasta")
  system("mkdir -p /storage/obj2-29March2019/RCSB/fltr_pdb")
  
  if(num_chain == 2){
    chain_len1 <- length(fasta[[1]])
    chain_len2 <- length(fasta[[2]])
    
    
    if(chain_len1 > 10 && chain_len1 == chain_len2){
      if(any(c(fasta[[1]], fasta[[2]]) %in% c("x","n","X","N","I","i","U","u"))){
        
      }else if(all(rev(comp(fasta[[1]])) == fasta[[2]])) {
        print(paste("This DNA is selected: ",fasta_file[n], sep = ""))
        system(paste("cp /storage/obj2-29March2019/RCSB/FASTA/",fasta_file[n], " /storage/obj2-29March2019/RCSB/fltr_fasta", sep = "" ))
        pr_name <- paste(tolower(strsplit(fasta_file[n],"[.]")[[1]][1]),".pdb", sep = "")
        system(paste("cp /storage/obj2-29March2019/RCSB/PDB/",pr_name, " /storage/obj2-29March2019/RCSB/fltr_pdb/", sep = "" ))
      } else {
        print( "Yupee! This DNA sequence is removed.")
      }
      
    }else{
      print(paste("This DNA is not selected: ",fasta_file[n], sep = ""))
    }
  }
  # Sys.sleep(0.05)
  return(num_chain)
})



..# To filter DNA sequences from NDB database -------------------------------

fasta_file <- dir("/storage/obj2-29March2019/NDB/FASTA")

lapply(1:length(fasta_file), function(n){
  fasta <- read.fasta(paste("/storage/obj2-29March2019/NDB/FASTA/",fasta_file[n], sep = ""))
  num_chain <- length(fasta)
  system("mkdir -p /storage/obj2-29March2019/NDB/fltr_fasta")
  system("mkdir -p /storage/obj2-29March2019/NDB/fltr_pdb")
  
  if(num_chain == 2){
    chain_len1 <- length(fasta[[1]])
    chain_len2 <- length(fasta[[2]])
    if(chain_len1 > 10 & chain_len1 == chain_len2){
      if(any(c(fasta[[1]], fasta[[2]]) %in% c("x","n","X","N","I","i","U","u"))){
        
      }else{
        print(paste("This DNA is selected: ",fasta_file[n], sep = ""))
        system(paste("cp /storage/obj2-29March2019/NDB/FASTA/",fasta_file[n], " /storage/obj2-29March2019/NDB/fltr_fasta/", sep = "" ))
        pr_name <- paste(tolower(strsplit(fasta_file[n],"[.]")[[1]][1]),".pdb", sep = "")
        system(paste("cp /storage/obj2-29March2019/NDB/PDB/",pr_name, " /storage/obj2-29March2019/NDB/fltr_pdb/", sep = "" ))
      }
      
    }else{
      print(paste("This DNA is not selected: ",fasta_file[n], sep = ""))
    }
  }
  # Sys.sleep(0.05)
  return(num_chain)
})


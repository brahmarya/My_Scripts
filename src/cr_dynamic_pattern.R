

# Remove the folders  -----------------------------------------------------

system("rm -rf /storage/obj2-29March2019/RCSB/dyna_pat")


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


# Reading ASA and dynaseq static feature files ----------------------------

fltr_dynamic <- read.table("E:/Academics/obj2-29March2019/RCSB/dynaseq_feature_dynamic.tsv")

fltr_ASA <- read.table("E:/Academics/obj2-29March2019/RCSB/all_dna.rsa")[,c(1,2,3,5,7,9,11,13)]

dna <-  sapply(rownames(fltr_dynamic), function(x) strsplit(x,"[].[]")[[1]][2])

# system("rm -rf E:/Academics/obj2-29March2019/RCSB/dyna_pat")

win_size <- c(1,3,5,7,9)

mclapply(1:5, function(w){
  win <- win_size[w];
  frnk <- w-1;
  vec <- c(-frnk:frnk);
  
  lapply(1:length(dna), function(n){
    naam <- strsplit(dna[n], ":")[[1]]
    mat <- fltr_ASA[fltr_ASA$V1 == tolower(naam[1]),]
    mat <- mat[mat[,3] == naam[2],]
    
    sat <- fltr_dynamic[str_detect(rownames(fltr_dynamic), dna[n]),]
    ptrn <- chain <- pattern <- dynafetr <- target <- c()
    
    lapply((frnk+1):(nrow(mat)-frnk), function(m){
      pat1 <- sapply(mat[(vec+m),2], function(x) substr(x, 2,2))
      b <- paste("pat_list[[pat1[",seq(1,win_size[w]),"]]]", sep = "")
      v1 <- unlist(lapply(1:win_size[w],function(n3){eval(parse(text=b[n3]))}))
      
      pat2 <- paste("[",dna[n],".pos.",vec+m,"]", sep = "")
      v2 <- as.vector(t(sat[pat2,]))
      
      names(v1) <- paste("pos", seq(1,(4*win_size[w])), sep = "")
      pattern <<- rbind(pattern,v1)
      
      names(v2) <- rep(colnames(fltr_dynamic),win_size[w])
      dynafetr <<- rbind(dynafetr,v2)
      
      ptrn <<- c(ptrn,paste(pat1,collapse =  ""))
      chain <<- c(chain, naam[2])
      
      target <<- rbind(target, mat[m,c(4,5,6,7,8)])
      
    })
    
    colnames(target) <- c("All-atoms", "Total-Side", "Main-Chain", "Non-polar", "All polar")
    newMat <- cbind(tolower(naam[1]),ptrn, chain, pattern, dynafetr, target)
    colnames(newMat)[1] <- "name"
    
    system(paste("mkdir -p E:/Academics/obj2-29March2019/RCSB/dyna_pat/win_",win, sep = ""))
    # write.table(newMat, file = paste("E:/Academics/obj2-29March2019/RCSB/pat/win_",win,"/",dna[n],".tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE )
    
    write.table(newMat, file = paste("E:/Academics/obj2-29March2019/RCSB/dyna_pat/win_",win,"/",naam[1],"_",naam[2],".tsv", sep = "", collapse = "_"), quote = FALSE, sep = "\t", row.names = FALSE )
  })
}, mc.cores = 1)

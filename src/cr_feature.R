
library(parallel)

data  <- read.table("/storage/ICMR/RCSB/all_dna.rsa")
dna <-  unique(sort(data[,1]))

AA <- c("A","C","T","G");
pat_list <- lapply(1:4, function(n){
  v <- rep(0,4);
  v[n]<-1;
  s <- v;
  return(s);
})
names(pat_list) <- AA;


load("/storage/ICMR/dictionary.static.Rdata")


dyna.name <- rownames(dynaseq.str)
names <- sapply(dyna.name, function(x) substr(x, 1,nchar(x)-1))

rownames(dynaseq.str) <- names
sort_name <- unique(sort(rownames(dynaseq.str)))
dynaseq <- dynaseq.str[sort_name,]


mclapply(1:length(dna), function(n){
  mat <- data[data$V1 == dna[n],]
  mat <- mat[mat[,3] == mat[1,3],]
  
  win <- 5
  frnk <- 2
  vec <- c(-frnk:frnk);
  ptrn <- chain <- pattern <- dynafetr <- target <- c()
  lapply(3:(nrow(mat)-2), function(m){
    pat1 <- sapply(mat[(vec+m),2], function(x) substr(x, 2,2))
    b <- paste("pat_list[[pat1[",seq(1,5),"]]]", sep = "")
    v1 <- unlist(lapply(1:5,function(n3){eval(parse(text=b[n3]))}))
    v2 <- dynaseq[paste(pat1,collapse =  ""),1:12]
    
    names(v1) <- paste("pos", seq(1,20), sep = "")
    pattern <<- rbind(pattern,v1)
    dynafetr <<- rbind(dynafetr,v2)
    
    ptrn <<- c(ptrn,paste(pat1,collapse =  ""))
    chain <<- c(chain, as.character(mat[m,3]))
    
    target <<- rbind(target, mat[m,c(5,7,9,11,13)])
  })
  newMat <- cbind(ptrn, chain, pattern, dynafetr, target)
  write.table(newMat, file = paste("/storage/ICMR/RCSB/pattern_files/",dna[n],".tsv", sep = ""), quote = FALSE, row.names = FALSE)
  return(0)
}, mc.cores = 16)





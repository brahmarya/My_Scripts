
library(ChIPseeker)
library(protr)
library(dynaseqR)
library(stringr)

load("/storage/obj2-29March2019/Dictionary/Static.win1.Rdata")
narrowpeak <- read.table("/storage/obj2-29March2019/data/human_genome.fa")

test <- narrowpeak[which(narrowpeak[,1] == "chr1"),]
all.max.asa <- all.mean.asa <- c()

system.time(all <- lapply(1:nrow(test), function(n) {
  if(str_detect(test[n,4], "N")){
    all.max.asa <<- rbind(all.max.asa, c(rep("NaN",5)))
    all.mean.asa <<- rbind(all.mean.asa, c(rep("NaN",5)))
  } else {
    frag <- toupper(substring(as.character(test[n,4]), 1:(str_length(test[n,4])-5+1),5:str_length(test[n,4])))
    
    all.max.asa <- apply(static.asa[frag,], 2, max)
    all.mean.asa <- apply(static.asa[frag,], 2, mean)
    
    # all.max.asa <<- rbind(all.max.asa, max.asa)
    # all.mean.asa <<- rbind(all.mean.asa, mean.asa)
  }
  
  return(c(all.max.asa,all.mean.asa))
})
)


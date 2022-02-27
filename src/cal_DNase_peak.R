

# To retrieve peak scores from all the files for selected coordinates --------



data <- read.table("/storage/obj2-29March2019/data/All_cell_intersect_narrow.bed")
head(data)
median <- apply(data[,c(2,3)],1,mean)

files <- list.files("/storage/obj2-29March2019/DNase Data", pattern = "Rep1.narrowPeak.gz", full.names = TRUE)

peak.score <- matrix(NA, nrow = length(median), ncol = 13)
lapply(1:length(files), function(n){
  data <- read.table(gzfile(files[n]))
  lapply(1:nrow(data), function(m){
    if(is.element("TRUE" , c(median %in% seq(data[m,2], data[m,3])))){
      peak.score[which(median %in% seq(data[m,2], data[m,3])),n+1] <<- data[m,7]
    }
  })
})

peak.score[,1] <- median
peak.score[,12] <- apply(peak.score[,2:11],1,mean, na.rm = TRUE)
peak.score[,13] <- apply(peak.score[,2:11],1,max, na.rm = TRUE)

colnames(peak.score) <- c("Median", "A549","Gm04503", "H7es", "Helas3", "Hepg2", "Hre", "Jurkat", "K562", "Mcf7", "Sknmc", "Mean", "Max" )

write.table(peak.score, file = "/storage/obj2-29March2019/data/All_cell_intersect_narrow_peak.tsv", sep = "\t", col.names = TRUE)

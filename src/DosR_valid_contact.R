library(dplyr)
library(BioCircos)


# Calculate distance frequency for Ra -------------------------------------

# cont <- read.table("/storage/Collab/Mohit/Contact map/h37rv_25kbdense.txt")
cont <- read.table("/storage/Collab/Mohit/R_Matrices_HiC/Matrices_Normalised/Ra_5kb_norm.tsv")
rownames(cont) <- paste(seq(0,4415,5),"Kb", sep = "")
colnames(cont) <- paste(seq(0,4415,5),"Kb", sep = "")
# napos <- unique(which(is.na(cont), arr.ind=TRUE)[,1])
# s <- apply(cont,1,sd, na.rm = T)
# napos <- which(s == 0, arr.ind=TRUE)
# cont <- cont[-napos, -napos]
# res <- chisq.test(cont)
# 

cont_mat <- cbind(expand.grid(Start=seq(0,length.out = 884, by = 5),End=seq(0,length.out = 884, by = 5),stringsAsFactors = FALSE),Diff=NA, Dist=NA, Frequency=NA)
cont_mat[,3] <- Mod(cont_mat[,1]-cont_mat[,2])
i<-1;
lapply(1:nrow(cont), function(n){
  lapply(1:nrow(cont), function(m){
    if(cont_mat[i,3] > 2210){
      cont_mat[i,4] <<- (2210 - (cont_mat[i,3] - 2210))
    } else {
      cont_mat[i,4] <<- cont_mat[i,3]
    }
    cont_mat[i,5] <<- cont[m,n];
    i<<-i+1;
    return(0)})
  return(0)})


dist <- unique(sort(cont_mat$Dist))

Out_Ra <- lapply(1:length(dist), function(n){
  res <- mean(cont_mat$Frequency[which(cont_mat$Diff == dist[n])], na.rm = TRUE)
  return(res)
})


names(Out_Ra) <- dist

# Calculate distance frequency for Rv -------------------------------------

cont <- read.table("/storage/Collab/Mohit/R_Matrices_HiC/Matrices_Normalised/Rv_5kb_norm.tsv")
rownames(cont) <- paste(seq(0,4410,5),"Kb", sep = "")
colnames(cont) <- paste(seq(0,4410,5),"Kb", sep = "")

cont_mat <- cbind(expand.grid(Start=seq(0,length.out = 883, by = 5),End=seq(0,length.out = 883, by = 5),stringsAsFactors = FALSE),Diff=NA, Dist=NA, Frequency=NA)
cont_mat[,3] <- Mod(cont_mat[,1]-cont_mat[,2])
i<-1;
lapply(1:nrow(cont), function(n){
   lapply(1:nrow(cont), function(m){
      if(cont_mat[i,3] > 2210){
         cont_mat[i,4] <<- (2210 - (cont_mat[i,3] - 2210))
      } else {
         cont_mat[i,4] <<- cont_mat[i,3]
      }
      cont_mat[i,5] <<- cont[m,n];
      i<<-i+1;
      return(0)})
   return(0)})


dist <- unique(sort(cont_mat$Dist))

Out_Rv <- lapply(1:length(dist), function(n){
   res <- mean(cont_mat$Frequency[which(cont_mat$Diff == dist[n])], na.rm = TRUE)
   return(res)
})


names(Out_Rv) <- dist



# DosR Matrix for Ra -------------------------------------------------------------

cont <- read.table("/storage/Collab/Mohit/R_Matrices_HiC/Matrices_Normalised/Ra_5kb_norm.tsv")
rownames(cont) <- seq(0,4415,5)
colnames(cont) <- seq(0,4415,5)

DosR <- read.table("/storage/Collab/Mohit/DosR/DosR Regulon gene coordinates new.tsv")[4]
DosR <- unique(DosR[,1])

DosR.Mat <- cont[as.character(DosR),as.character(DosR)]


cont_mat <- cbind(expand.grid(DosR,DosR),Diff=NA, Dist=NA, AvgFrequency=NA, Frequency=NA, Freq_diff=NA)
cont_mat[,3] <- Mod(cont_mat[,1]-cont_mat[,2])

i<-1;
lapply(1:nrow(DosR.Mat), function(n){
  lapply(1:nrow(DosR.Mat), function(m){
    if(cont_mat[i,3] > 2210){
      cont_mat[i,4] <<- (2210 - (cont_mat[i,3] - 2210))
    } else {
      cont_mat[i,4] <<- cont_mat[i,3]
    }
    cont_mat[i,6] <<- DosR.Mat[m,n];
    cont_mat[i,5] <<- Out_Ra[[as.character(cont_mat[i,4])]]
    cont_mat[i,7] <<- cont_mat[i,6] - cont_mat[i,5]
    i<<-i+1;
    return(0)})
  return(0)})

write.csv(cont_mat, file = "/storage/Collab/Mohit/DosR/Ra_DosR Regulon contacts.csv")


cont_mat <- cont_mat[which(!duplicated(apply(cont_mat,1,sum))), ]
cont_mat <- na.omit(cont_mat)

top_Ra <- cont_mat %>% filter(Dist > 25 & Freq_diff > median(cont_mat$Freq_diff))


# DosR Matrix for Rv -------------------------------------------------------------

cont <- read.table("/storage/Collab/Mohit/R_Matrices_HiC/Matrices_Normalised/Rv_5kb_norm.tsv")
rownames(cont) <- seq(0,4410,5)
colnames(cont) <- seq(0,4410,5)

DosR <- read.table("/storage/Collab/Mohit/DosR/DosR Regulon gene coordinates new.tsv")[4]
DosR <- unique(DosR[,1])

DosR.Mat <- cont[as.character(DosR),as.character(DosR)]


cont_mat <- cbind(expand.grid(DosR,DosR),Diff=NA, Dist=NA, AvgFrequency=NA, Frequency=NA, Freq_diff=NA)
cont_mat[,3] <- Mod(cont_mat[,1]-cont_mat[,2])

i<-1;
lapply(1:nrow(DosR.Mat), function(n){
   lapply(1:nrow(DosR.Mat), function(m){
      if(cont_mat[i,3] > 2210){
         cont_mat[i,4] <<- (2210 - (cont_mat[i,3] - 2210))
      } else {
         cont_mat[i,4] <<- cont_mat[i,3]
      }
      cont_mat[i,6] <<- DosR.Mat[m,n];
      cont_mat[i,5] <<- Out_Rv[[as.character(cont_mat[i,4])]]
      cont_mat[i,7] <<- cont_mat[i,6] - cont_mat[i,5]
      i<<-i+1;
      return(0)})
   return(0)})

write.csv(cont_mat, file = "/storage/Collab/Mohit/DosR/Rv_DosR Regulon contacts.csv")


cont_mat <- cont_mat[which(!duplicated(apply(cont_mat,1,sum))), ]
cont_mat <- na.omit(cont_mat)

top_Rv <- cont_mat %>% filter(Dist > 25 & Freq_diff > median(cont_mat$Freq_diff))



# Circos Plot -------------------------------------------------------------


tb_size <- list("Tuberclosis" = 4425)

# BioCircos(genome = tb_size, genomeFillColor = "darkblue", chrPad = 0, displayGenomeBorder = FALSE, genomeTicksDisplay = FALSE, genomeLabelDy = 0)

links_chr_1 = rep("Tuberclosis", nrow(top))
links_chr_2 = rep("Tuberclosis", nrow(top))

links_pos_1 = top_Ra[,1]
links_pos_2 = top_Ra[,2]
# links_labels = paste("Link ", seq(1,nrow(high_cont)), sep = "")

tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 0.95, borderSize = 0, fillColors = "#EEFFEE")  

# tracklist = tracklist + BioCircosLinkTrack('myLinkTrack', links_chr_1, links_pos_1, links_pos_1, links_chr_2, links_pos_2, links_pos_2 , maxRadius = 0.95, displayLabel = F, width = 0.3 )

links_chr_3 = rep("Tuberclosis", nrow(high_cont))
links_chr_4 = rep("Tuberclosis", nrow(high_cont))

links_pos_3 = top_Rv[,1]
links_pos_4 = top_Rv[,2]
# links_labels = paste("Link ", seq(1,nrow(high_cont)), sep = "")

tracklist = tracklist + BioCircosLinkTrack("testLink", links_chr_1, links_pos_1, links_pos_1, links_chr_2, links_pos_2, links_pos_2, axisPadding = 6, color = "#0000FF", width = 0.5, displayLabel = F, maxRadius = 0.95)

tracklist = tracklist + BioCircosLinkTrack("testLink2", links_chr_3, links_pos_3, links_pos_3, links_chr_4, links_pos_4, links_pos_4, axisPadding = 6, displayLabel = F, color = "#FF6666",width = 0.7, maxRadius = 0.95)



BioCircos(genome= tb_size,tracklist, genomeFillColor = "Spectral", chrPad = 0, displayGenomeBorder = FALSE, genomeTicksScale = 100, genomeLabelTextSize = "8pt", genomeLabelDy = 15)



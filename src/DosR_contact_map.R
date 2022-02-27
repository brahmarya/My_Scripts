library(dplyr)
library(BioCircos)
# To generate contact map of DosR from H37Rv and H37Ra contact map --------


# data <- read.table("/storage/Collab/Mohit/Contact map/h37rv_25kbdense.txt") # Old Contact Matrix
data <- read.table("/storage/Collab/Mohit/R_Matrices_HiC/Matrices_Normalised/Ra_5kb_norm.tsv")
rownames(data) <- seq(0,4415,5)
colnames(data) <- seq(0,4415,5)
# DosR <- read.table("/storage/Collab/Mohit/DosR/DosR Regulon gene coordinates.tsv")[5]
DosR <- read.table("/storage/Collab/Mohit/DosR/DosR Regulon gene coordinates new.tsv")[4]
DosR <- unique(sort(DosR[,1]))

DosR.Mat <- data[as.character(DosR),as.character(DosR)]

write.table(DosR.Mat, file = "/storage/Collab/Mohit/Contact map/DosR_new.tsv", sep = "\t")

library(BioCircos)


data <- DosR.Mat
cont_mat <- cbind(expand.grid(as.numeric(rownames(DosR.Mat)),as.numeric(rownames(DosR.Mat))),0)
i<-1;
lapply(1:nrow(data), function(n){
  lapply(1:nrow(data), function(m){
    cont_mat[i,3] <<- data[m,n];
    i<<-i+1;
    return(0)})
  return(0)})

# cont_3q <- cont_mat[which(cont_mat[,3] > 56),]
# d <- cont_3q[duplicated(cont_3q[,3]),]

cont_mat <- cont_mat[order(cont_mat[,3], decreasing = TRUE), ]
d <- cont_mat

high_cont <- d[1:10,]
mod_cont <- d[11:30,]

tb_size <- list("Tuberclosis" = 4425)

# BioCircos(genome = tb_size, genomeFillColor = "darkblue", chrPad = 0, displayGenomeBorder = FALSE, genomeTicksDisplay = FALSE, genomeLabelDy = 0)

links_chr_1 = rep("Tuberclosis", nrow(high_cont))
links_chr_2 = rep("Tuberclosis", nrow(high_cont))

links_pos_1 = high_cont[,1]
links_pos_2 = high_cont[,2]
# links_labels = paste("Link ", seq(1,nrow(high_cont)), sep = "")

tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 0.95, borderSize = 0, fillColors = "#EEFFEE")  

# tracklist = tracklist + BioCircosLinkTrack('myLinkTrack', links_chr_1, links_pos_1, links_pos_1, links_chr_2, links_pos_2, links_pos_2 , maxRadius = 0.95, displayLabel = F, width = 0.3 )

links_chr_3 = rep("Tuberclosis", nrow(high_cont))
links_chr_4 = rep("Tuberclosis", nrow(high_cont))

links_pos_3 = mod_cont[,1]
links_pos_4 = mod_cont[,2]
# links_labels = paste("Link ", seq(1,nrow(high_cont)), sep = "")

tracklist = tracklist + BioCircosLinkTrack("testLink", links_chr_1, links_pos_1, links_pos_1, links_chr_2, links_pos_2, links_pos_2, axisPadding = 6, color = "#0000FF", width = 0.5, displayLabel = F, maxRadius = 0.95)

tracklist = tracklist + BioCircosLinkTrack("testLink2", links_chr_3, links_pos_3, links_pos_3, links_chr_4, links_pos_4, links_pos_4, axisPadding = 6, displayLabel = F, color = "#FF6666",width = 0.7, maxRadius = 0.95)



BioCircos(genome= tb_size,tracklist, genomeFillColor = "Spectral", chrPad = 0, displayGenomeBorder = FALSE, genomeTicksScale = 100, genomeLabelTextSize = "8pt", genomeLabelDy = 15, zoom = TRUE)








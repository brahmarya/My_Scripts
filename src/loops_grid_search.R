


system("hicConvertFormat --matrices /storage/Collab/Mohit/10Kb/hicMatrix/replicateMergedCorrected_ICE_ra.h5 --inputFormat h5 -o /storage/Collab/Mohit/10Kb/hicMatrix/replicateMergedCorrected_ICE_ra --outputFormat ginteractions")
system("hicConvertFormat --matrices /storage/Collab/Mohit/10Kb/hicMatrix/replicateMergedCorrected_ICE_rv.h5 --inputFormat h5 -o /storage/Collab/Mohit/10Kb/hicMatrix/replicateMergedCorrected_ICE_rv --outputFormat ginteractions")

library(GenomicRanges)
library(InteractionSet)

ra <- read.delim("/storage/Collab/Mohit/10Kb/hicMatrix/replicateMergedCorrected_ICE_ra.tsv", header = FALSE)[7]
rv <- read.delim("/storage/Collab/Mohit/10Kb/hicMatrix/replicateMergedCorrected_ICE_rv.tsv", header = FALSE)[7]

pval <- c(0.05, 0.04, 0.03, 0.02, 0.01)
ra_pit <- quantile(ra, na.rm = TRUE, probs = seq(0.1,1,0.1))[c(6,7,8,9)]
rv_pit <- quantile(rv, na.rm = TRUE, probs = seq(0.1,1,0.1))[c(6,7,8,9)]

Mat <- c()

lapply(1:length(pval), function(n1){
   lapply(1:length(ra_pit), function(n2){
      cmd_ra <- paste("hicDetectLoops -m /storage/Collab/Mohit/10Kb/hicMatrix/replicateMergedCorrected_ICE_ra.h5 -o /storage/Collab/Mohit/10Kb/Loops/Ra_loops.bedgraph --maxLoopDistance 2000000 --windowSize 5 --peakWidth 3 --pValuePreselection ", pval[n1], " --pValue ", pval[n1], " --peakInteractionsThreshold ", ra_pit[n2], sep = "")
      
      cmd_rv <- paste("hicDetectLoops -m /storage/Collab/Mohit/10Kb/hicMatrix/replicateMergedCorrected_ICE_rv.h5 -o /storage/Collab/Mohit/10Kb/Loops/Rv_loops.bedgraph --maxLoopDistance 2000000 --windowSize 5 --peakWidth 3 --pValuePreselection ", pval[n1], " --pValue ", pval[n1], " --peakInteractionsThreshold ", rv_pit[n2], sep = "")
      
      system(cmd_ra)
      system(cmd_rv)
      
      ra_loop <- nrow(read.table("/storage/Collab/Mohit/10Kb/Loops/Ra_loops.bedgraph"))
      rv_loop <- nrow(read.table("/storage/Collab/Mohit/10Kb/Loops/Rv_loops.bedgraph"))
      
      mat <- c("Ra_pval" = pval[n1], "Ra_pit" = ra_pit[n2], "Rv_pval" = pval[n1], "Rv_pit" = rv_pit[n2], "Ra_loop" = ra_loop, "Rv_loop" = rv_loop)
      Mat <<- rbind(Mat, mat)
      
   })
})



write.table(Mat, file = "/storage/Collab/Mohit/10Kb/loops_grid_ICE.tsv", sep = "\t", quote = FALSE)



step <- c(10000,20000)
minBound <- c(10000,20000)
minDepth <- c(40000,50000)
maxDepth <- c(50000,60000,70000)
fdr <- c(0.01,0.04,0.05)

Mat <- c()
system("rm /storage/Collab/Mohit/HiC_New/Ra_hic/TAD/TADs-5kb_boundaries.bed")
system("rm /storage/Collab/Mohit/HiC_New/Rv_hic/TAD/TADs-5kb_boundaries.bed")

lapply(step, function(i){
   lapply(minBound, function(j){
      lapply(minDepth, function(k){
         lapply(maxDepth, function(l){
            lapply(fdr, function(m){
               Ra_tad <- paste("hicFindTADs --matrix /storage/Collab/Mohit/5Kb/hicMatrix/replicateMergedCorrected_ICE_ra.h5 --minDepth ", k, " --maxDepth ", l, " --numberOfProcessors 14 --step ", i, " --outPrefix /storage/Collab/Mohit/HiC_New/Ra_hic/TAD/TADs-5kb  --minBoundaryDistance ", j, " --correctForMultipleTesting fdr --threshold ", m)
               system(Ra_tad)
               if(file.info("/storage/Collab/Mohit/HiC_New/Ra_hic/TAD/TADs-5kb_boundaries.bed")$size > 0){
                  lapply(step, function(n){
                     lapply(minBound, function(o){
                        lapply(minDepth, function(p){
                           lapply(maxDepth, function(q){
                              lapply(fdr, function(r){
                                 Rv_tad <- paste("hicFindTADs --matrix /storage/Collab/Mohit/5Kb/hicMatrix/replicateMergedCorrected_ICE_rv.h5 --minDepth ", p, " --maxDepth ", q, " --numberOfProcessors 14 --step ", n, " --outPrefix /storage/Collab/Mohit/HiC_New/Rv_hic/TAD/TADs-5kb  --minBoundaryDistance ", o, " --correctForMultipleTesting fdr --threshold ", r)
                                 system(Rv_tad)
                                 
                                 if(file.info("/storage/Collab/Mohit/HiC_New/Rv_hic/TAD/TADs-5kb_boundaries.bed")$size > 0){
                                    if(k < l){
													if(p < q){
                                   Ra_score <- read.table("/storage/Collab/Mohit/HiC_New/Ra_hic/TAD/TADs-5kb_boundaries.bed")[5]
                                   Rv_score <- read.table("/storage/Collab/Mohit/HiC_New/Rv_hic/TAD/TADs-5kb_boundaries.bed")[5]
                                   p_val <- t.test(Ra_score, Rv_score)[3]
											  Ra_tad <- nrow(Ra_score)
											  Rv_tad <- nrow(Rv_score)

                                   
                                   mat <- c("Ra_step" = i, "Ra_minBoundary" = j, "Ra_minDepth" = k, "Ra_maxDepth" = l, "Ra_fdr" = m, "Rv_step" = n, "Rv_minBoundary" = o, "Rv_minDepth" = p, "Rv_maxDepth" = q, "Rv_fdr" = r, p_val[1], "Ra_TAD" = Ra_tad, "Rv_TAD" = Rv_tad)
                                   Mat <<- rbind(Mat, mat)
                                   show(mat) 
                                 }
										}
									}
                                 
                                 
                              })
                           })
                        })
                     })
                  })
                  
               }

            })
         })
      })
   })
})


write.table(Mat, file = "/storage/Collab/Mohit/Results/TAD_score_matrix_ICE_5.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

mydyna <- function (sequences = sequences, names = names) 
{
  library("e1071")
  path <- find.package("dynaseqR")
  load(paste(path, "/data/dynaseq-ws5.svm", sep = ""))
  load(paste(path, "/data/encoder.Rdata", sep = ""))
  if (is.null(sequences)) {
    sequences <- read.table("sequences.txt")
    sequences <- as.character(sequences[, 2])
  }
  sequences <- as.matrix(sequences)
  sqnames <- paste("sequence.", c(1:nrow(sequences)), sep = "")
  rownames(sequences) <- sqnames
  nbrsize = 2
  encoding <- c(0, 0, 0, 1)
  encoding <- rbind(encoding, c(0, 0, 1, 0))
  encoding <- rbind(encoding, c(0, 1, 0, 0))
  encoding <- rbind(encoding, c(1, 0, 0, 0))
  rownames(encoding) <- c("A", "C", "G", "T")
  res_all <- c()
  myrownames <- c()
  names <- as.matrix(names)
  
  for (seqid in c(1:nrow(sequences))) {
    allinputs <- c()
    mynames <- c()
    sequence <- sequences[seqid, 1]
    for (winpos in c(3:(nchar(sequence) - 2))) {
      mypat <- c()
      mynames <- c(mynames, paste("[", names[seqid], 
                                  ".pos.", winpos, "]", sep = ""))
      for (nbr in seq(-1 * nbrsize, nbrsize)) {
        tmppat <- substr(sequence, winpos + nbr, winpos + 
                           nbr)
        tmppat <- encoding[tmppat, ]
        mypat <- c(mypat, tmppat)
      }
      allinputs <- rbind(allinputs, mypat)
    }
    params <- rownames(dynaseq[[61]])
    pheader = c()
    for (i in c(1:nrow(dynaseq[[61]]))) {
      for (j in c(1:(ncol(dynaseq[[61]]) - 1))) {
        pheader <- c(pheader, paste(params[i], "-bin-", 
                                    j, sep = ""))
      }
    }
    res <- c()
    for (binid in c(1:60)) {
      pred <- predict(dynaseq[[binid]], allinputs)
      res <- rbind(res, signif(pred, 3))
    }
    rownames(res) <- pheader
    colnames(res) <- mynames
    res <- t(res)
    res_all <- rbind(res_all, res)
    myrownames <- c(myrownames, mynames)
  }
  rownames(res_all) <- myrownames
  return(res_all)
}
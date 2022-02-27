

# load required libraries  ------------------------------------------------

library(protr)
library(dynaseqR)



# read filtered fasta file ------------------------------------------------

fasta_data <- readFASTA("/storage/obj2-29March2019/RCSB/fltr_fasta.fasta")

chain <- names(fasta_data)
ch <- sapply(chain, function(x) strsplit(x,"[|]")[[1]][1])

# calculate dynaseq features for each sequence ----------------------------

ftr_dynamic <- mydyna(sequences = fasta_data, names = ch)

write.table(ftr_dynamic, file = "/storage/obj2-29March2019/RCSB/dynaseq_feature_dynamic.tsv", quote = FALSE)


ftr_static <- dyna2static(ftr_dynamic)

write.table(ftr_static, file = "/storage/obj2-29March2019/RCSB/dynaseq_feature_static.tsv", quote = FALSE)

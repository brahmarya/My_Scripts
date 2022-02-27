

data <- read.table("/storage/Collab/Mohit/data/Operon/H37Ra_operon.txt", header = TRUE)
h <- hist(data$Median, breaks = 45, plot = FALSE)
density <- h$density;
count <- h$counts

start <- seq(0,4419977, 100000)
end <- c(seq(100000, 4419977, 100000),4419977)

operon_den <- as.data.frame(cbind("NC_009525.1", start, end, density, count))

write.table(operon_den, file = "/storage/Collab/Mohit/data/Operon/H37Ra_operon_den.bed", col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")


data <- read.table("/storage/Collab/Mohit/data/Operon/H37Rv_operon.txt", header = TRUE)
h <- hist(data$Median, breaks = 45, plot = FALSE)
density <- h$density;
count <- h$counts
 
start <- seq(0, 4411532, 100000)
end <- c(seq(100000, 4411532, 100000),4411532)

operon_den <- as.data.frame(cbind("NC_000962.3", start, end, density, count))

write.table(operon_den, file = "/storage/Collab/Mohit/data/Operon/H37Rv_operon_den.bed", col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")

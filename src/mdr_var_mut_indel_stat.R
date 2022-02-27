
path <- dir("/storage/Collab/Mohit/MDR_Data", full.names = TRUE)
file <- sapply(path,function(x) strsplit(x,"[/-]")[[1]][6])
fullPath <- dir("/storage/Collab/Mohit/MDR_Data", full.names = TRUE)

mut_mat <- c()
indel_mat <- c()
M_Data <- I_Data <- c()

lapply(1:28, function(n){
   
   Path <- paste(fullPath[n],"/Rv_freebayes.vcf.gz", sep = "")
   data <- read.table(Path, stringsAsFactors = FALSE)
   mut_data <- data[-c(which(nchar(data[,4]) > 1)), ]
   indel_data <- data[c(which(nchar(data[,4]) > 1)), ]
   
   mut <- c("A/C","A/G","A/T","C/A","C/G","C/T","G/A","G/C","G/T","T/A","T/C","T/G")
   
   t <- table(paste(mut_data[,4], mut_data[,5], sep ="/"))
   mut_mat <<- rbind(mut_mat, t[mut])
   
   indel_mat <<- rbind(indel_mat, c(nrow(mut_data),nrow(indel_data)))
   
   
   M_Data <<- rbind(M_Data, mut_data)
   I_Data <<- rbind(I_Data, indel_data)
  
   return(0) 
})
   
rownames(mut_mat) <- file

rownames(indel_mat) <- file
colnames(indel_mat) <- c("Mutations", "Indels")
   
write.table(mut_mat, file = "/storage/Collab/Mohit/Results/Point_mutations_MDR.tsv", sep = "\t", quote = FALSE)
write.table(indel_mat, file = "/storage/Collab/Mohit/Results/Indel_mutations_MDR.tsv", sep = "\t", quote = FALSE)



# Plots -------------------------------------------------------------------

barplot(mut_mat, col = c("red", "blue"), xlab = "Mutations", ylab = "Mutation Frequency", ylim = c(0, 700))
par(xpd=TRUE)
legend("topright", inset=c(-0.01,-0.04), c("NE24", "NE30"), fill = c("red", "blue"))

barplot(t(indel_mat), beside = T, col = c("red","blue"), xlab = "MDR Variant", ylab = "Indel Frequncy", cex.names = 0.67, ylim = c(0, 3500))


gap.barplot(t(indel_mat), col = rep(c("red","blue"), 28), xlab = "Clinical Variant", ylab = "Indel Frequency", gap=c(3500, 12000), ylim = c(0,12000), ytics = c(0,500,1000,1500,2000,2500,15000,16000,17000,18000,19000), xaxt='n')
axis.break(2, 3500, breakcol="snow", style="gap")
axis.break(2, 3500*(1+0.02), breakcol="black", style="slash")
axis.break(4, 3500*(1+0.02), breakcol="black", style="slash")
axis(2, at=3500)
axis(1, at=seq(1.5,56, 2), file , las = 2, cex.axis = 0.67)







data1 <- read.table("E:/Academics/obj2-29March2019/Results/Seq/win_3_err.tsv")
data2 <- read.table("E:/Academics/obj2-29March2019/Results/Seq/win_5_err.tsv")
data3 <- read.table("E:/Academics/obj2-29March2019/Results/Seq/win_7_err.tsv")
data4 <- read.table("E:/Academics/obj2-29March2019/Results/Seq/win_9_err.tsv")

seq_mae <- cbind(NA,data1[,2], data2[,2], data3[,2], data4[,2])
seq_cor <- cbind(NA,data1[,5], data2[,5], data3[,5], data4[,5])

data1 <- read.table("E:/Academics/obj2-29March2019/Results/DynaSeq/win_1_err.tsv")
data2 <- read.table("E:/Academics/obj2-29March2019/Results/DynaSeq/win_3_err.tsv")
data3 <- read.table("E:/Academics/obj2-29March2019/Results/DynaSeq/win_5_err.tsv")
data4 <- read.table("E:/Academics/obj2-29March2019/Results/DynaSeq/win_7_err.tsv")
data5 <- read.table("E:/Academics/obj2-29March2019/Results/DynaSeq/win_9_err.tsv")

dynaseq_mae <- cbind(data1[,2], data2[,2], data3[,2], data4[,2], data5[,2])
dynaseq_cor <- cbind(data1[,5], data2[,5], data3[,5], data4[,5], data5[,5])

data1 <- read.table("E:/Academics/obj2-29March2019/Results/Hybrid/win_1_err.tsv")
data2 <- read.table("E:/Academics/obj2-29March2019/Results/Hybrid/win_3_err.tsv")
data3 <- read.table("E:/Academics/obj2-29March2019/Results/Hybrid/win_5_err.tsv")
data4 <- read.table("E:/Academics/obj2-29March2019/Results/Hybrid/win_7_err.tsv")
data5 <- read.table("E:/Academics/obj2-29March2019/Results/Hybrid/win_9_err.tsv")

hybrid_mae <- cbind(data1[,2], data2[,2], data3[,2], data4[,2], data5[,2])
hybrid_cor <- cbind(data1[,5], data2[,5], data3[,5], data4[,5], data5[,5])


# All Atom error ----------------------------------------------------------

# layout(matrix(1:2,ncol =2,byrow = TRUE),c(4,1), TRUE)
main_chain <- cbind(seq_mae[1,], dynaseq_mae[1,], hybrid_mae[1,])
barplot(main_chain, beside = T, col = c("red", "green", "yellow","blue", "orange"), main = "All Atom", xlab = "Model")
axis(1, at = c(3.5,9.5,15.5), labels = c("Seq Model", "DynaSeq Model","Hybrid Model"), cex.axis = 0.8,tick = F)
# legend(18,8,c("Win 1", "Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")


# Total Side error ----------------------------------------------------------

# layout(matrix(1:2,ncol =2,byrow = TRUE),c(4,1), TRUE)
main_chain <- cbind(seq_mae[2,], dynaseq_mae[2,], hybrid_mae[2,])
barplot(main_chain, beside = T, col = c("red", "green", "yellow","blue", "orange"), main = "Total Side", xlab = "Model")
axis(1, at = c(3.5,9.5,15.5), labels = c("Seq Model", "DynaSeq Model","Hybrid Model"), cex.axis = 0.8,tick = F)
# legend(18,8,c("Win 1", "Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")

# Main Chain error ----------------------------------------------------------

# layout(matrix(1:2,ncol =2,byrow = TRUE),c(4,1), TRUE)
main_chain <- cbind(seq_mae[3,], dynaseq_mae[3,], hybrid_mae[3,])
barplot(main_chain, beside = T, col = c("red", "green", "yellow","blue", "orange"), main = "Main Chain", xlab = "Model")
axis(1, at = c(3.5,9.5,15.5), labels = c("Seq Model", "DynaSeq Model","Hybrid Model"), cex.axis = 0.8,tick = F)
# legend(18,8,c("Win 1", "Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")

# Non Polar error ----------------------------------------------------------

# layout(matrix(1:2,ncol =2,byrow = TRUE),c(4,1), TRUE)
main_chain <- cbind(seq_mae[4,], dynaseq_mae[4,], hybrid_mae[4,])
barplot(main_chain, beside = T, col = c("red", "green", "yellow","blue", "orange"), main = "Non Polar", xlab = "Model")
axis(1, at = c(3.5,9.5,15.5), labels = c("Seq Model", "DynaSeq Model","Hybrid Model"), cex.axis = 0.8,tick = F)
# legend(18,8,c("Win 1", "Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")

# All Polar error ----------------------------------------------------------

# layout(matrix(1:2,ncol =2,byrow = TRUE),c(4,1), TRUE)
main_chain <- cbind(seq_mae[5,], dynaseq_mae[5,], hybrid_mae[5,])
barplot(main_chain, beside = T, col = c("red", "green", "yellow","blue", "orange"), main = "All Polar", xlab = "Model")
axis(1, at = c(3.5,9.5,15.5), labels = c("Seq Model", "DynaSeq Model","Hybrid Model"), cex.axis = 0.8,tick = F)
# legend(18,8,c("Win 1", "Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")


plot.new()
legend("center",c("Win 1", "Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")



#####################################################################################################################################################################



# All Atom error ----------------------------------------------------------

# layout(matrix(1:2,ncol =2,byrow = TRUE),c(4,1), TRUE)
main_chain <- cbind(seq_cor[1,], dynaseq_cor[1,], hybrid_cor[1,])
barplot(main_chain, beside = T, col = c("red", "green", "yellow","blue", "orange"), main = "All Atom", xlab = "Model")
axis(1, at = c(3.5,9.5,15.5), labels = c("Seq Model", "DynaSeq Model","Hybrid Model"), cex.axis = 0.8,tick = F)
# legend(18,8,c("Win 1", "Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")


# Total Side error ----------------------------------------------------------

# layout(matrix(1:2,ncol =2,byrow = TRUE),c(4,1), TRUE)
main_chain <- cbind(seq_cor[2,], dynaseq_cor[2,], hybrid_cor[2,])
barplot(main_chain, beside = T, col = c("red", "green", "yellow","blue", "orange"), main = "Total Side", xlab = "Model")
axis(1, at = c(3.5,9.5,15.5), labels = c("Seq Model", "DynaSeq Model","Hybrid Model"), cex.axis = 0.8,tick = F)
# legend(18,8,c("Win 1", "Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")

# Main Chain error ----------------------------------------------------------

# layout(matrix(1:2,ncol =2,byrow = TRUE),c(4,1), TRUE)
main_chain <- cbind(seq_cor[3,], dynaseq_cor[3,], hybrid_cor[3,])
barplot(main_chain, beside = T, col = c("red", "green", "yellow","blue", "orange"), main = "Main Chain", xlab = "Model")
axis(1, at = c(3.5,9.5,15.5), labels = c("Seq Model", "DynaSeq Model","Hybrid Model"), cex.axis = 0.8,tick = F)
# legend(18,8,c("Win 1", "Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")

# Non Polar error ----------------------------------------------------------

# layout(matrix(1:2,ncol =2,byrow = TRUE),c(4,1), TRUE)
main_chain <- cbind(seq_cor[4,], dynaseq_cor[4,], hybrid_cor[4,])
barplot(main_chain, beside = T, col = c("red", "green", "yellow","blue", "orange"), main = "Non Polar", xlab = "Model")
axis(1, at = c(3.5,9.5,15.5), labels = c("Seq Model", "DynaSeq Model","Hybrid Model"), cex.axis = 0.8,tick = F)
# legend(18,8,c("Win 1", "Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")

# All Polar error ----------------------------------------------------------

# layout(matrix(1:2,ncol =2,byrow = TRUE),c(4,1), TRUE)
main_chain <- cbind(seq_cor[5,], dynaseq_cor[5,], hybrid_cor[5,])
barplot(main_chain, beside = T, col = c("red", "green", "yellow","blue", "orange"), main = "All Polar", xlab = "Model")
axis(1, at = c(3.5,9.5,15.5), labels = c("Seq Model", "DynaSeq Model","Hybrid Model"), cex.axis = 0.8,tick = F)
# legend(18,8,c("Win 1", "Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")


plot.new()
legend("center",c("Win 1", "Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")



data <- read.table("E:/Academics/obj2-29March2019/Results/DynaSeq/win_1_predscore.tsv", sep = "\t")



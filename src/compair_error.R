

data1 <- read.table("E:/Academics/obj2-29March2019/Results/DynaSeq/win_1_err.tsv")
data2 <- read.table("E:/Academics/obj2-29March2019/Results/DynaSeqPair/err.tsv")
data3 <- read.table("E:/Academics/obj2-29March2019/Results/DynaSeqStep/err.tsv")

data <- cbind(data1[,5], data2[,5], data3[,5])
barplot(t(data), col = c("red", "green","blue"),main = "R^2 Comparison", sub = "DynaSeq vs DynaSeqPair vs DynaSeqStep", ylab = "R^2 Coeeficient", beside = TRUE, ylim = c(0,0.50), names.arg = rownames(data1), cex.names = 0.71, legend.text = c("Dynaseq", "Dynaseqpair", "Dynaseqstep"), args.legend = list(x = 20, cex = 0.5))


data <- cbind(data1[,2], data2[,2], data3[,2])
barplot(t(data), col = c("red", "green","blue"),main = "MAE Comparison", sub = "DynaSeq vs DynaSeqPair vs DynaSeqStep", ylab = "Mean Absolute Error", beside = TRUE, ylim = c(0,50), names.arg = rownames(data1), cex.names = 0.71, legend.text = c("Dynaseq", "Dynaseqpair", "Dynaseqstep"), args.legend = list(x = 20, cex = 0.5))



data1 <- read.table("E:/Academics/obj2-29March2019/Results/Hybrid/win_1_err.tsv")
data2 <- read.table("E:/Academics/obj2-29March2019/Results/HybridPair/err.tsv")
data3 <- read.table("E:/Academics/obj2-29March2019/Results/HybridStep/err.tsv")

data <- cbind(data1[,5], data2[,5], data3[,5])
barplot(t(data), col = c("red", "green","blue"), main = "R^2 Comparison", sub = "Hybrid vs HybridPair vs HybridStep", ylab = "R^2 Coeeficient", beside = TRUE, ylim = c(0,0.9), names.arg = rownames(data1), cex.names = 0.71, legend.text = c("Hybrid", "Hybridpair", "Hybridstep"), args.legend = list(x = 20, cex = 0.6))


data <- cbind(data1[,2], data2[,2], data3[,2])
barplot(t(data), col = c("red", "green","blue"), main = "MAE Comparison", sub = "Hybrid vs HybridPair vs HybridStep", ylab = "Mean Absolute Error", beside = TRUE, ylim = c(0,50), names.arg = rownames(data1), cex.names = 0.71, legend.text = c("Hybrid", "Hybridpair", "Hybridstep"), args.legend = list(x = 20, cex = 0.6))


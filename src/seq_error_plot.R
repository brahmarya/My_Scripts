

data1 <- read.table("E:/Academics/obj2-29March2019/Results/Seq/win_3_err.tsv")
data2 <- read.table("E:/Academics/obj2-29March2019/Results/Seq/win_5_err.tsv")
data3 <- read.table("E:/Academics/obj2-29March2019/Results/Seq/win_7_err.tsv")
data4 <- read.table("E:/Academics/obj2-29March2019/Results/Seq/win_9_err.tsv")


pdf("E:/Academics/obj2-29March2019/Results/Seq model error plot.pdf")

# Plot for SD ---------------------------------------------------

nf <- layout(matrix(1:6,ncol = 3,byrow = TRUE), TRUE)
# layout.show(nf)
# par(oma = c(1,1,1,1))


a <- min(c(data1[,1],data2[,1],data3[,1],data4[,1]))
b <- max(c(data1[,1],data2[,1],data3[,1],data4[,1]))

plot(data1[,1], type = "b", col = "red", xlab = "Atom Types", ylab = "Stndard Deviation", xlim = c(1,5),ylim = c(a,b), xaxt = "n", main = "Standard Deviation")
lines(data2[,1], type = "b", col = "green")
lines(data3[,1], type = "b", col = "yellow")
lines(data4[,1], type = "b", col = "blue")
# legend("bottomright",c("Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue"), xpd = NA, bty = "n")
axis(1, at = c(1,2,3,4,5), labels = c("All_Atom","Total_Side","Main_Chain","Non_Polar","All_Polar"), cex.axis = 0.8)



# Plot for MAE  -----------------------------------------------------------

a <- min(c(data1[,2],data2[,2],data3[,2],data4[,2]))
b <- max(c(data1[,2],data2[,2],data3[,2],data4[,2]))

plot(data1[,2], type = "b", col = "red", xlab = "Atom Types", ylab = "MAE", xlim = c(1,5),ylim = c(a,b), xaxt = "n", main = "Mean Absolute Error")
lines(data2[,2], type = "b", col = "green")
lines(data3[,2], type = "b", col = "yellow")
lines(data4[,2], type = "b", col = "blue")
# legend("bottomright",c("Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")
axis(1, at = c(1,2,3,4,5), labels = c("All_Atom","Total_Side","Main_Chain","Non_Polar","All_Polar"), cex.axis = 0.8)


# Plot for MSE -----------------------------------------------------------

a <- min(c(data1[,3],data2[,3],data3[,3],data4[,3]))
b <- max(c(data1[,3],data2[,3],data3[,3],data4[,3]))

plot(data1[,3], type = "b", col = "red", xlab = "Atom Types", ylab = "MSE", xlim = c(1,5),ylim = c(a,b), xaxt = "n", main = "Mean Square Error")
lines(data2[,3], type = "b", col = "green")
lines(data3[,3], type = "b", col = "yellow")
lines(data4[,3], type = "b", col = "blue")
# legend("bottomright",c("Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")
axis(1, at = c(1,2,3,4,5), labels = c("All_Atom","Total_Side","Main_Chain","Non_Polar","All_Polar"), cex.axis = 0.8)


# Plot for RMSE -----------------------------------------------------------

a <- min(c(data1[,4],data2[,4],data3[,4],data4[,4]))
b <- max(c(data1[,4],data2[,4],data3[,4],data4[,4]))

plot(data1[,4], type = "b", col = "red", xlab = "Atom Types", ylab = "RMSE", xlim = c(1,5),ylim = c(a,b), xaxt = "n", main = "Root Mean Square Error")
lines(data2[,4], type = "b", col = "green")
lines(data3[,4], type = "b", col = "yellow")
lines(data4[,4], type = "b", col = "blue")
# legend("bottomright",c("Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")
axis(1, at = c(1,2,3,4,5), labels = c("All_Atom","Total_Side","Main_Chain","Non_Polar","All_Polar"), cex.axis = 0.8)



# Plot for RSQ -----------------------------------------------------------

a <- min(c(data1[,5],data2[,5],data3[,5],data4[,5]))
b <- max(c(data1[,5],data2[,5],data3[,5],data4[,5]))

plot(data1[,5], type = "b", col = "red", xlab = "Atom Types", ylab = "RSQ", xlim = c(1,5),ylim = c(a,b), xaxt = "n", main = "Coefficient of Determination")
lines(data2[,5], type = "b", col = "green")
lines(data3[,5], type = "b", col = "yellow")
lines(data4[,5], type = "b", col = "blue")
# legend("bottomright",c("Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue","orange"), xpd = NA, bty = "n")
axis(1, at = c(1,2,3,4,5), labels = c("All_Atom","Total_Side","Main_Chain","Non_Polar","All_Polar"), cex.axis = 0.8)



plot.new()
legend("center",c("Win 3", "Win 5", "Win 7", "Win 9"), fill = c("red","green","yellow","blue"), xpd = NA, bty = "n")


dev.off()


# ggplots of errors -------------------------------------------------------

ggplot(b,aes(Win,rmse, fill=Win)) + geom_col(width=0.3) + ggtitle(label = "Root Mean Square Erorr") + theme_classic() + theme(plot.title = element_text(size = 20,colour = "blue", hjust = 0.5))

# ggplot() + geom_point(data = data, aes(x = All.atoms, y = X),size=4)

smoothScatter(data$All.atoms~data$X, col = c("red","green"), cex = 3, xlab = "Actual All_atoms ASA", ylab = "Predicted All_atoms ASA", main = "Scatter plot Actual vs Predicted ASA", cex.axis = 0.8)



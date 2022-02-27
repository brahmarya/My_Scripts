

# Load required packages --------------------------------------------------

library(glmnet)
library(caret)
library(Metrics)
library(dplyr)

files <- dir("E:/Academics/obj2-29March2019/RCSB/pat_step/", full.names = TRUE);
data <- c();
lapply(1:length(files), function(x){
  df <- read.table(files[x], sep = "\t", header = TRUE)
  data <<- rbind(data,df )
  return(0)
})

data <- na.omit(data)
all_pred <- c()
lapply(1:nrow(data), function(n){
  train_data <- data[-n,]
  train <- train_data[!(duplicated(train_data[,12:35]) | duplicated(train_data[,12:35], fromLast = TRUE)),]
  
  x <- as.matrix(train[,12:35])
  y <- as.matrix(train[,(ncol(data)-4):ncol(data)])
  
  test_in <- as.matrix(data[n,12:35])
  test_t <- data[n,(ncol(data)-4):ncol(data)]
  
  cvmfit <- cv.glmnet(x, y, family = "mgaussian")
  
  mpred <- predict(cvmfit, newx = test_in, s = cvmfit$lambda.min)   
  all_pred <<- rbind(all_pred, mpred)
})

actual <- as.matrix(data[,(ncol(data)-4):ncol(data)])

score <- cbind(actual,all_pred);
write.table(score,file = "E:/Academics/obj2-29March2019/Results/DynaSeqStep/predscore.tsv", sep = "\t");

rmse <- unlist(lapply(1:5, function(n){r <- rmse(actual[,n], all_pred[,n]);return(r)}))
mse <- unlist(lapply(1:5, function(n){r <- mse(actual[,n], all_pred[,n]);return(r)}))
mae <- unlist(lapply(1:5, function(n){r <- mae(actual[,n], all_pred[,n]);return(r)}))
sd <- unlist(lapply(1:5, function(n){r <- sqrt(var(actual[,n]) + var(all_pred[,n]));return(r)}))
rsq <- unlist(lapply(1:5, function(n){r<- cor(actual[,n], all_pred[,n])^2;return(r)}))


err <- cbind(sd,mae,mse,rmse,rsq)
row.names(err) <- c("All_Atom","Total_Side","Main_Chain","Non_Polar","All_Polar")

write.table(err, file = "E:/Academics/obj2-29March2019/Results/DynaSeqStep/err.tsv", sep = "\t")





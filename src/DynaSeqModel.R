

# Load required packages --------------------------------------------------

library(glmnet)
library(caret)
library(Metrics)
library(dplyr)

win <- c(1,3,5,7,9)

lapply(1:5, function(m){
  data <- read.table(paste("/storage/obj2-29March2019/RCSB/pat/win_",win[m],"_A.tsv", sep = ""), sep = "\t", header = TRUE)
  data <- data[,((win[m]*4)+4):ncol(data)]
  data <- na.omit(data)
  all_pred <- c()
  lapply(1:nrow(data), function(n){
    train_data <- data[-n,]
    train <- train_data[!(duplicated(train_data[,1:(ncol(data) - 5)]) | duplicated(train_data[,1:(ncol(data) - 5)], fromLast = TRUE)),]
    
    x <- as.matrix(train[,1:(ncol(data) - 5)])
    y <- as.matrix(train[,(ncol(data)-4):ncol(data)])
    
    test_in <- as.matrix(data[n,1:(ncol(data) - 5)])
    test_t <- data[n,(ncol(data)-4):ncol(data)]
    
    cvmfit <- cv.glmnet(x, y, family = "mgaussian")
    saveRDS(cvmfit, file = paste("/storage/obj2-29March2019/Model/DynaSeq/win_",win[m],".RDS", sep = ""))
    mpred <- predict(cvmfit, newx = test_in, s = cvmfit$lambda.min)  
    all_pred <<- rbind(all_pred, mpred)
  })
  
  actual <- as.matrix(data[,(ncol(data)-4):ncol(data)])
  
  score <- cbind(actual,all_pred);
  write.table(score,paste("/storage/obj2-29March2019/Results/DynaSeq/win_",win[m],"_predscore.tsv", sep = ""), sep = "\t");
  
  rmse <- unlist(lapply(1:5, function(n){r <- rmse(actual[,n], all_pred[,n]);return(r)}))
  mse <- unlist(lapply(1:5, function(n){r <- mse(actual[,n], all_pred[,n]);return(r)}))
  mae <- unlist(lapply(1:5, function(n){r <- mae(actual[,n], all_pred[,n]);return(r)}))
  sd <- unlist(lapply(1:5, function(n){r <- sqrt(var(actual[,n]) + var(all_pred[,n]));return(r)}))
  rsq <- unlist(lapply(1:5, function(n){r<- cor(actual[,n], all_pred[,n])^2;return(r)}))
  
  
  err <- cbind(sd,mae,mse,rmse,rsq)
  row.names(err) <- c("All_Atom","Total_Side","Main_Chain","Non_Polar","All_Polar")
  
  write.table(err, file = paste("/storage/obj2-29March2019/Results/DynaSeq/win_",win[m],"_err.tsv", sep = ""), sep = "\t")
})





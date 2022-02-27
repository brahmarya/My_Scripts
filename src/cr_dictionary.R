

load("/storage/ICMR/dictionary.static.Rdata")

dynaseq.str <- dynaseq.str[1:1024,1:12]

a <- c("A","T","G","C")
b <- expand.grid(a,a,a,a,a)
c <- paste(b$Var1, b$Var2, b$Var3, b$Var4, b$Var5, sep = "")
rownames(dynaseq.str) <- c


model <- readRDS("/storage/obj2-29March2019/Model/DynaSeq/win_1.RDS")
mpred <- predict(model, newx = dynaseq.str, s = model$lambda.min)
static.asa <- mpred[,,1]
save(static.asa, file = "/storage/obj2-29March2019/Dictionary/Static.win1.Rdata")

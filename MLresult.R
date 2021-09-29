#This code is designed to give a comprehensive results with multiple machine learning methods. It takes the standard traindata and
# valdata as the inputs, and returns with the Sensitivity, Specificity as well as the accuracy of training and validation datasets.
# @Input: traindata and valdata, are the standard input, with the rows indicated the samples, the first column is the type ]
#         of the samples, and the other columns are the variables
# @output: a vector. 

library(randomForest)
library(pROC)
library(e1071)
library(nnet)
library(MASS)
library(mda)
library(xgboost)
library(catboost)
library(caret)
MLresult = function(traindata, valdata){
  #Logistic regression
  #train
  glm.fit = glm(Type ~. , data = traindata, family ="binomial")
  logistic.rocobj  = roc(traindata$Type, as.vector(t(predict(glm.fit))),smooth = FALSE)
  logistic.rocdata = data.frame(Sens = logistic.rocobj$sensitivities, Spec = logistic.rocobj$specificities)
  logistic.rocdata$Combine = logistic.rocdata[,1] + logistic.rocdata[,2]
  seq.bestcutoff = which( logistic.rocdata$Combine == max(logistic.rocdata$Combine))[1]
  bestcutoff = logistic.rocobj$thresholds[seq.bestcutoff]
  predType = cut(as.vector(t(predict(glm.fit))), breaks =c(-Inf, bestcutoff, Inf), labels=c(0,1))
  cTab = table(traindata$Type, predType, dnn=c("actual", "predicted"))
  Sens.log.tr = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.log.tr = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.log.tr = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  #test
  response = as.vector(t(predict(glm.fit, valdata[,-1])))
  predType = cut(response,breaks =c(-Inf, bestcutoff, Inf), labels=c(0,1))
  cTab = table(valdata$Type, predType, dnn=c("actual", "predicted"))
  Sens.log.vl = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.log.vl = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.log.vl = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  #RandomForest
  #train
  traindata$Type = as.factor(traindata$Type)
  rf = randomForest(x = traindata[,-1], y = traindata[,1], ntree=500, keep.forest=TRUE)
  cTab = rf$confusion
  Sens.rf.tr = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.rf.tr = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.rf.tr = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  #test
  pred = predict(rf, valdata[,-1])
  cTab = table(valdata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.rf.vl = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.rf.vl = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.rf.vl = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  #SVM
  #train
  svm.model = svm(Type~. , data = traindata)
  cTab = table(traindata[,1], fitted(svm.model), dnn=c("Acutal","Predicted"))
  Sens.svm.tr = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.svm.tr = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.svm.tr = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  #test
  pred = predict(svm.model, valdata)
  cTab = table(valdata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.svm.vl = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.svm.vl = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.svm.vl = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  #NaiveBayes
  #train
  bayes.model = naiveBayes(Type~. , data = traindata)
  cTab = table(traindata[,1], predict(bayes.model,traindata), dnn=c("Acutal","Predicted"))
  Sens.NB.tr = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.NB.tr = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.NB.tr = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  #test
  pred = predict(bayes.model, valdata)
  cTab = table(valdata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.NB.vl = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.NB.vl = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.NB.vl = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  #Neural Network
  #train
  nn.model = nnet(Type~., size=3,maxit = 10000, decay=0.01, rang = 0.05, data = traindata, trace = FALSE)
  cTab = table(traindata[,1], predict(nn.model, traindata, type="class"), dnn=c("Acutal","Predicted"))
  Sens.NN.tr = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.NN.tr = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.NN.tr = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  #test
  pred = predict(nn.model, valdata, type="class")
  cTab = table(valdata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.NN.vl = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.NN.vl = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.NN.vl = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  #Linear Discriminant Analysis
  #train
  ldamod = lda(traindata[,-1], traindata[,1])
  pred   = predict(ldamod, traindata[,-1])$class
  cTab   = table(traindata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.LDA.tr = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.LDA.tr = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.LDA.tr = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  #test
  pred   = predict(ldamod, valdata[,-1])$class
  cTab   = table(valdata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.LDA.vl = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.LDA.vl = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.LDA.vl = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  #Mixture Discriminant Analysis
  #train
  mdamod = mda(traindata$Type~., data = traindata)
  pred   = predict(mdamod, traindata[,-1])
  cTab   = table(traindata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.MDA.tr = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.MDA.tr = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.MDA.tr = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  #test
  pred   = predict(mdamod, valdata[,-1])
  cTab   = table(valdata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.MDA.vl = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.MDA.vl = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.MDA.vl = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  #Flexible Discriminant Analysis
  #train
  fdamod = fda(traindata$Type~. , data = traindata)
  pred   = predict(fdamod, traindata[,-1])
  cTab   = table(traindata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.FDA.tr = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.FDA.tr = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.FDA.tr = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  #test
  pred   = predict(fdamod, valdata[,-1])
  cTab   = table(valdata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.FDA.vl = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.FDA.vl = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.FDA.vl = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  # Gradient boosting machine
  #train 
  Type = as.factor(ifelse(traindata$Type == 1, "Case", "Control"))
  objControl <- trainControl(method='cv', number=3, returnResamp='none', summaryFunction = twoClassSummary, classProbs = TRUE)
  gbmmod = train(x= traindata[,-1], y = Type, method='gbm',  trControl=objControl,  metric = "ROC",preProc = c("center", "scale"))
  pred = predict(gbmmod, traindata[,-1])
  pred = ifelse(pred == "Case", 1, 0)
  cTab   = table(traindata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.GBM.tr = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.GBM.tr = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.GBM.tr = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  #test
  pred   = predict(gbmmod, valdata[,-1])
  pred   = ifelse(pred == "Case", 1, 0)
  cTab   = table(valdata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.GBM.vl = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.GBM.vl = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.GBM.vl = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  
  #XGboost
  #train
  Type = as.numeric(as.character(traindata$Type))
  bstSparse <- xgboost(data = as.matrix(traindata[,-1]), label =Type, max.depth = 6, eta = 0.3, nthread = 8, nrounds = 100, objective = "binary:logistic")
  pred = predict(bstSparse, as.matrix(traindata[,-1]))
  pred = as.numeric(pred > 0.5)
  cTab = table(traindata$Type, pred, dnn=c("actual", "predicted"))
  Sens.xgb.tr = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.xgb.tr = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.xgb.tr = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  #test
  pred   = predict(bstSparse, as.matrix(valdata[,-1]))
  pred = as.numeric(pred > 0.5)
  cTab   = table(valdata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.xgb.vl = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.xgb.vl = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.xgb.vl = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  #catboost
    
  #train
  train_pool = catboost.load_pool(data = as.matrix(traindata[,-1]), label = as.numeric(as.character(traindata$Type)))
  params <- list(iterations=200,
                 learning_rate=0.08,
                 depth=10,
                 loss_function='RMSE',
                 eval_metric='RMSE',
                 random_seed = 55,
                 od_type='Iter',
                 metric_period = 50,
                 od_wait=20,
                 use_best_model=FALSE)
  catboost_model <- catboost.train(learn_pool = train_pool, test_pool = NULL,params)
  pred <- catboost.predict(catboost_model, train_pool)
  pred = as.numeric(pred > 0.5)
  cTab = table(traindata$Type, pred, dnn=c("actual", "predicted"))
  Sens.catb.tr = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.catb.tr = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.catb.tr = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  #test
  val_pool = catboost.load_pool(data = as.matrix(valdata[,-1]), label = as.numeric(as.character(valdata$Type)))
  pred <- catboost.predict(catboost_model, val_pool)
  pred = as.numeric(pred > 0.5)
  cTab   = table(valdata[,1], pred, dnn=c("Acutal","Predicted"))
  Sens.catb.vl = cTab[2,2]/(cTab[2,2] + cTab[2,1])
  Spec.catb.vl = cTab[1,1]/(cTab[1,1] + cTab[1,2])
  Accu.catb.vl = (cTab[1,1] + cTab[2,2])/sum(cTab[1:2,1:2])
  
  
  return( c(Sens.log.tr, Spec.log.tr, Accu.log.tr, Sens.log.vl, Spec.log.vl, Accu.log.vl, 
            Sens.rf.tr,  Spec.rf.tr,  Accu.rf.tr,  Sens.rf.vl,  Spec.rf.vl,  Accu.rf.vl ,
            Sens.svm.tr, Spec.svm.tr, Accu.svm.tr, Sens.svm.vl, Spec.svm.vl, Accu.svm.vl,
            Sens.NB.tr , Spec.NB.tr,  Accu.NB.tr , Sens.NB.vl , Spec.NB.vl , Accu.NB.vl,
            Sens.NN.tr , Spec.NN.tr,  Accu.NN.tr , Sens.NN.vl , Spec.NN.vl,  Accu.NN.vl,
            Sens.LDA.tr, Spec.LDA.tr, Accu.LDA.tr, Sens.LDA.vl, Spec.LDA.vl, Accu.LDA.vl,
            Sens.MDA.tr, Spec.MDA.tr, Accu.MDA.tr, Sens.MDA.vl, Spec.MDA.vl, Accu.MDA.vl,
            Sens.FDA.tr, Spec.FDA.tr, Accu.FDA.tr, Sens.FDA.vl, Spec.FDA.vl, Accu.FDA.vl,
            Sens.GBM.tr, Spec.GBM.tr, Accu.GBM.tr, Sens.GBM.vl, Spec.GBM.vl, Accu.GBM.vl,
            Sens.xgb.tr, Spec.xgb.tr, Accu.xgb.tr, Sens.xgb.vl, Spec.xgb.vl, Accu.xgb.vl,
            Sens.catb.tr, Spec.catb.tr, Accu.catb.tr, Sens.catb.vl, Spec.catb.vl, Accu.catb.vl
            ))
}


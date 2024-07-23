    setwd("/data5/shihong/fengmh/data/CWJ/BC")
    T_sce<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/癌细胞亚群/sub_cancer_umap手动.rds")
    T_sce[["RNA"]] <- as(T_sce[["RNA"]], "Assay")
    t_expr_data <- t(as.matrix(T_sce@assays$RNA@data))
    dim(t_expr_data)
    T_sce_markers<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/癌细胞亚群/markers_subcancer_手动.rds")
    t_expr<-t_expr_data[,which(colnames(t_expr_data) %in% rownames(T_sce_markers))]
    dim(t_expr)
    T_sce <- SetIdent(T_sce, value = T_sce$group)
    levels(Idents(T_sce))
    inTrain<-createDataPartition(y= Idents(T_sce) ,p=0.3,list=F)
    test_expr <-t_expr[inTrain,]
    train_expr <-t_expr[-inTrain,]

    test_y <- Idents(T_sce)[inTrain]
    train_y  <- Idents(T_sce)[-inTrain] 
    save(test_expr,train_expr,test_y,train_y,file="/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/SVM/k折交叉验证/input2.Rdata")

    ############################random forest
    setwd("/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能")
    library(randomForest)
    library(caret)
    rf_imp<-read.csv(file="SVM/k折交叉验证/rf_imp.csv")
    load(file="/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/SVM/k折交叉验证/rf_rfcv_model_rgcv.Rdata")
    rownames(rf_imp)<-rf_imp[,1]
    rf_imp<-rf_imp[,-1]
    rf_rfcv_model$error.cv[which(rf_rfcv_model$error.cv==min(rf_rfcv_model$error.cv))]
    rf_imp_kfold<-rf_imp[1:549,]
    load(file="/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/SVM/k折交叉验证/input2.Rdata")

    predictor_data = train_expr[,match(rownames(rf_imp_kfold),colnames(train_expr))]
    target = train_y 

    rf_k_fold=randomForest(x=predictor_data, y=target,
                        importance = TRUE, ntree = 200)
    save(rf_k_fold,file="SVM/k折交叉验证/rf_k_fold_new.Rdata")
    load("SVM/k折交叉验证/rf_k_fold_new.Rdata")
    test_outputs <- predict(rf_k_fold,newdata = test_expr,type="prob")
    test_expr[1:4,1:4]
    summary(test_outputs)
    head(test_outputs)
    pred_y = colnames(test_outputs)[apply(test_outputs, 1, which.max)]
    pred_y = factor(pred_y,levels = levels(test_y))
    confusion_matrix <- confusionMatrix(pred_y, test_y)

    accuracy <- confusion_matrix$overall["Accuracy"]
    pdf('SVM/k折交叉验证/RF_confusionMatrix_tree200_kfold.pdf',width = 10)
    p<-gplots::balloonplot(table(pred_y,test_y),main="RandomForest Confusion Matrix")
    print(p)
    dev.off()
    pdf('/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/RF_confusionMatrix.pdf',width = 10)
    p<-gplots::balloonplot(table(pred_y,test_y),main="RandomForest Confusion Matrix")
    print(p)
    dev.off()
    write.csv(test_outputs,file="/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/RF_pred_y.csv")

    #######################LASSO
    library(tidyverse)
    library(glmnet)
    library(caret)
    library(ggplot2)
    setwd("/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能")
    set.seed(12345)

    fit<-readRDS(file="SVM/LASSO/lasso_model_new.rds")
    load(file="/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/SVM/k折交叉验证/input2.Rdata")
    x <- train_expr
    y <- train_y
    cvfit = cv.glmnet(x, y,nfold=10,
                    family = "multinomial", type.measure = "class")
    cvfit<-readRDS("SVM/LASSO/10_folds_LASSO_model_new.rds")
    test_outputs <-  predict(cvfit,  as.matrix(test_expr) , 
                            type="response", s="lambda.1se") 
    head( test_outputs ) 
    pred_y = colnames(test_outputs)[apply(test_outputs, 1, which.max)]
    pred_y = factor(pred_y,levels = levels(test_y))
    pdf('/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/lasso_ConfusionMatrix.pdf',width = 10)
    p<-gplots::balloonplot(table(pred_y,test_y),main="Lasso Confusion Matrix")
    print(p)
    dev.off()
    confusion_matrix <- confusionMatrix(pred_y, test_y)
    accuracy <- confusion_matrix$overall["Accuracy"]
    pdf('SVM/k折交叉验证/LASSO_ROC_kfold.pdf',width=10)
    p<-plot(1-confusion_matrix$byClass[,2],confusion_matrix$byClass[,2],type="l",col="red",lty=1,xlab = "1-Specificity",ylab = "Sensitivities",lwd=2)
    print(p)
    dev.off()

    {
        levels(pred_y)<-c(0,1,2,3,4)
        pred_y[which(pred_y=="Ta")]<-0
        pred_y[which(pred_y=="T1")]<-1
        pred_y[which(pred_y=="T2")]<-2
        pred_y[which(pred_y=="T3")]<-3
        pred_y[which(pred_y=="T4")]<-4
        levels(test_y)<-c(0,1,2,3,4)
        test_y[which(test_y=="Ta")]<-0
        test_y[which(test_y=="T1")]<-1
        test_y[which(test_y=="T2")]<-2
        test_y[which(test_y=="T3")]<-3
        test_y[which(test_y=="T4")]<-4
        test_outputs<-test_outputs[,,1]
        write.csv(test_outputs,file="/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/LASSO_pred_y.csv")
        write.csv(test_y,file="/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/test_y.csv")
    }
    cvfit$lambda.min

    myCoefs <- coef(cvfit, s="lambda.min")
    str(myCoefs[1])
    str(myCoefs[2])
    str(myCoefs[3])
    str(myCoefs[4])
    str(myCoefs[5])
    Ta_ceofs<-myCoefs[1]
    T1_ceofs<-myCoefs[2]
    T2_ceofs<-myCoefs[3]
    T3_ceofs<-myCoefs[4]
    T4_ceofs<-myCoefs[5]

    lasso_fea_T2 <- T2_ceofs$T2 %>% as.matrix(.) %>% .[which(.[,1]!=0),] %>% .[order(., decreasing = T)] %>% as.matrix(.)
    write.csv(lasso_fea_T2,"/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/lasso/feature_lasso_T2.csv")
    lasso_fea_T3 <- T3_ceofs$T3 %>% as.matrix(.) %>% .[which(.[,1]!=0),] %>% .[order(., decreasing = T)] %>% as.matrix(.)
    write.csv(lasso_fea_T3,"/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/lasso/feature_lasso_T3.csv")
    lasso_fea_T4 <- T4_ceofs$T4 %>% as.matrix(.) %>% .[which(.[,1]!=0),] %>% .[order(., decreasing = T)] %>% as.matrix(.)
    write.csv(lasso_fea_T4,"/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/lasso/feature_lasso_T4.csv")
    lasso_fea_T1 <- T1_ceofs$T1 %>% as.matrix(.) %>% .[which(.[,1]!=0),] %>% .[order(., decreasing = T)] %>% as.matrix(.)
    write.csv(lasso_fea_T1,"/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/lasso/feature_lasso_T1.csv")
    lasso_fea_Ta <- Ta_ceofs$Ta %>% as.matrix(.) %>% .[which(.[,1]!=0),] %>% .[order(., decreasing = T)] %>% as.matrix(.)
    write.csv(lasso_fea_Ta,"/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/lasso/feature_lasso_Ta.csv")

    #################Xgboost
    library(xgboost)
    setwd("/data5/shihong/fengmh/data/CWJ/BC/first_text/整体seurat对象机器学习_TCGA验证")
    library(Seurat)
    library(dplyr)
    library(caret)
    set.seed(123)
    cancer_umap<-readRDS("/data5/shihong/fengmh/data/CWJ/BC/单细胞浸润V非浸润数据/HC/R/T1T2T3T4_2/癌细胞亚群/sub_cancer_umap手动.rds")
    cancer_umap[["RNA"]] <- as(cancer_umap[["RNA"]], "Assay")
    expr_data <- t(as.data.frame(cancer_umap@assays$RNA@data))
    group<-cancer_umap@meta.data$group
    raw_expr<-cbind(as.data.frame(cancer_umap@meta.data$group),expr_data)
    sample <- sample(c(TRUE, FALSE), nrow(raw_expr), replace=TRUE, prob=c(0.7,0.3))
    train  <- raw_expr[sample, ]
    test   <- raw_expr[!sample, ]
    train_y<-train[,1]
    train_x<-train[,-1]
    test_y<-test[,1]
    test_x<-test[,-1]
    levels(train_y)<-c(0,1,2,3,4)
    train_y[which(train_y=="Ta")]<-0
    train_y[which(train_y=="T1")]<-1
    train_y[which(train_y=="T2")]<-2
    train_y[which(train_y=="T3")]<-3
    train_y[which(train_y=="T4")]<-4
    levels(test_y)<-c(0,1,2,3,4)
    test_y[which(test_y=="Ta")]<-0
    test_y[which(test_y=="T1")]<-1
    test_y[which(test_y=="T2")]<-2
    test_y[which(test_y=="T3")]<-3
    test_y[which(test_y=="T4")]<-4
    xgb_train = xgb.DMatrix(data = as.matrix(train_x), label = (as.numeric(train_y)-1),nthread=12)
    xgb_test = xgb.DMatrix(data = as.matrix(test_x), label = (as.numeric(test_y)-1),nthread=12)
    watchlist = list(train=xgb_train, test=xgb_test)
    wacth_model = xgb.train(data = xgb_train, max.depth = 3, watchlist=watchlist, nrounds = 500, nthread=20)
    # [44]    train-rmse:0.440112     test-rmse:0.463950 
    # [45]    train-rmse:0.439322     test-rmse:0.463443 
    # [46]    train-rmse:0.438338     test-rmse:0.463206 
    # [47]    train-rmse:0.437841     test-rmse:0.463175 
    # [48]    train-rmse:0.437225     test-rmse:0.463260 
    # [49]    train-rmse:0.436797     test-rmse:0.463315 
    # [50]    train-rmse:0.436456     test-rmse:0.463548 
    # saveRDS(wacth_model,file="xgboost_output/xgboost_watch_model.rds")
    params <- list(objective = "multi:softprob",num_class =5, eval_metric = "mlogloss", eta = 0.1, max_depth = 3)
    nrounds <- 47
    xgb_model <- xgboost(params = params, data = xgb_train, nrounds = nrounds, nthread=20)
    saveRDS(xgb_model,file="xgboost_output/xgboost_cancer_model.rds")
    xgb_model<-readRDS("xgboost_output/xgboost_cancer_model.rds")
    pre<-predict(xgb_model,newdata = xgb_test)
    pred_prob_matrix <- matrix(pre, nrow = nrow(test_x), byrow = TRUE)
    rownames(pred_prob_matrix)<-rownames(test_x)
    colnames(pred_prob_matrix)<-c("0","1","2","3","4")
    write.csv(pred_prob_matrix,file="/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/xgboost_pred_y.csv")
    pred_xgboost = colnames(pred_prob_matrix)[apply(pred_prob_matrix, 1, which.max)]
    pred_xgboost = factor(pred_xgboost,levels = levels(test_y))

    tab = table(test_y,pre,dnn=c("true","pre"))
    tab = table(test_y,pred_xgboost,dnn=c("true","pre"))
    caret::confusionMatrix(tab)
    importance <- xgb.importance(feature_names = colnames(train_x), model = xgb_model)
    importance<-as.data.frame(importance)
    write.csv(importance,file="xgboost_importance.csv")
    pred_y<-pred_xgboost
    pdf('/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/xgboost-performance.pdf',width = 10)
    p<-gplots::balloonplot(table(pred_y,test_y),main="RandomForest Confusion Matrix")
    print(p)
    dev.off()
    {
        write.csv(pred_y,file="/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/xgboost_pred_y.csv")
        write.csv(test_y,file="/data5/shihong/fengmh/data/CWJ/BC/first_text/单细胞功能/机器学习补图/xgboost_test_y.csv")
    }

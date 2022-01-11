cat("
 Dscription: This R script implements Breiman's random forest algorithm (based
              on Breiman and Cutler's original Fortran code) for classification.
               
       Input: (1) Training data
              (2) Metadata of training data
              (3) Test data (The data you want to predict)
              (4) Metadata of test data
              (5) The factor category (i.e. status) both in the two metadata you 
                  provided
              
              Note: The non-overlaping variables between training and test data will
              be screen out for modeling. 
              
      Output: (1) Model performance in training data: Kappa statistics, ROC analysis, 
                      probability calculation.
              (2) The number of important variables required for maximum accurary of RF model.
              (3) The predicted probability of certain state in the samples of test data. 
             \n")
cat("
       Usage: setwd(\"$PWD\"); soucre(\"Raman.randomforest.pred-prob.simple.R\") \n\n")            
#################################################################
# Function:  RandomForest for classification of certain states of samples  
# Last update: 2021-11-01, Rongze Chen, Gongchao Jing
#################################################################
# install necessary libraries

## Clean R environment
rm(list=ls())
setwd('./')
sourcedir <- paste("/home/gene/jinggc/RamEX/databases/RBCS")
source(paste(sourcedir, "Raman_RF_util.R", sep = "/"))

#--------------------------------------------------
#opts<-list(
#                prefix="Ramandata", # The prefix of output filename related to training data
#                filename="Alldata.txt",
#                metadata.filename="Metadata.txt",
#                group.type="Group",
#                               
#                pred.prefix=NA, # The prefix of output filename related to test data. If you don't have test data, please omit.
#                pred.filename=NA, 
#                pred.metadata.filename=NA,
#                pred.group.type=NA,
#                alpha=0.0001,
#                device="png"
#                )
args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
        make_option(c("-p", "--prefix_of_training_output"),type="character", help="Input_the_prefix_training"),
        make_option(c("-i", "--filename_of_training"),type="character", help="Input_the_training_raman_data"),
        make_option(c("-m", "--metadata_of_training"),type="character", help="Input_the_training_metadata"),
        make_option(c("-g", "--grouptype_of_training"),type="character", default="Group", help="Input_the_training_group"),
        make_option(c("-o", "--prefix_of_prediction_output"),type="character", help="Input_the_prefix_prediction"),
        make_option(c("-f", "--filename_of_prediction"),type="character", default="NA", help="Input_the_prediction_raman_data [default %default]"),
        make_option(c("-d", "--metadata_of_prediction"),type="character", default="NA", help="Input_the_prediction_metadata [default %default]"),
        make_option(c("-t", "--grouptype_of_prediction"),type="character", default="NA", help="Input_the_prediction_group [default %default]"),
        make_option(c("-a", "--alpha"),type="double", default=0.0001, help="Input_the_prediction_group [default %default]"),
        make_option(c("-v", "--device"),type="character", default="png", help="Input_the_prediction_group [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)
#--------------------------------------------------
rf.opts<-list( ntree=1000, errortype = 'oob', verbose=FALSE, nfolds=3, rf.cv=FALSE)
#--------------------------------------------------
filename<-opts$filename_of_training
metadata.filename<-opts$metadata_of_training
pred.filename<-opts$filename_of_prediction
pred.metadata.filename<-opts$metadata_of_prediction
alpha<-opts$alpha
device<-opts$device

#--------------------------------------------------
# Load packages used in this analysis
#--------------------------------------------------
p <- c("optparse","ade4","pheatmap","lattice","ggplot2","squash","randomForest","pROC","reshape2","fmsb")
usePackage <- function(p) {
    if (!is.element(p, installed.packages()[,1]))
        install.packages(p, dep = TRUE, repos = "http://cran.us.r-project.org")
    suppressWarnings(suppressMessages(invisible(require(p, character.only = TRUE))))
}
invisible(lapply(p, usePackage))
#-------------------------------
# data format
#-------------------------------
# matrix: 
#         row.names	Sample_id
#         col.names	Varibles
# For example, data should be organized like this:
# Sample_id	group	V1	V2	etc...
# sample_0001	A	6	25
# sample_0002	B	9	32
# etc...
#------------------------------e

if(!is.na(as.logical(opts$prefix_of_prediction_output))){
outpath<-"./"
#dir.create(outpath<-paste("./",opts$prefix_of_training_output,"-",opts$prefix_of_prediction_output,".",opts$grouptype_of_training,"/",sep=""))
con <- file(paste(outpath,opts$prefix_of_training_output,"-",opts$prefix_of_prediction_output,".",opts$grouptype_of_training,".log",sep=""))}else{
outpath<-"./"
#dir.create(outpath<-paste("./",opts$prefix_of_training_output,".",opts$grouptype_of_training,"/",sep=""))
con <- file(paste(outpath,opts$prefix_of_training_output,".",opts$grouptype_of_training,".log",sep=""))
}

#sink(con, append=TRUE)
#sink(con, append=TRUE, type="message")

#-------------------------------
# Training data input
#-------------------------------
    g<-read.table(filename,header=T,row.names=1)
    g<-g[order(rownames(g)),]
    print(paste("The number of variables : ", ncol(g) ,sep=""))
    #-------------------------------filtering taxa with zero variance
    g<-g[,which(apply(g,2,var)!=0)]
    print(paste("The number of variables (removed variables with zero variance) : ", ncol(g) ,sep=""))
    gmat<-data.matrix(g)

#-------------------------------
# Metadata of training data input
#-------------------------------
    metadata<-read.table(metadata.filename,header=T,sep="\t",row.names=1)
    metadata<-metadata[order(rownames(metadata)),]
    if(length(metadata)==nrow(g)){group<-metadata; names(group)<-rownames(metadata)}else{group<-metadata[,opts$grouptype_of_training]; names(group)<-rownames(metadata)}
if(!is.na(as.logical(opts$prefix_of_prediction_output))){
#-------------------------------
# Test data input
#-------------------------------
    pred.g<-read.table(pred.filename,header=T,row.names=1)
    pred.g<-pred.g[order(rownames(pred.g)),]
    print(paste("The number of variables of prediction data: ", ncol(pred.g) ,sep=""))
    #-------------------------------filtering taxa with X% zero
    NonZero.p<-1
    pred.g<-pred.g[,which(colSums(pred.g==0)<NonZero.p*nrow(pred.g))]
    print(paste("The number of variables (removed variables containing over ", NonZero.p," zero) of prediction data: ", ncol(pred.g) ,sep=""))
    #-------------------------------
    pred.gmat<-data.matrix(pred.g)
#-------------------------------
# Metadata of test data input
#-------------------------------
    pred.metadata<-read.table(pred.metadata.filename,header=T,sep="\t",row.names=1)
    pred.metadata<-pred.metadata[order(rownames(pred.metadata)),]
    if(length(pred.metadata)==nrow(pred.g)){pred.group<-pred.metadata; names(pred.group)<-rownames(pred.metadata)}else{pred.group<-pred.metadata[,opts$grouptype_of_prediction]; names(pred.group)<-rownames(pred.metadata)}
#-------------------------------
# To get shared taxa between training and prediction data set
#-------------------------------
    gmatO<-gmat[,colnames(gmat) %in% colnames(pred.gmat)]
    pred.gmatO<-pred.gmat[,colnames(pred.gmat) %in% colnames(gmat)]
    print(paste("The number of variables for prediction in test data: ", ncol(gmatO) ,sep=""))
    #--------------------------------------------------
    x<-gmatO
    y<-group
    }else{
    x<-gmat
    y<-group
    }
#--------------------------------------------------
# Ten-folds cross validation in traing data
#--------------------------------------------------
    if(rf.opts$errortype != 'oob'){
    result<-rf.cross.validation(x,y,nfolds=rf.opts$nfolds, verbose=rf.opts$verbose,ntree=rf.opts$ntree)
    print(paste("The ", rf.opts$nfolds, "-folds cross validation for randomforest classification",sep=""))
    print(result$confusion.matrix)
    kappa.result.cv<-Kappa.test(result$confusion.matrix)
    acc.cv<-sum(diag(result$confusion.matrix))/sum(result$confusion.matrix)
    cat(paste("Cohen's kappa statistics for agreement in CV: ", kappa.result.cv$Result$estimate," -- ",kappa.result.cv$Judgement,"\n",sep=""))
    cat(paste("Accuracy in CV: ",acc.cv ,"\n",sep=""))
    
    #---------------------------   
    sink(paste(outpath,opts$prefix_of_training_output,".",opts$grouptype_of_training,".CV",rf.opts$nfolds,".Prob.xls",sep="")); cat("\t"); write.table(result$probabilities,sep="\t",quote=FALSE); sink()
    #----------------------------probability_boxplot
    pdf(paste(outpath,opts$prefix_of_training_output,".",opts$grouptype_of_training,".CV",rf.opts$nfolds,".probability_boxplot.pdf",sep=""), height=8, width=4)
    par(mar=c(4,4,3,3))
    boxplot(result$probabilities[,1]~y, ylab=paste("Probability of ",colnames(result$probabilities)[1],sep=""),ylim=c(0,1),outline=FALSE)
    abline(h=0.5,lty=2)
    stripchart(result$probabilities[,1]~y,method="jitter",jitter=.05,vertical=T,add=T,pch=20) 
    dev.off()
    
    if(nlevels(group)==2){
    sens.cv<-diag(result$confusion.matrix)[1]/sum(result$confusion.matrix[,1])
    spec.cv<-diag(result$confusion.matrix)[2]/sum(result$confusion.matrix[,2])
    cat(paste("Sensitivity in CV: ", sens.cv ,"\n",sep=""))
    cat(paste("Specificity in CV: ", spec.cv ,"\n",sep=""))
    #--------------------------------------------------
    #  ROC plot using "ROCR" package
    #--------------------------------------------------
    pred<-prediction(result$probabilities[,2],y)
    perf<-performance(pred,"tpr","fpr")
    auc.tmp <- performance(pred,"auc"); auc <- as.numeric(auc.tmp@y.values)
    cat(paste("AUC in CV: ", auc ,"\n",sep="")) 
    #----------------------------
    pdf(paste(outpath,opts$prefix_of_training_output,".CV",rf.opts$nfolds,".Prob.ROCR_plot.pdf",sep=""), height=4, width=4)
    par(mar=c(4,4,3,3))
    plot(perf,main="Training-CV",col=2,lwd=2)
    abline(a=0,b=1,lwd=2,lty=2,col="gray")
    text(0.8,0.2,paste("AUC: ", formatC(auc,digits=2,format="f"),sep=""))
    dev.off()
    #----------------------------
    pdf(paste(outpath,opts$prefix_of_training_output,".CV",rf.opts$nfolds,".ROCR.accuracy_VS_cutoff.plot.pdf",sep=""), height=4, width=4)
    par(mar=c(4,4,3,3))
    acc.perf <- performance(pred, "acc")
    plot(acc.perf, avg= "vertical")
    dev.off()
    #--------------------------------------------------
    #  ROC plot using "pROC" package
    #--------------------------------------------------
    pdf(paste(outpath,opts$prefix_of_training_output,".CV",rf.opts$nfolds,".Prob.pROC.ci.pdf",sep=""),width=5,height=5)
    rocobj <- plot.roc(y, result$probabilities[,2],main="", percent=TRUE,ci=TRUE) # print the AUC (will contain the CI)
    ciobj <- ci.se(rocobj,specificities=seq(0, 100, 5)) # over a select set of specificities
    plot(ciobj, type="shape", col="#1c61b6AA", print.auc=TRUE) # plot as a blue shape
    text(50,20,paste("AUC = ",formatC(rocobj$auc,digits=2,format="f"),sep=""),pos=4)
    ci.lower<-formatC(rocobj$ci[1],digits=2,format="f")
    ci.upper<-formatC(rocobj$ci[3],digits=2,format="f")
    text(50,10,paste("95% CI: ",ci.lower,"-",ci.upper,sep=""),pos=4)
    dev.off()
    }
    }
#--------------------------------------------------
# Out-of-bag classification in traing data and prediction in test data
#--------------------------------------------------
#--------------------------------------------------
# (1) Out-of-bag classification in traing data 
#--------------------------------------------------
    y
    oob.result <- rf.out.of.bag(x, y, verbose=rf.opts$verbose, ntree=rf.opts$ntree)
    print(paste("The out-of-bag (3-fold CV) result for randomforest classification",sep=""))
    #----------------------------
    print(oob.result$confusion.matrix)
    kappa.result.oob<-Kappa.test(oob.result$confusion.matrix)
    cat(paste("Cohen's kappa statistics for agreement in out-of-bag classification: ", kappa.result.oob$Result$estimate," -- ",kappa.result.oob$Judgement,"\n",sep=""))
    #---------------------------   
    sink(paste(outpath,opts$prefix_of_training_output,".",opts$grouptype_of_training,".oob.Prob.xls",sep="")); cat("\t"); write.table(oob.result$probabilities,sep="\t",quote=FALSE); sink()
    if(nlevels(group)==2){
    #--------------------------------------------------
    #  ROC plot using "pROC" package
    #--------------------------------------------------
    Cairo(file =paste(outpath,opts$prefix_of_training_output,".CV",rf.opts$nfolds,".Prob.pROC.ci.", device, sep=""), unit="in", dpi=300,  height = 7, width= 7, type=device, bg="white", units="in");
    #pdf(paste(outpath,opts$prefix_of_training_output,".CV",rf.opts$nfolds,".Prob.pROC.ci.pdf",sep=""),width=5,height=5)
    rocobj <- plot.roc(y, oob.result$probabilities[,2],main="", percent=TRUE,ci=TRUE) # print the AUC (will contain the CI)
    ciobj <- ci.se(rocobj,specificities=seq(0, 100, 5)) # over a select set of specificities
    plot(ciobj, type="shape", col="#1c61b6AA", print.auc=TRUE) # plot as a blue shape
    text(50,20,paste("AUC = ",formatC(rocobj$auc,digits=2,format="f"),sep=""),pos=4)
    ci.lower<-formatC(rocobj$ci[1],digits=2,format="f")
    ci.upper<-formatC(rocobj$ci[3],digits=2,format="f")
    text(50,10,paste("95% CI: ",ci.lower,"-",ci.upper,sep=""),pos=4)
    dev.off()
    }
#--------------------------------------------------
# (2) Prediction in test data
#--------------------------------------------------
    if(!is.na(as.logical(opts$prefix_of_prediction_output))){
    #----------------------------
    predicted.prob<-predict(oob.result$rf.model,pred.gmatO,type="prob")      
    sink(paste(outpath,opts$prefix_of_prediction_output,".",opts$grouptype_of_prediction,".oob.Pred_Prob.xls",sep="")); cat("\t"); write.table(predicted.prob,sep="\t",quote=FALSE); sink()
    #----------------------------probability_boxplot
    pdf(paste(outpath,opts$prefix_of_prediction_output,".",opts$grouptype_of_prediction,".oob.prob_boxplot.pdf",sep=""), height=8, width=4)
    par(mar=c(4,4,3,3))
    boxplot(predicted.prob[,1]~pred.group,xlab=opts$prefix_of_prediction_output, ylab=paste("Probability of ",colnames(predicted.prob)[1],sep=""),ylim=c(0,1),outline=FALSE)
    abline(h=0.5,lty=2)
    stripchart(predicted.prob[,1]~pred.group,method="jitter",jitter=.05,vertical=T,add=T,pch=20) 
    dev.off()
    
    if(nlevels(group)==2){
    #----------------------------probability_xyplot
    pdf(paste(outpath,opts$prefix_of_prediction_output,".",opts$grouptype_of_prediction,".oob.prob_xyplot.pdf",sep=""), height=6, width=6)
    par(mar=c(4,4,3,3))
    plot(predicted.prob[,1],predicted.prob[,2],xlab=paste("Probability of ",colnames(predicted.prob)[1],sep=""),
                                               ylab=paste("Probability of ",colnames(predicted.prob)[2],sep=""),
                                               main=opts$prefix_of_prediction_output, xlim=c(0,1),ylim=c(0,1),col="white",pch=21,bg="black")
    abline(v=0.5,h=0.5,lty=2)
    dev.off()
    #----------------------------
    if(nlevels(pred.group)==2){
    #--------------------------------------------------
    #  ROC plot using "ROCR" package
    #--------------------------------------------------
    pred<-prediction(predicted.prob[,2],pred.group)
    perf<-performance(pred,"tpr","fpr")
    auc.tmp <- performance(pred,"auc"); auc <- as.numeric(auc.tmp@y.values)
    #----------------------------
    pdf(paste(outpath,opts$prefix_of_prediction_output,".",opts$grouptype_of_prediction,".oob.ROCR_plot.pdf",sep=""), height=4, width=4)
    par(mar=c(4,4,3,3))
    plot(perf,main="Predicted",col=2,lwd=2)
    abline(a=0,b=1,lwd=2,lty=2,col="gray")
    text(0.8,0.2,paste("AUC: ", formatC(auc,digits=2,format="f"),sep=""))
    dev.off()
    #----------------------------
    pdf(paste(outpath,opts$prefix_of_training_output,".",opts$grouptype_of_prediction,".oob.ROCR.accuracy_VS_cutoff.plot.pdf",sep=""), height=4, width=4)
    par(mar=c(4,4,3,3))
    acc.perf <- performance(pred, "acc")
    plot(acc.perf, avg= "vertical")
    dev.off()
    }
    }
    #--------------------------------------------------
    #  probability_boxplot of both train and test data
    #--------------------------------------------------
    if(opts$grouptype_of_prediction==opts$grouptype_of_prediction){
    comb.group<-as.factor(c(as.character(group),as.character(pred.group)))
    #comb.group<-factor(comb.group,levels=c(levels(comb.group)[2],levels(comb.group)[3],levels(comb.group)[1]))
    comb.prob<-rbind(oob.result$probabilities,predicted.prob)
    pdf(paste(outpath,opts$prefix_of_training_output,"-",opts$prefix_of_prediction_output,".",opts$grouptype_of_training,".probability_boxplot.pdf",sep=""), height=8, width=4)
    par(mar=c(4,4,3,3))
    boxplot(comb.prob[,1]~comb.group, ylab=paste("Probability of ",colnames(comb.prob)[1],sep=""),ylim=c(0,1),outline=FALSE)
    abline(h=0.5,lty=2)
    stripchart(comb.prob[,1]~comb.group,method="jitter",jitter=.05,vertical=T,add=T,pch=20) 
    dev.off()
    
    sink(paste(outpath,opts$prefix_of_training_output,"-",opts$prefix_of_prediction_output,".",opts$grouptype_of_training,".Prob.xls",sep="")); cat("\t"); write.table(comb.prob,sep="\t",quote=FALSE); sink()
    }
    }
    
#--------------------------------------------------
#   RF importances of features
#--------------------------------------------------
    imps<-oob.result$importances
    imps<-imps[order(imps,decreasing = TRUE)]
    imps.cutoff<-imps[which(imps>0)]
    sink(paste(outpath,opts$prefix_of_training_output,".",opts$grouptype_of_training,".imps_All_sorted.xls",sep=""));cat("\t");write.table(imps,quote=FALSE,sep="\t");sink()
    sink(paste(outpath,opts$prefix_of_training_output,".",opts$grouptype_of_training,".imps_cutoff_sorted.xls",sep=""));cat("\t");write.table(imps.cutoff,quote=FALSE,sep="\t");sink()
	
#--------------------------------------------------
#   Univariate BetweenGroup test
#--------------------------------------------------
    test.results<-BetweenGroup.test(x,y,p.adj.method = "fdr")
    cat("hereeeeeeee\n")
    if(nlevels(group)==2){
    sink(paste(outpath,opts$prefix_of_training_output,".",opts$grouptype_of_training,".imps_wilcox.",nrow(test.results),".xls",sep=""))
    }else{
    sink(paste(outpath,opts$prefix_of_training_output,".",opts$grouptype_of_training,".imps_kruscal.",nrow(test.results),".xls",sep=""))}
    cat("\t");write.table(test.results,quote=FALSE,sep="\t");sink()
    logP=-log(test.results[,ncol(test.results)-2],10)
    logP_cutoff=-log(alpha, 10)
#--------------------------------------------------
#   univariate AUC
#--------------------------------------------------    
    AUCs<-apply(x,2,function(x) auc(y,x))
    
#--------------------------------------------------
#   Comparisons of feature importance
#--------------------------------------------------        
    ImpP<-data.frame(RF_imp=oob.result$importances,Uni_logP=logP, Uni_AUC=AUCs)
    sink(paste(outpath,opts$prefix_of_training_output,".",opts$grouptype_of_training,".imps_comparisons.xls",sep=""));cat("\t");write.table(ImpP,quote=FALSE,sep="\t");sink()
#--------------------------------------------------
#   Mean spectral Plot 
#--------------------------------------------------     	
	N<-100
    imps.IfTopN<-factor(rank(-oob.result$importances)<N)
    imps.IflogP_cutoff<-factor(logP>logP_cutoff); sigN<-summary(imps.IflogP_cutoff)[2]
    meanResults<-data.frame(Shifts=as.numeric(gsub("X","",rownames(test.results))),test.results[1:nlevels(y)],Imps=oob.result$importances, Imps.IfTopN=imps.IfTopN, Imps.IflogP_cutoff=imps.IflogP_cutoff)
    meanSpec<-melt(meanResults,id.var=c("Shifts","Imps","Imps.IfTopN", "Imps.IflogP_cutoff"),variable.name="Group")
    p<-ggplot(meanSpec,aes(x=Shifts,y=value,group=Group))+
        geom_line(aes(color=Group), alpha=0.8)+theme_bw()+
        xlab("Raman shift (cm-1)")+ylab("Intensity")+
        geom_point(data=subset(meanSpec,Imps.IflogP_cutoff==TRUE & !duplicated(Shifts)),aes(x=Shifts, alpha=Imps))+
        #geom_vline(data=subset(meanSpec,Imps.IflogP_cutoff==TRUE & !duplicated(Shifts)),aes(xintercept =Shifts, alpha=Imps))+
        geom_text(data=subset(meanSpec,Imps.IflogP_cutoff==TRUE & !duplicated(Shifts)),aes(x=Shifts,y=value,label =Shifts),angle = 90,hjust=1, vjust=1,size=0.5,show_guide  = F)
    ggsave(filename=paste(outpath,opts$prefix_of_training_output,".Sig",sigN,".SpectrumPlot.ggplot.", device, sep=""),plot=p, width=10, height=4, device=device)
    #ggsave(filename=paste(outpath,opts$prefix_of_training_output,".Top",N,".SpectrumPlot.ggplot.pdf",sep=""),plot=p,width=12, height=4)
   
#--------------------------------------------------
if(rf.opts$rf.cv){
#--------------------------------------------------
#--------------------------------------------------
# Estimate the minErr of RF model and the top discriminatory taxa
#--------------------------------------------------
    results <- replicate(10, rfcv(x, y, step=0.9,cv.fold=10), simplify=FALSE)
    err.cv <- sapply(results, "[[", "error.cv")
    #--------------------------------------------------matplot
    pdf(paste(outpath,opts$prefix_of_training_output,".CV.MeanErrPlot.pdf",sep=""), height=6, width=6)
    matplot(results[[1]]$n.var, cbind(rowMeans(err.cv), err.cv), type="p",log="x",
            col=c(2, rep("grey60", ncol(err.cv))), 
            pch=c(19, rep(20, ncol(err.cv))), 
            xlab="Number of variables", ylab="CV Error")
    lines(results[[1]]$n.var,rowMeans(err.cv),col=2)
    breaks<-axTicks(side=1)
    dev.off()
    #--------------------------------------------------ggplot
    ErrSumm<-data.frame(NumVars=results[[1]]$n.var,MeanErr=apply(err.cv,1,mean),MedianErr=apply(err.cv,1,median),SdErr=apply(err.cv,1,sd),SeErr=apply(err.cv,1,function(x) sd(x)/sqrt(length(x))))
    #-------------------------------------------------MeanErr
    p<-ggplot(ErrSumm,aes(x=NumVars,y=MeanErr))+
       geom_point(size=2)+geom_line(alpha=0.8)+
       xlab("Number of peaks")+ylab("Ten-fold CV error")+
       geom_errorbar(aes(ymin=MeanErr-SeErr,ymax=MeanErr+SeErr),width=0.05)+
       coord_flip()+
       scale_x_continuous(trans = "log",breaks=breaks)+
       geom_hline(aes(yintercept=min(MeanErr)),linetype="longdash",alpha=0.2)
    ggsave(filename=paste(outpath,opts$prefix_of_training_output,".CV.MeanErrPlot.ggplot.pdf",sep=""),plot=p,width=4,height=6)
    
    #-------------------------------------------------
    minErr_n_imps<-names(which.min(rowMeans(err.cv)))
    imps.minErr<-imps[1:minErr_n_imps]
    #-------------------------------------------------
    sink(paste(outpath,opts$prefix_of_training_output,".",opts$grouptype_of_training,".imps_minErr_sorted.xls",sep=""));cat("\t");write.table(imps.minErr,quote=FALSE,sep="\t");sink()
    
}


#-------------------------------
sink(type="message")
sink() 
#-------------------------------


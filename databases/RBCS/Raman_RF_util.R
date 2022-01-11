
# install necessary libraries
p <- c("optparse", "RColorBrewer", "permute", "ggplot2", "reshape2", "dplyr", "Cairo")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://mirrors.opencas.cn/cran/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))


# runs cross-validation 
# if predict.fun is NULL, uses S3 predict method
# if nfolds > length(y) or nfolds==-1, uses leave-one-out cross-validation
# ...: additional parameters for train.fun
#
# value:
# y: true values
# predicted: cv predicted values
# probabilities: cv predicted class probabilities (or NULL if unavailable)
# confusion.matrix: confusion matrix (true x predicted)
# nfolds: nfolds
# params: list of additional parameters
# importances: importances of features as predictors

# Get probability of mislabeling by several measures
# returns matrix of p(alleged), max(p(others)), p(alleged) - max(p(others))
"get.mislabel.scores" <- function(y,y.prob){
    result <- matrix(0,nrow=length(y),ncol=3)
    # get matrices containing only p(other classes), and containing only p(class)
    mm <- model.matrix(~0 + y)
    y.prob.other.max <- apply(y.prob * (1-mm),1,max)
    y.prob.alleged <- apply(y.prob * mm, 1, max)
    result <- cbind(y.prob.alleged, y.prob.other.max, y.prob.alleged - y.prob.other.max)
    rownames(result) <- rownames(y.prob)
    colnames(result) <- c('P(alleged label)','P(second best)','P(alleged label)-P(second best)')
    return(result)
}

# Get balanced folds where each fold has close to overall class ratio
"balanced.folds" <- function(y, nfolds=10){
    folds = rep(0, length(y))
    classes = levels(y)
    # size of each class
    Nk = table(y)
    # -1 or nfolds = len(y) means leave-one-out
    if (nfolds == -1 || nfolds == length(y)){
        invisible(1:length(y))
    }
    else{
    # Can't have more folds than there are items per class
    nfolds = min(nfolds, max(Nk))
    # Assign folds evenly within each class, then shuffle within each class
        for (k in 1:length(classes)){
            ixs <- which(y==classes[k])
            folds_k <- rep(1:nfolds, ceiling(length(ixs) / nfolds))
            folds_k <- folds_k[1:length(ixs)]
            folds_k <- sample(folds_k)
            folds[ixs] = folds_k
        }
        invisible(folds)
    }
}


"rf.cross.validation" <- function(x, y, nfolds=10, verbose=verbose, ...){
    if(nfolds==-1) nfolds <- length(y)
    folds <- balanced.folds(y,nfolds=nfolds)
    result <- list()
    result$y <- as.factor(y)
    result$predicted <- result$y
    result$probabilities <- matrix(0, nrow=length(result$y), ncol=length(levels(result$y)))
    rownames(result$probabilities) <- rownames(x)
    colnames(result$probabilities) <- levels(result$y)
    result$importances <- matrix(0,nrow=ncol(x),ncol=nfolds)
    result$errs <- numeric(length(unique(folds)))

    # K-fold cross-validation
    for(fold in sort(unique(folds))){
        if(verbose) cat(sprintf('Fold %d...\n',fold))
        foldix <- which(folds==fold)
        model <- randomForest(x[-foldix,], factor(result$y[-foldix]), importance=TRUE, do.trace=verbose, ...)
        newx <- x[foldix,]
        if(length(foldix)==1) newx <- matrix(newx,nrow=1)
        result$predicted[foldix] <- predict(model, newx)
        probs <- predict(model, newx, type='prob')
        result$probabilities[foldix,colnames(probs)] <- probs
        result$errs[fold] <- mean(result$predicted[foldix] != result$y[foldix])
        result$importances[,fold] <- model$importance[,'MeanDecreaseAccuracy']
    }

    result$nfolds <- nfolds
    result$params <- list(...)
    result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
    return(result)
}

# Runs standard random forests with out-of-bag error estimation
# This is merely a wrapper that extracts relevant info
# Return values are the same as rf.cross.validation
"rf.out.of.bag" <- function(x,y, verbose=verbose, ...){
    rf.model <- randomForest(x,y,keep.inbag=TRUE,importance=TRUE,do.trace=verbose,...)
    result <- list()
    result$rf.model <- rf.model
    result$probabilities <- get.oob.probability.from.forest(rf.model,x)
    result$y <- y
    result$predicted <- rf.model$predicted
    result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
    result$params <- list(ntree=rf.opts$ntree)
    result$errs <- as.numeric(result$predicted != result$y)
    result$importances <- rf.model$importance[,'MeanDecreaseAccuracy']
    return(result)
}

# get probability of each class using only out-of-bag predictions from RF
"get.oob.probability.from.forest" <- function(model,x){
    # get aggregated class votes for each sample using only OOB trees
    votes <- get.oob.votes.from.forest(model,x)
    # convert to probs
    probs <- sweep(votes, 1, apply(votes, 1, sum), '/')
    rownames(probs) <- rownames(x)
    colnames(probs) <- model$classes
    
    return(invisible(probs))
}

# get votes for each class using only out-of-bag predictions from RF
"get.oob.votes.from.forest" <- function(model,x){
    # get aggregated class votes for each sample using only OOB trees
    votes <- matrix(0, nrow=nrow(x), ncol=length(model$classes))
   
    rf.pred <- predict(model, x, type="vote",predict.all=T)
    for(i in 1:nrow(x)){
        # find which trees are not inbag for this sample
        outofbag <- model$inbag[i,]==0
        # get oob predictions for this sample
        votes[i,] <- table(factor(rf.pred$individual[i,][outofbag],levels=model$classes))
    }
    rownames(votes) <- rownames(x)
    colnames(votes) <- model$classes
    
    return(invisible(votes))
}

# prints random forests results file
"save.rf.results" <- function(result, rf.opts, feature.ids, outdir){
    save.rf.results.summary(result, rf.opts, outdir=outdir)
    save.rf.results.probabilities(result, outdir=outdir)
    save.rf.results.mislabeling(result, outdir=outdir)
    save.rf.results.importances(result, feature.ids=feature.ids, outdir=outdir)
    save.rf.results.confusion.matrix(result, outdir=outdir)
}

# Print "summary" file
"save.rf.results.summary" <- function(result, rf.opts, filename='summary.xls', outdir){
    err <- mean(result$errs)
    err.sd <- sd(result$errs)
    baseline.err <- 1-max(table(y))/length(y)
    filepath <- sprintf('%s/%s',outdir, filename)
    sink(filepath)
    cat(sprintf('Model\tRandom Forest\n'))
    cat(sprintf('Error type\t%s\n',result$error.type))
    if(rf.opts$errortype == 'oob' || rf.opts$errortype == 'cvloo'){
        cat(sprintf('Estimated error\t%.5f\n',err))
    } else {
        cat(sprintf('Estimated error (mean +/- s.d.)\t%.5f +/- %.5f\n',err,err.sd))
    }
    cat(sprintf('Baseline error (for random guessing)\t%.5f\n',baseline.err))
    cat(sprintf('Ratio baseline error to observed error\t%.5f\n',baseline.err / err))
    cat(sprintf('Number of trees\t%d\n',result$params$ntree))
    sink(NULL)
}

# Print "probabilities" file
"save.rf.results.probabilities" <- function(result, filename='cv_probabilities.xls', outdir){
    filepath <- sprintf('%s/%s',outdir,filename)
    sink(filepath)
    cat('SampleID\t')
    write.table(result$probabilities,sep='\t',quote=F)
    sink(NULL)
}

# Print "mislabeling" file
"save.rf.results.mislabeling" <- function(result, filename='mislabeling.xls', outdir){
    filepath <- sprintf('%s/%s',outdir,filename)
    sink(filepath)
    cat('SampleID\t')
    write.table(get.mislabel.scores(result$y,result$probabilities),sep='\t',quote=F)
    sink(NULL)
}

"save.rf.results.importances" <- function(result,feature.ids, filename='feature_importance_scores.xls', outdir){
    filepath <- sprintf('%s/%s',outdir,filename)
    if(is.null(dim(result$importances))){
        imp <- result$importances
        imp.sd <- rep(NA,length(imp))
    } else {
        imp <- rowMeans(result$importances)
        imp.sd <- apply(result$importances, 1, sd)
    }
    output.table <- cbind(imp, imp.sd)
    rownames(output.table) <- feature.ids
    output.table <- output.table[sort(imp,dec=T,index=T)$ix,]
    colnames(output.table) <- c('Mean_decrease_in_accuracy','Standard_deviation')

    sink(filepath)
    cat('Feature_id\t')
    write.table(output.table,sep='\t',quote=F)
    sink(NULL)
}

# Print "confusion matrix" file
"save.rf.results.confusion.matrix" <- function(result, filename='confusion_matrix.xls', outdir){
    filepath <- sprintf('%s/%s',outdir,filename)
    
    # add class error column to each row
    x <- result$confusion.matrix
    class.errors <- rowSums(x * (1-diag(nrow(x)))) / rowSums(x)
    output <- cbind(result$confusion.matrix, class.errors)
    colnames(output)[ncol(output)] <- "Class error"
    sink(filepath)
    cat('True\\Predicted\t')
    write.table(output,quote=F,sep='\t')
    sink(NULL)
}

BetweenGroup.test <-function(data, group, p.adj.method="bonferroni",paired=FALSE){
# p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")

n_group<-nlevels(group)
        if(!is.numeric(n_group) | n_group==1)
        stop("group must be a numeric and up to two levels\n")
if(n_group==2){
               output1<-matrix(NA,ncol=9,nrow=ncol(data))
               rownames(output1)<-colnames(data)
               colnames(output1)<-c(paste("mean_",levels(group)[1],sep=""),paste("mean_",levels(group)[2],sep=""),
                                   paste("sd_",levels(group)[1],sep=""),paste("sd_",levels(group)[2],sep=""),"Var.test","T-test","Wilcoxon.test",
                                   paste("T-test_",p.adj.method,sep=""),paste("Wilcoxon.test_",p.adj.method,sep=""))
               for(i in 1:ncol(data))
               {
               output1[i,1]<-mean(data[which(group==levels(group)[1]),i])
               output1[i,2]<-mean(data[which(group==levels(group)[2]),i])
               output1[i,3]<-sd(data[which(group==levels(group)[1]),i])
               output1[i,4]<-sd(data[which(group==levels(group)[2]),i])
               output1[i,5]<-var.test(data[,i]~group)$p.value
               if(output1[i,5]<0.01)
               output1[i,6]<-t.test(data[,i]~group,paired=paired)$p.value
               else
               output1[i,6]<-t.test(data[,i]~group, var.equal=T,paired=paired)$p.value
               output1[i,7]<-wilcox.test(data[,i]~group, paired=paired, conf.int=TRUE, exact=FALSE, correct=FALSE)$p.value
               output1[i,8]<-NA
               output1[i,9]<-NA
               }
               output1[,8]<-p.adjust(output1[,6], method = p.adj.method, n = ncol(data))
               output1[,9]<-p.adjust(output1[,7], method = p.adj.method, n = ncol(data))
               
               return(data.frame(output1))
}else{
      output2<-matrix(NA,ncol=n_group+5,nrow=ncol(data))
      rownames(output2)<-colnames(data)
      colnames.output2<-array(NA)
      for(j in 1:ncol(output2)){
      if(j<=n_group){
      colnames.output2[j]<-c(paste("mean_",levels(group)[j],sep=""))
      }else{
      colnames.output2[(n_group+1):(n_group+5)]<-c("Var.test","Oneway-test","Kruskal.test",
                                                    paste("Oneway-test_",p.adj.method,sep=""),paste("Kruskal.test_",p.adj.method,sep=""))
                                                    }
      }
      colnames(output2)<-colnames.output2
      for(i in 1:ncol(data))
      {
      for(j in 1:n_group)
      {
      output2[i,j]<-mean(data[which(group==levels(group)[j]),i])
      }
      output2[i,(n_group+1)]<-bartlett.test(data[,i]~group)$p.value
      if(output2[i,(n_group+1)]<0.01)
      output2[i,(n_group+2)]<-oneway.test(data[,i]~group)$p.value
      else
      output2[i,(n_group+2)]<-oneway.test(data[,i]~group, var.equal=T)$p.value
      output2[i,(n_group+3)]<-kruskal.test(data[,i]~group)$p.value
      output2[i,(n_group+4)]<-NA
      output2[i,(n_group+5)]<-NA
      }
      output2[ ,(n_group+4)]<-p.adjust(output2[,(n_group+2)], method = p.adj.method, n = ncol(data))
      output2[ ,(n_group+5)]<-p.adjust(output2[,(n_group+3)], method = p.adj.method, n = ncol(data))
      return(data.frame(output2))
      }
      
      
}


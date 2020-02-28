MergeLists=function(mylist){
  return(mapply(c, mylist[[1]], mylist[[2]]))
}


# dataX=proteins; dataY=as.numeric(covars$case); ncomp=1; Nrepeat=niter; xseq=myxseq

CalibratesPLSDA=function(dataX, dataY, ncomp=1, Nrepeat=100, name="1", xseq=NULL){
  # only for univ Y
  TmpSummary <- NULL
  SelectedX<-NULL
  comp_error=NULL
  
  if (is.null(xseq)){
    xseq=1:ncol(dataX)
  }
  all_ids=1:length(dataY)
  for (comp in 1:ncomp){
    for(i in xseq){
      TmpKeepX <- c(SelectedX,i)
      TmpsPLS <- splsda(dataX,dataY,keepX=TmpKeepX,ncomp=comp,mode='regression')
      error=NULL
      print(i)
      pb=txtProgressBar(style=3)
      for (k in 1:Nrepeat){
        folds=MergeLists(lapply(split(1:length(dataY), f=dataY), FUN=function(x){split(x, rep(1:5, length.out=length(x)))}))
        for (l in 1:length(folds)){
          test_ids=folds[[l]]
          train_ids=all_ids[!all_ids%in%test_ids]
          dataX_train=dataX[train_ids,]
          dataY_train=dataY[train_ids]
          dataX_test=dataX[test_ids,]
          dataY_test=dataY[test_ids]
          
          sPLS_train <- splsda(dataX_train,dataY_train,keepX=TmpKeepX,ncomp=comp,mode='regression')
          predicted=predict(sPLS_train, newdata=dataX_test, dist="max.dist")
          mytable=table(dataY_test, predicted$class$max.dist)
          misclassif_rate=1-sum(diag(mytable))/sum(mytable)
          error=c(error, misclassif_rate)
        }
        # TmpPerf=perf.splsda(TmpsPLS, validation = "Mfold", dist = "max.dist",
        #                     # folds=lapply(split(sample(seq(1, nrow(dataX), by=2)), rep(1:5, length=length(dataY)/2)), FUN=function(x){return(c(x, (x+1)))}),
        #                     folds=MergeLists(lapply(split(1:length(dataY), f=dataY), FUN=function(x){split(x, rep(1:5, length.out=length(x)))})),
        #                     nrepeat=1, progressBar = FALSE)
        # error=c(error, TmpPerf$error.rate$overall)
        setTxtProgressBar(pb, k/Nrepeat)
      }
      cat("\n")
      
      TmpSummary=rbind(TmpSummary, c(comp, i, mean(error)))
      TmpSummary=as.data.frame(TmpSummary, stringsAsFactors = FALSE)
      colnames(TmpSummary)=c("Ncomp", "NVar", "error")
    }
    SelectedX=c(SelectedX, which.min(TmpSummary$error[TmpSummary$Ncomp==comp]))
    comp_error=c(comp_error, min(TmpSummary$error[TmpSummary$Ncomp==comp]))
  }
  
  if (ncomp>1){
    names(comp_error)=paste0("Comp", seq(1,ncomp))
    res=list(Summary=TmpSummary, CompMinError=comp_error, 
             NComp=which.min(comp_error), NVar=SelectedX[1:which.min(comp_error)])
  } else {
    res=list(Summary=TmpSummary, MinError=comp_error, NVar=SelectedX[1:which.min(comp_error)])
  }
  
  return(res)
}


CalibrategPLSDA=function(dataX, dataY, ncomp=1, Nrepeat=100, name="1", Xgroups){
  TmpSummary=NULL
  comp_error=NULL
  SelectedX=NULL
  
  all_ids=1:length(dataY)
  for (comp in 1:ncomp){
    for (NGroups in 1:(length(Xgroups)+1)){
      error=NULL
      TmpGroups=c(SelectedX, NGroups)
      pb=txtProgressBar(style=3)
      for (i in 1:Nrepeat){
        folds=MergeLists(lapply(split(1:length(dataY), f=dataY), FUN=function(x){split(x, rep(1:5, length.out=length(x)))}))
        for (l in 1:length(folds)){
          test_ids=folds[[l]]
          train_ids=all_ids[!all_ids%in%test_ids]
          dataX_train=dataX[train_ids,]
          dataY_train=dataY[train_ids]
          dataX_test=dataX[test_ids,]
          dataY_test=dataY[test_ids]
          
          gPLS_train <- gPLSda(dataX_train, dataY_train, ncomp=comp, ind.block.x=Xgroups, keepX=TmpGroups)
          predicted=predict(gPLS_train, newdata=dataX_test, dist="max.dist")
          mytable=table(dataY_test, predicted$class$max.dist)
          misclassif_rate=1-sum(diag(mytable))/sum(mytable)
          error=c(error, misclassif_rate)
        }
        setTxtProgressBar(pb, i/Nrepeat)
      }
      print(NGroups)
      TmpSummary=rbind(TmpSummary, c(comp, NGroups, mean(error)))
      TmpSummary=as.data.frame(TmpSummary, stringsAsFactors = FALSE)
      colnames(TmpSummary)=c("Ncomp", "NGroups", "error")
    }
    SelectedX=c(SelectedX, which.min(TmpSummary$error[TmpSummary$Ncomp==comp]))
    comp_error=c(comp_error, min(TmpSummary$error[TmpSummary$Ncomp==comp]))
  }
  
  
  if (ncomp>1){
    names(comp_error)=paste0("Comp", seq(1,ncomp))
    res=list(Summary=TmpSummary, CompMinError=comp_error,
             NComp=which.min(comp_error), NGroup=SelectedX[1:which.min(comp_error)])
  } else {
    res=list(Summary=TmpSummary, MinError=comp_error, NGroup=SelectedX[1:which.min(comp_error)])
  }
  return(res)
}


CalibratesgPLSDA=function(dataX, dataY, ncomp=1, Nrepeat=100, name="1", Xgroups, Alpha=seq(0.9, 0, by=-0.1)){
  TmpSummary=NULL
  SelectedGroups=NULL
  SelectedAlpha=NULL
  comp_error=NULL
  
  all_ids=1:length(dataY)
  for (comp in 1:ncomp){
    for (NGroups in 1:(length(Xgroups)+1)){
      print(NGroups)
      for (alpha in Alpha){
        error=NULL
        TmpGroups=c(SelectedGroups, NGroups)
        TmpAlpha=c(SelectedAlpha, alpha)
        TmpsgPLSDA <- sgPLSda(dataX, dataY, ncomp = comp, ind.block.x = Xgroups, keepX = TmpGroups, alpha.x = TmpAlpha)
        pb=txtProgressBar(style=3)
        for (i in 1:Nrepeat){
          folds=MergeLists(lapply(split(1:length(dataY), f=dataY), FUN=function(x){split(x, rep(1:5, length.out=length(x)))}))
          for (l in 1:length(folds)){
            test_ids=folds[[l]]
            train_ids=all_ids[!all_ids%in%test_ids]
            dataX_train=dataX[train_ids,]
            dataY_train=dataY[train_ids]
            dataX_test=dataX[test_ids,]
            dataY_test=dataY[test_ids]
            
            sgPLS_train <- sgPLSda(dataX_train, dataY_train, ncomp=comp, 
                                  ind.block.x=Xgroups, keepX=TmpGroups, alpha.x=TmpAlpha)
            predicted=predict(sgPLS_train, newdata=dataX_test, dist="max.dist")
            mytable=table(dataY_test, predicted$class$max.dist)
            misclassif_rate=1-sum(diag(mytable))/sum(mytable)
            error=c(error, misclassif_rate)
          }
          
          setTxtProgressBar(pb, i/Nrepeat)
        }
        print(alpha)
        TmpSummary=rbind(TmpSummary, c(comp, NGroups, alpha, mean(error)))
        TmpSummary=as.data.frame(TmpSummary, stringsAsFactors = FALSE)
        colnames(TmpSummary)=c("Ncomp", "NGroups", "alpha", "error")
      }
    }
    SelectedGroups=c(SelectedGroups, TmpSummary$NGroups[TmpSummary$Ncomp==comp][which.min(TmpSummary$error[TmpSummary$Ncomp==comp])])
    SelectedAlpha=c(SelectedAlpha, TmpSummary$alpha[TmpSummary$Ncomp==comp][which.min(TmpSummary$error[TmpSummary$Ncomp==comp])])
    comp_error=c(comp_error, min(TmpSummary$error[TmpSummary$Ncomp==comp]))
  }
  
  if (ncomp>1){
    names(comp_error)=paste0("Comp", seq(1,ncomp))
    res=list(Summary=TmpSummary, CompMinError=comp_error, 
             NComp=which.min(comp_error), 
             NGroup=SelectedGroups[1:which.min(comp_error)], alpha=SelectedAlpha[1:which.min(comp_error)])
  } else {
    res=list(Summary=TmpSummary, MinError=comp_error, 
             NGroup=TmpSummary$NGroups[which.min(TmpSummary$error)], alpha=TmpSummary$alpha[which.min(TmpSummary$error)])
  }
  return(res)
}


PlotCalib=function(res, type='sPLSDA', ncomp_selected=1, main=NULL){
  if (type=='gPLSDA'){
    name='groups'
    res$Summary$NVar=res$Summary$NGroups
  } else {
    name='variables'
  }
  
  if (type=='sgPLSDA'){
    par(mar=c(5,7,3,1))
    plot(1:sum(res$Summary$Ncomp==ncomp_selected),
         res$Summary$error[res$Summary$Ncomp==ncomp_selected],
         type="l", xaxt="n", xlab="", ylab="Misclassification rate",
         main="sgPLS-DA calibration", cex.lab=1.5)
    axis(side=1, at=1:sum(res$Summary$Ncomp==ncomp_selected),
         labels=res$Summary$alpha[res$Summary$Ncomp==ncomp_selected])
    tmp=c(which(!duplicated(res$Summary$NGroups[res$Summary$Ncomp==ncomp_selected]))-0.5, sum(res$Summary$Ncomp==ncomp_selected)+0.5)
    abline(v=tmp, lty=2, col="grey")
    axis(side=1, at=tmp, labels=NA, line=2.5)
    axis(side=1, at=apply(rbind(tmp[-1], tmp[-length(tmp)]),2,mean),
         labels=unique(res$Summary$NGroups), line=2.5, tick=FALSE)
    points(which.min(res$Summary$error[res$Summary$Ncomp==ncomp_selected]),
           res$Summary$error[res$Summary$Ncomp==ncomp_selected][which.min(res$Summary$error[res$Summary$Ncomp==ncomp_selected])],
           pch=19, col="red")
    mtext("Penalty", side=1, line=1, at=0, adj=1)
    mtext("Number of groups", side=1, line=3.5, at=0, adj=1)
  } else {
    par(mar=c(5,5,3,1))
    plot(res$Summary$NVar[res$Summary$Ncomp==ncomp_selected],
         res$Summary$error[res$Summary$Ncomp==ncomp_selected],
         type="l", xaxt="n", xlab=paste0("Number of ", name), ylab="Misclassification rate",
         main=ifelse(is.null(main), yes=paste0("Calibration of the number of ", name), no=main),
         sub=paste("Number of components =", ncomp_selected), cex.lab=1.5)
    # axis(side = 1, at = 1:nrow(res$Summary),
    #labels=paste(res$Summary$NGroups, "-", res$Summary$alpha))
    points(which.min(res$Summary$error[res$Summary$Ncomp==ncomp_selected]),
           res$Summary$error[res$Summary$Ncomp==ncomp_selected][which.min(res$Summary$error[res$Summary$Ncomp==ncomp_selected])],
           pch=19, col="red")
    points(res$Summary$NVar[res$Summary$Ncomp==ncomp_selected],
           res$Summary$error[res$Summary$Ncomp==ncomp_selected], pch=19, cex=0.5)
    axis(side = 1, at = res$Summary$NVar[res$Summary$Ncomp==ncomp_selected])
    points(res$Summary$NVar[which.min(res$Summary$error[res$Summary$Ncomp==ncomp_selected])],
           res$Summary$error[res$Summary$Ncomp==ncomp_selected][which.min(res$Summary$error[res$Summary$Ncomp==ncomp_selected])],
           pch=19, col="red")
  }
}


StabilityPlot=function(X, Y, NIter=100, plot=TRUE){
  MyStab=NULL
  
  pb=txtProgressBar(style=3)
  for (NVar in seq((ncol(X)-1), 1, by=-1)){
    setTxtProgressBar(pb, 1-(NVar-1)/(ncol(X)))
    TmpStab=NULL
    for (k in 1:NIter){
      s=sample(seq(1, nrow(X)), size=0.8*nrow(X), replace = FALSE) # subsampling procedure
      X_sub=X[s,]
      rownames(X_sub)=s
      Y_sub=Y[s]
      
      TmpsPLS=splsda(X_sub, Y_sub, keepX=NVar, ncomp=1, mode='regression')
      TmpStab=rbind(TmpStab, TmpsPLS$loadings$X[,1])
    }
    MyStab=rbind(MyStab, (apply(TmpStab, 2, FUN=function(x){sum(x!=0)})/NIter))
  }
  
  MyStab=rbind(rep(1, ncol(X)), MyStab)
  rownames(MyStab)=seq(ncol(X), 1, by=-1)
  
  if (plot){
    pheatmap(MyStab, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE)
  }
  
  return(MyStab)
}



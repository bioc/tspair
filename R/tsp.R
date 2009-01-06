tspcalc <- function(dat,grp){
  if(class(dat) == "ExpressionSet"){
    genenames <- as.character(1:dim(exprs(dat))[1])
    if(!is.null(featureNames(dat))){genenames <- featureNames(dat)}
    if(length(grp) == 1){
      labels <- as.character(unique(pData(dat)[,grp]))
      grp <- make.consecutive.int(pData(dat)[,grp])
      dat <- exprs(dat)
      rownames(dat) <- genenames
      if(max(grp) != 1){stop("TSPs can only be calculated for variables with two classes")}
      tsp <- ts.pair(dat,grp)
      tsp$labels <- labels
      return(tsp)
    }
    else{
        labels <- as.character(unique(grp))
        grp <- make.consecutive.int(grp)
        dat <- exprs(dat)
        rownames(dat) <- genenames
        if(max(grp) != 1){stop("TSPs can only be calculated for variables with two classes")}
        tsp <- ts.pair(dat,grp)
        tsp$labels <- labels
        return(tsp)
      }
  }
  else{
    labels <- as.character(unique(grp))
    grp <- make.consecutive.int(grp)
    if(is.null(rownames(dat))){rownames(dat) <- as.character(1:dim(dat)[1])}
    if(max(grp) != 1){stop("TSPs can only be calculated for variables with two classes")}
    tsp <- ts.pair(dat,grp)
    tsp$labels <- labels
  }
  return(tsp)
}

### C Version of the TS-Pair Calculation Function

ts.pair <- function(dat, grp){
  labels <- as.character(unique(grp))
  out <- .Call("ctspair",
               as.double(dat),
               as.double(grp),PACKAGE="tspair")
  index <- matrix(out[[1]],ncol=2,byrow=T)
  index <- matrix(index[,2:1],ncol=2)
  tspscore <- out[[2]]
  tspdat <- dat[c(index[,1],index[,2]),]
  
  # Determine if there are ties, and if so, calculate the tiebreaker
  score <- NULL
  if(out[[3]] > 1){
    score <- rep(0,out[[3]])
    R <- apply(dat,2,rank)
    for(i in 1:out[[3]]){
      r1 <- R[index[i,1],]
      r2 <- R[index[i,2],]
      score[i] <- abs(mean(r1[grp==1] - r2[grp==1]) - mean(r1[grp==0] - r2[grp==0]))
    }
  }

  
  tsp <- list(index = index,tspscore=tspscore,score=score,grp=grp,tspdat=tspdat,labels=labels)
  class(tsp) <- "tsp"
  return(tsp) 
}


### Optimized R Version of the TS-Pair Calculation Function 

# ts.pair <- function(dat,grp){
#  # Set up the hat matrix
#  n <- dim(dat)[2]
#  m <- dim(dat)[1]
#  cdat <- dat
#  labels <- as.character(unique(grp))
#  
#  xx <- cbind(rep(1,n),grp)
#  hatxx <- solve(t(xx) %*% xx) %*% t(xx)
#  dat <- t(dat)
#
#  # Initialize the pairs
#  out <- c(0,0,0)
#
#  # Loop through and do the regressions m at a time
#  for(i in 1:(m-1)){
#    yy <- dat > dat[,1]
#    fit <- abs(hatxx %*% yy)[2,]	  
#    mx <- max(fit)
#    ind <- which(fit == mx) + (i-1)
#    out <- rbind(out,cbind(rep(i,length(ind)),ind,rep(mx,length(ind))))
#    dat <- dat[,-1]
#  }
#
#  # Find the maximum and return the indices
#  tspind <- which(out[,3] == max(out[,3]))
#
#  # Determine if there are ties, and if so, calculate the tiebreaker
#  score <- NULL
#  if(length(tspind) > 1){
#    score <- rep(0,length(tspind))
#    for(i in 1:length(tspind)){
#      r1 <- rank(cdat[out[tspind[i],1],])
#      r2 <- rank(cdat[out[tspind[i],2],])
#      score[i] <- abs(mean(r1[grp==1] - r2[grp==1]) - mean(r1[grp==0] - r2[grp==0]))
#    }
#  }
#  
#  tsp <- list(index = out[tspind,1:2], tspscore = out[tspind,3][1], score = score, grp = grp , tspdat = cdat[c(out[tspind,1],out[tspind,2]),],labels=labels)
#  class(tsp) <- "tsp"
#  return(tsp)
#}




tspplot <- function(tspobj){
  ntsp <- 1
  if(!is.null(tspobj$score)){ntsp <- length(tspobj$score)}
  grp <- tspobj$grp
  labels <- tspobj$labels
  dat <- tspobj$tspdat
  index <- tspobj$index
  
  cat("Number of TSPs: ",ntsp,"\n")

  if(ntsp == 1){
    cat("TSP 1 \n")
    par(mar=c(4,4,4,4))
    plot(dat[1,],dat[2,],xlab=paste("Gene:",rownames(dat)[1],"Expression"),ylab=paste("Gene:",rownames(dat)[2],"Expression"),main=paste("Groups:",labels[1],"= Red |",labels[2],"= Blue; Score:",round(tspobj$tspscore,3)),type="n")
    points(dat[1,grp==0],dat[2,grp==0],col="red",pch=19)
    points(dat[1,grp==1],dat[2,grp==1],col="blue",pch=19)
    abline(c(0,1),lwd=2)
  }
  else{
    for(i in 1:ntsp){
      par(mar=c(4,4,4,4))
      plot(dat[i,],dat[(i + ntsp),],xlab=paste("Gene:",rownames(dat)[i],"Expression"),ylab=paste("Gene:",rownames(dat)[(i+ntsp)],"Expression"),type="n")
      mtext(paste("Groups:",labels[1],"= Red |",labels[2],"= Blue; Score:",round(tspobj$tspscore,3)),line=1)
      points(dat[i,grp==0],dat[(i+ntsp),grp==0],col="red",pch=19)
      points(dat[i,grp==1],dat[(i+ntsp),grp==1],col="blue",pch=19)
      abline(c(0,1),lwd=2)
      readline(paste("TSP",i,": Hit return for next TSP.\n"))
    }
  }
}

plot.tsp <- function(x,y,...){tspplot(x)}

print.tsp <- function(x,...){
  tspobj <- x
  ntsp <- 1
  if(!is.null(tspobj$score)){ntsp <- length(tspobj$score)}
  if(ntsp > 5){
    cat(c("tsp object with:",ntsp, "TSPs\n"))
    cat(c("TSP Score:", round(tspobj$tspscore,2)))
  }
  else{
    cat(c("tsp object with:",ntsp, "TSPs\n"))
    cat(c("Pair:\t\tTSP Score\t Tie-Breaker\t\tIndices\n"))
    if(ntsp==1){
      cat(c("TSP",1,":","\t",round(tspobj$tspscore,2),"\t\t","NA","\t\t\t\t",tspobj$index,"\n"))
    }
    else{
      for(i in 1:length(tspobj$score)){
        cat(c("TSP",i,":","\t",round(tspobj$tspscore,2),"\t\t",round(tspobj$score[i],2),"\t\t\t\t",tspobj$index[i,],"\n"))
      }
    }
  }
}

summary.tsp <- function(object,select=NULL,printall=FALSE,...){
  tspobj <- object
  grp <- tspobj$grp
  ntsp <- 1
  if(!is.null(tspobj$score)){ntsp <- length(tspobj$score)}
  if(ntsp == 1){printall <- TRUE}
  grplabels <- character(length(grp))
  grplabels[grp==0] <- tspobj$labels[1]
  grplabels[grp==1] <- tspobj$labels[2]
  if(!is.null(select)){
    cat(paste("Data for TSP:", select,"\n"))
    z <- tspobj$tspdat[select,] < tspobj$tspdat[(select + ntsp),]
    print(table(z,grplabels,dnn=list(paste("1(Gene",rownames(tspobj$tspdat)[select]," < Gene", rownames(tspobj$tspdat)[(select + ntsp)],")"),"Group Labels")))
  }
  if(is.null(select) & printall==FALSE){
    cat(paste("There are",ntsp,"TSPs\n\n"))
    for(i in 1:ntsp){
      cat(paste("Data for TSP:", i,"\n"))
      z <- tspobj$tspdat[i,] < tspobj$tspdat[(i + ntsp),]
      print(table(z,grplabels,dnn=list(paste("1(Gene",rownames(tspobj$tspdat)[i]," < Gene", rownames(tspobj$tspdat)[(i + ntsp)],")"),"Group Labels")))
      cat("\n\n")
      readline("Hit return for next TSP")
    }
  }
  if(is.null(select) & printall==TRUE){
    cat(paste("There are",ntsp,"TSPs\n\n"))
    for(i in 1:ntsp){
      cat(paste("Data for TSP:", i,"\n"))
      z <- tspobj$tspdat[i,] < tspobj$tspdat[(i + ntsp),]
      print(table(z,grplabels,dnn=list(paste("1(Gene",rownames(tspobj$tspdat)[i]," < Gene", rownames(tspobj$tspdat)[(i + ntsp)],")"),"Group Labels")))
      cat("\n\n")
    }
  }
}

tspsig <- function(dat,grp,B=50,seed=NULL){
  if(class(dat) == "ExpressionSet"){
    if(length(grp) == 1){
      grp <- make.consecutive.int(pData(dat)[,grp])
      if(max(grp) != 1){stop("TSPs can only be calculated for variables with two classes")}
      dat <- exprs(dat)
    }
    else{
      grp <- make.consecutive.int(grp)
      if(max(grp) != 1){stop("TSPs can only be calculated for variables with two classes")}
      dat <- exprs(dat)
    }
  }
  grp <- make.consecutive.int(grp)
  if(!is.null(seed)){set.seed(seed)}
  score <- ts.pair(dat,grp)$tspscore
  score0 <- rep(0,B)
  prog <- progressBar()
  for(i in 1:B){
    score0[i] <- ts.pair(dat,sample(grp))$tspscore
    prog <- progressBar(i/B,prog)
  }
  p <- (sum(score0 > score) + 1)/(length(score) + 1)
  return(list(p = p,nullscores = score0))
}


predict.tsp <- function(object,dat=NULL,select=NULL,...){
  tspobj <- object
  ntsp <- 1
  if(!is.null(tspobj$score)){ntsp <- length(tspobj$score)}
  grp <- tspobj$grp
  grplabels <- character(length(grp))
  grplabels[grp==0] <- tspobj$labels[1]
  grplabels[grp==1] <- tspobj$labels[2]

  if(is.null(select) & ntsp > 1){select <- which.max(tspobj$score)}
  if(is.null(select) & ntsp == 1){select <- 1}
  
  if(select < 1 | select > ntsp){stop(paste("TSP selection must be between 1 and",ntsp))}
  
  if(is.null(dat)){
    predict <- character(length(grp))
    z <- tspobj$tspdat[select,] < tspobj$tspdat[(select + ntsp),]
    predict[z == 0] <- tspobj$labels[which.max(table(z,grplabels)[1,])]
    predict[z == 1] <- tspobj$labels[which.max(table(z,grplabels)[2,])]
    return(predict)
  }else{
    tspnames <- rownames(tspobj$tspdat)[c(select,(select+ntsp))]
    if(class(dat) == "matrix"){
      if(is.null(rownames(dat))){
        cat("No rownames found, using indices \n")
        tspnames <- as.numeric(tspnames)
        if(any(!(tspnames %in% 1:dim(dat)[1]))){stop("Rownames of new data do not include the TSP names")}
      }else{
        if(any(!(tspnames %in% rownames(dat)))){stop("Rownames of new data do not include the TSP names")}
      }
      predict <- character(dim(dat)[2])
      z <- tspobj$tspdat[select,] < tspobj$tspdat[(select + ntsp),]
      w <- dat[tspnames[1],] < dat[tspnames[2],]
      predict[w == 0] <- tspobj$labels[which.max(table(z,grplabels)[1,])]
      predict[w == 1] <- tspobj$labels[which.max(table(z,grplabels)[2,])]
      return(predict)
    }
    if(class(dat) == "ExpressionSet"){
      if(is.null(featureNames(dat))){
        cat("No featureNames info found, using indices \n")
        tspnames <- as.numeric(tspnames)
        if(any(!(tspnames %in% 1:dim(exprs(dat))[1]))){stop("Rownames of new data do not include the TSP names")}
        dat <- exprs(dat)
      }else{
        if(any(!(tspnames %in% featureNames(dat)))){stop("Rownames of new data do not include the TSP names")}
        genenames <- featureNames(dat)
        dat <- exprs(dat)
        rownames(dat) <- genenames
      }
      predict <- character(dim(dat)[2])
      z <- tspobj$tspdat[select,] < tspobj$tspdat[(select + ntsp),]
      w <- dat[tspnames[1],] < dat[tspnames[2],]
      predict[w == 0] <- tspobj$labels[which.max(table(z,grplabels)[1,])]
      predict[w == 1] <- tspobj$labels[which.max(table(z,grplabels)[2,])]
      return(predict)
    }
  } 
}



make.consecutive.int <- function(y) {
  oldWarn = getOption("warn")
  ## Turn off warnings.
  options(warn = -1)

  if(is.null(y)) {return(NULL)}
  
  if(!is.vector(y))
    y = as.vector(as.character(y))

  out <- as.integer(as.factor(as.character(y)))-1
  
  options(warn = oldWarn)
  
  return(out)
}

  
attr(ts.pair, "source") <- NULL
attr(tspplot,"source") <- NULL
attr(tspcalc,"source") <- NULL
attr(tspsig,"source") <- NULL


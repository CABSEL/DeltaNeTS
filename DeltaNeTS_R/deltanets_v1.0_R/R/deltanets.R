#' DeltaNeTS
#' 
#' The main function for DeltaNeTS for inferring gene perturbations for each sample. 
#' @import glmnet foreach doParallel
#' 
#' @param lfc The numeric matrix or data.frame of log2FC data. Each row represents a gene and each column represents a sample.
#' @param slope The slope matrix from log2FC data. This matrix can be obtained using the function \code{generateSlope()}. If the data are not time-series, set slope to an empty matrix (i.e. \code{slope}=\code{NULL}).
#' @param glist The list of genes following the same order as the one in \code{lfc}. If glist is \code{NULL}, the default is just a numbering from 1 to n. 
#' @param grn The list of TF-gene interactions (optional). The gene symbol in \code{grn} should be the same style of the one in \code{glist}. If glist is not available, use indices in \code{grn} instead of gene symbol. If grn is used, ridge regression is used. Otherwise, lasso regression is used to infer A and P.
#' @param kfold The number of folds used in the k-fold cross validation.
#' @param par A Boolean variable \code{TRUE} or \code{FALSE} indicating whether to use parallel computing. The default is \code{FALSE} (no parallel computation).
#' @param numCores The number of CPU cores to be used for parallel computing. This parameter is considered only if par is \code{TRUE}. The default is 4.
#' 
#' @return a list of two matrices
#' \item{P}{The matrix of gene perturbations. Each row corresponds to a gene following the same order as the one in the log2FC data, while each column corresponds to samples as defined in the tobject.}
#' \item{A}{The n*n matrix of the inferred gene regulatory network. The (i,j)-th element of the matrix corresponds to the regulatory impact of gene j on gene i. The rwos and columns of the matrix correspond to genes following the same order as the one in the log2FC data.}
#' 
#' @export
deltanets <- function(lfc,slope=NULL, glist=NULL, grn=NULL, kfold=10,par=FALSE,numCores=4){
  
  norm_vec <- function(x) sqrt(sum(x^2))
  lfc <- lfc/apply(lfc, 1, norm_vec)*sqrt(dim(lfc)[2]-1)
  n <- dim(lfc)[1]
  m <- dim(lfc)[2]

  if(!is.null(slope)){
    slope <- slope/apply(slope, 1, norm_vec)*sqrt(dim(slope)[2]-1)
    ms <- dim(slope)[2]
  }
  if(is.null(glist)){
    glist <- 1:n
  }
  if(!is.null(grn)){
    grn = as.matrix(grn)
    v = length(intersect(unique(matrix(grn,ncol=1)),glist))
    if (v==0){stop("The gene symbol style in grn should be the same as in glist.")}
    grn_idx = cbind(match(grn[,1],glist),match(grn[,2],glist))
    grn_idx = grn_idx[which(apply(is.na(grn_idx),1,sum)==0),] ## remove edges having no match with gene symbols in glist
    grn_idx = grn_idx[which(grn_idx[,1]-grn_idx[,2]!=0),] ## remove self-loop
    grn_usage = TRUE
    dgi = unique(grn_idx[,2])
    alpha = 0 ## ridge regression when grn is used
  }else{
    grn_usage = FALSE
    dgi = 1:n
    alpha = 1 ## lasso regression when grn is not used
  }
  
  
  if(!is.null(slope)){
    X0 <- cbind(lfc,slope)
    Im <- rbind(diag(m),matrix(0,nrow = ms,ncol = m))
  }else{
    X0 <- lfc
    Im <- diag(m)
  }
  
  
  #DeltaNeTS with Lasso 
  
  if(par){
    library(doParallel)
    library(foreach)
    cl <- makeCluster(numCores,outfile='')
    registerDoParallel(cl)
    
    
    Beta <- foreach(j=1:length(dgi),.combine = cbind,.packages = "glmnet")%dopar%{
      cat(sprintf("DeltaNeTS is running with parallel computing...(%4d/%d)\n",j,length(dgi)))
      
      if(grn_usage){tfi = grn_idx[which(grn_idx[,2]==dgi[j]),1]}else{
        tfi <- 1:n
        tfi <- tfi[-dgi[j]] ## remove self-loop when grn is not applied
        }
      
      X <- t(lfc[tfi,])
      if(!is.null(slope)){ 
        if(length(tfi)>1) X <- rbind(X,t(slope[tfi,]))
        else{
          X <- c(X,t(slope[tfi,]))
        }
      }
      
      X <- X0[tfi,]
      y <- X0[dgi[j],]
      if(length(tfi)>1){
        Xin <- cbind(t(X),Im)
      }else Xin <- cbind(X,Im)
      
      
      cvfit <- cv.glmnet(Xin, y, alpha = alpha, intercept = FALSE,nfolds = kfold)
      coef.cvfit <- coef.cv.glmnet(cvfit,s='lambda.min')
      bj <- vector(mode =  "integer",length = dim(Xin)[2])
      bj[coef.cvfit@i] <- coef.cvfit@x
      r <- vector(length = (n+m))
      r[c(tfi,(n+1):(n+m))] <- bj
      r
    }
  }else{

    Beta <- matrix(0,nrow = (n+m), ncol = length(dgi))
    for (j in 1:length(dgi)) {
      cat(sprintf("DeltaNeTS is running...(%4d/%d)\n",j,length(dgi)))
      if(grn_usage){
        tfi = grn_idx[which(grn_idx[,2]==dgi[j]),1]
        }else{
        tfi <- 1:n
        tfi <- tfi[-dgi[j]]
      }
      
      X <- X0[tfi,]
      y <- X0[dgi[j],]
      if(length(tfi)>1){
      Xin <- cbind(t(X),Im)
      }else Xin <- cbind(X,Im)
      
      cvfit <- cv.glmnet(Xin, y, family='gaussian',alpha = alpha, intercept = FALSE, nfolds = kfold)
      coef.cvfit <- coef.cv.glmnet(cvfit,s='lambda.min')
      bj <- vector(mode =  "integer",length = dim(Xin)[2])
      bj[coef.cvfit@i] <- coef.cvfit@x
      Beta[c(tfi,(n+1):(n+m)),j] <- bj
    }
  }
  Beta_sc <- matrix(0,nrow = (n+m), ncol = n)
  Beta_sc[,dgi] <- Beta
  
  A <- t(Beta_sc[1:n,])
  P <- t(Beta_sc[(n+1):(n+m),])
  return(list(P=P,A=A))
}
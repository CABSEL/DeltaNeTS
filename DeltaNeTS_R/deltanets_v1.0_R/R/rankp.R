#' rank genes from P matrix
#' 
#' The function for ranking genes per sample, based on absolute magnitudes in estimates of perturbation impacts
#' @param DeltaNeTSobj The object returned by deltanet function.
#' @param glist The list of genes corresponding to the rows of P matrix in DeltaNeTSobj. If glist is not given, gene symbols are just numbers (e.g. G1, G2, G3, ...)
#' @return a list of rankOfGenes and rankedGenes
#' \item{rankOfGenes}{The matrix of gene ranks. Rows and columns correspond to genes and samples in the same order as the one in P matrix, respectively.A maximum rank is assigned for the genes with zero values in the P matrix}
#' \item{rankedGenes}{A list of ranked gene lists. Each element in the list includes a ranked gene list for the corresponding sample. The genes closer to the front are the gene targets with the higher confidence. }
#' 
#' @export
rankp <- function(DeltaNeTSobj, glist=NULL){

  P = result$P
  n = dim(P)[1]
  
  if(is.null(glist)){glist = paste("G",1:n,sep="") }


oi <- apply(abs(P),2,function(x) order(x,decreasing=TRUE))
rankOfGenes <- apply(oi,2,order)
rankOfGenes[which(P==0)]=n

oi[which(P==0)]=0
oilist = apply(oi,2,function(x) x[x>0])

rankedGenes = lapply(oilist,function(x) glist[x])

return(list(rankOfGenes=rankOfGenes, rankedGenes=rankedGenes))
}


#' 
#' 
#' @param data Output from \code{read.GWASpoly}
#' @param models Vector of model names
#' @param traits Vector trait names (by default, all traits)
#' @param max.iter maximum number of iteration
#' @param force.PCs if true it adds PC1 and PC2 to the models, if false, the PCs are only added if model BIC reduce 
#' @param verbose TRUE/FALSE whether to suppress output charting progress
#' 
#' @import bigstep
#' 
#' @export
GWASpoly_glm <- function(data,
                         models,
                         traits=NULL,
                         max.iter = 100,
                         force.PCs = FALSE,
                         verbose=F) {
  
  # Checks
  if(!is.null(traits)){
    if(!all(colnames(data@pheno)[-1] %in% traits))
      stop(paste("Variable name:", traits[which(!colnames(data@pheno) %in% traits)],"not found in phenotypic data"))
  } else traits <- colnames(data@pheno)[-1]
  colnames(data@pheno)[1] <- "id"
  
  # Read geno file, dosage codification
  marker.data <- data@geno
  
  # Guarantee row order match with BLUES - individuals in row and markers in column
  X <- as.matrix(marker.data[match(data@pheno$id,rownames(marker.data)),])
  m <- ncol(X) # number of markers
  
  scores <- list()
  for(i in 2:ncol(data@pheno)){
    for(j in 1:length(models)){ # Expand for other models
      pheno <- data@pheno[,i]  
      prep1 <- prepare_data(y=pheno,X)
      data1 <- reduce_matrix(prep1,minpv=0.1) # remove uncorrelated variables
      
      # Add variables while the BIC is reducing
      data2 <- fast_forward(data1,crit = bic)
      
      iter=1
      while (iter < max.iter & !setequal(data1$model,data2$model)) {  # loop goes for max.iter iterations or selected variables are the same
        suppressMessages(data1 <- multi_backward(data2,crit=bic)) # remove variables as long as they reduce the BIC
        suppressMessages(data2 <- fast_forward(data1,crit=bic))   # add variables as long as they reduce the BIC
        iter <- iter+1
      }
      print(paste("Total Iterations:",iter))
      
      # markers in the baseline model (bm)
      bm.qtn <- match(unique(data2$model),colnames(X))
      
      # Principal Coordinates (from multi-dimensional scaling)
      dm <- dist(X)
      mds <- cmdscale(dm,k=2)
      
      data3 <- data.frame(y=data2$y,pc1=mds[,1],pc2=mds[,2],X)

      #data3[,which(colnames(data3) == "Chr13_022193410")]
      
      # This model add PC1, PC2 and all variables
      mf <- paste("y",paste(c("1+pc1+pc2",data2$model),collapse="+"),sep="~")
      bm12 <- lm(as.formula(mf),data3)
      
      # Compare BIC with and without PCs
      if(!force.PCs){
        mf <- paste("y",paste(c("1+pc1",data2$model),collapse="+"),sep="~")
        bm1 <- lm(as.formula(mf),data3)
        mf <- paste("y",paste(c("1+pc2",data2$model),collapse="+"),sep="~")
        bm2 <- lm(as.formula(mf),data3)
        mf <- paste("y",paste(c("1",data2$model),collapse="+"),sep="~")
        bm <- lm(as.formula(mf),data3)
        all.mod <- list("PCs not included"=bm, "PC1 included" = bm1, "PC2 included"=bm2, "PC1 and PC2 included"=bm12)
        
        best <- which.min(sapply(all.mod, BIC))
        print(paste(names(best), "in the model based on BIC changes."))
        
        bm <- all.mod[[best]]
      } else bm <- bm12
      
      
      scope <- colnames(X)[bm.qtn] 
      if(any(is.na(bm$coefficients))) scope <- scope[-which(scope %in% names(bm$coefficients)[which(is.na(bm$coefficients))])]
      nonscope <- colnames(X)[-which(colnames(X) %in% scope)]
      #p-values for markers in the bm computed by backwards elimination
      drop1.ans <- drop1(bm,scope=scope,test="F")
      
      #p-values for other markers computed by forward entry
      add1.ans <- add1(bm,scope=nonscope,test="F")
      
      pval1 <- numeric(m)
      names(pval1) <- colnames(X)
      pval1[match(rownames(drop1.ans)[-1],names(pval1))] <- drop1.ans$`Pr(>F)`[-1]
      pval1[match(rownames(add1.ans)[-1],names(pval1))] <- add1.ans$`Pr(>F)`[-1]
      
      score.blink <- as.data.frame(-log10(pval1))
      colnames(score.blink) <- models[j]
    }
    scores[[i-1]] <- score.blink
  }
  
  names(scores) <- colnames(data@pheno)[-1]
  
  return(new("GWASpoly.fitted",
             map=data@map,
             pheno=data@pheno,
             geno=data@geno,
             ploidy=data@ploidy,
             scores=scores))  
}

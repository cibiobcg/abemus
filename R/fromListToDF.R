#' fromListToDF
#'
#' @author Davide Prandi
#' @param inputlist list to merge
#' @return data.frame
fromListToDF <- function(inputList){
  if (is.null(inputList)){return(NULL)}
  #check if some is null and remove
  nullPositions <- which(sapply(inputList,is.null))
  if (length(nullPositions) > 0){
    inputList <- inputList[-nullPositions]
  }
  firstEl <- inputList[[1]]
  inputList <- lapply(inputList, function(x){ matrix(unlist(x), ncol=ncol(x))} )
  outDF <- as.data.frame(do.call(rbind, inputList),stringsAsFactors=F)
  colnames(outDF) <-names(firstEl)

  for(idx in c(1:ncol(outDF))){
    if (class(firstEl[[idx]]) == "logical"){
      if (is.na(firstEl[[idx]])){
        class(outDF[[idx]]) <- "numeric"
      }else if (outDF[1,idx] == "TRUE" ||outDF[1,idx] == "FALSE"  ){
        outDF[,idx] <- as.logical(outDF[,idx]) * 1
      }
      class(outDF[[idx]]) <- "numeric"
    }else{
      class(outDF[[idx]]) <- class(firstEl[[idx]])
    }
  }
  return(outDF)
}

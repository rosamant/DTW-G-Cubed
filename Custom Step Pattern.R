stepPattern <- function(v,norm=NA) {
  obj <- NULL;
  if(is.vector(v)) {
    obj <- matrix(v,ncol=4,byrow=TRUE);
  } else if(is.matrix(v)) {
    obj <- v;
  } else {
    stop("stepPattern constructor only supports vector or matrix");
  }
  class(obj)<-"stepPattern";
  attr(obj,"npat") <- max(obj[,1]);
  attr(obj,"norm") <- norm;
  return(obj);
}

asymmetricP1.1 <- stepPattern(c(
  1, 1 , 2 , -1 ,
  1, 0 , 1 , 1 ,
  1, 0 , 0 , 1 ,
  2, 1 , 1 , -1 ,
  2, 0 , 0 ,  1 ,
  3, 2 , 1 , -1 ,
  3, 1 , 0 ,  2 ,
  3, 0 , 0 ,  2
),"N")

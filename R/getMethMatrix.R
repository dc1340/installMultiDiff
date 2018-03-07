#' Get the raw methylation matrix from a methBase object
#'
#' Pulls out raw methylation values for all samples. Missing data is filled with NAs.
#' @param mbase A methylBase object
#' @param long.rownames Inserts rownames which give chr and pos. Default=F
#' @keywords methylBase,methylKit
#' @export

getMethMatrix<-function(mbase, long.rownames=F){
  meth_matrix=matrix(nrow=nrow(mbase), ncol=length(mbase@sample.ids))
  mdata=getData(mbase)
  meth_matrix=mdata[ , mbase@numCs.index]/mdata[  , mbase@coverage.index]
  colnames(meth_matrix)=mbase@sample.ids
  if (long.rownames){
    rownames(meth_matrix)=with(getData(mbase), paste(chr, start, strand, sep='.'))
  }
  return (meth_matrix)
}

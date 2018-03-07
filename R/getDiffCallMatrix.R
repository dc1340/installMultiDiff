getDiffCallMatrix <- function(multiDiff, meth.thresh=10, q.thresh=0.01, long.rownames=F){ 
  
  pos_dat=multiDiff[[1]]
  cur_terms=names(multiDiff[[3]])
  
  all_sites=with(pos_dat, paste(chr, start, strand, sep='.'))
  
  #diffMatrix=getDiffMatrix(multiDiff)  
  
  call_mat=matrix(nrow=nrow(pos_dat), ncol=length(cur_terms))
  colnames(call_mat)=c( cur_terms)
  rownames(call_mat)=rownames(pos_dat)
  
  for (i in 1:length(cur_terms)){
    cur_dat=multiDiff[[2]][ , , i]
    call_mat[ , i]=abs(cur_dat[ , 1])>=meth.thresh & (cur_dat[ , 3]<=q.thresh)
  }
  if (long.rownames){
    rownames(call_mat)=with(multiDiff[[1]], paste(chr, start, strand, sep='.'))
  }
  
  return(call_mat)
}

getSummaryDiffCall <-function( multiDiff, meth.thresh=10, q.thresh=0.01){
  coefs=names(multiDiff[[3]])
  sum_mat=matrix(data=0, nrow=3, ncol=length(coefs))
  colnames(sum_mat)=coefs
  rownames(sum_mat)=c('-1', '0', '1')
  
  if (nrow(multiDiff[[1]])==0){
    return(sum_mat)
  }
  
  call_mat=getDiffCallMatrix(multiDiff, meth=meth.thresh, q=q.thresh)
  diff_mat=getDiffMatrix(multiDiff)
  tmp_mat=call_mat * sign(diff_mat)

  


  sum_mat[ 1, ]=colSums(tmp_mat==-1, na.rm = T)
  sum_mat[ 3, ]=colSums(tmp_mat==1,na.rm = T)
  sum_mat[ 2, ]=rep(nrow(call_mat), length(coefs))-(sum_mat[ 1, ]+sum_mat[ 3, ])
  
  return(sum_mat)
}
 
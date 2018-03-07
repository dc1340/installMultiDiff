filterMultiOnMulti <-function(baseMulti, filtMulti, inv=F, coef, meth.thresh=25, q.thresh=0.01){

  filt_sites=getDiffCallMatrix(filtMulti,
                               meth.thresh = meth.thresh,
                               q.thresh=q.thresh)[ , coef]
  
  filt_sites=with(filtMulti[[1]], paste(chr,start,end,strand, sep=".")
                  )[ filt_sites ]
  
  base_sites=with(baseMulti[[1]], paste(chr,start,end,strand , sep="."))
  
  keep_sites=base_sites %in% filt_sites
  
  if (inv) keep_sites=!keep_sites
  
  outMulti=baseMulti
  outMulti[[1]]=outMulti[[1]][ keep_sites, ]
  outMulti[[2]]=outMulti[[2]][ keep_sites, , ]
  
  return(outMulti)
}
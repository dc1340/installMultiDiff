convertMethBaseToBSeqData<-function(methBase){
  
  sample_count=length(methBase@sample.ids)
  list_dat=list(1:sample_count)
  
  cur_dat=getData(methBase)
  
  for(i in 1:sample_count){
    
                      
    list_dat[[i]]=data.frame(chr=cur_dat$chr,
                             pos=cur_dat$start,
                             N=cur_dat[, methBase@coverage.index[[i]]],
                             X=cur_dat[, methBase@numCs.index[[i]]]
                             )
  }
  
  BSobj <- makeBSseqData( list_dat,
                           methBase@sample.ids )
  
}

convertMultiDSStoMultiDiff<- function(multiDSS){
  output=list(pos.dat, glm.dat, coef.index, slot_data)  
}

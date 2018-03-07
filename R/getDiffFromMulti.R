getDiffFromMulti<-function(multiDiff, coef=NA){
  
  cur_index=as.integer(multiDiff[[3]][ coef])
  
  if(nrow(multiDiff[[1]])==1){
    cur_dat=multiDiff[[2]][ , cur_index]
    x=data.frame(multiDiff[[1]],pvalue=cur_dat[[2]] ,
                 qvalue=cur_dat[[3]],meth.diff=cur_dat[[1]],stringsAsFactors=F )
  } else{
    cur_dat=multiDiff[[2]][ , , cur_index]
    x=data.frame(multiDiff[[1]],pvalue=cur_dat[ , 2] ,
                 qvalue=cur_dat[ , 3],meth.diff=cur_dat[ , 1],stringsAsFactors=F )
  }
  

  rownames(x)=rownames(multiDiff[[1]])
  # make a dataframe and return it
  
  ## Make the new methylDiff object
  
  slot_data=multiDiff[[4]]
  output=new("methylDiff",x,
             sample.ids=slot_data$sample.ids,
             assembly=slot_data$assembly,
             context=slot_data$context,
             treatment=slot_data$treatment,
             destranded=slot_data$destranded,
             resolution=slot_data$resolution)
  

  
  return(output)
}
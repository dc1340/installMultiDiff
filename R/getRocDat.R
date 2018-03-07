getRocDat<-function( multiDiff, site_states,base_q.thresh=1,
                     base_meth.thresh=0, add.q.roc=F){
  cur_terms=names(multiDiff[[3]])
  num_sites=nrow(multiDiff[[1]])
  
  roc_dat=data.frame(matrix(nrow = 203*( length(cur_terms)+1), ncol=5) )
  colnames(roc_dat)=c('Term', 'type', 'thresh', 'FPR' , 'TPR' )
  diff_mat=abs(getDiffMatrix(multiDiff))
  q_mat=multiDiff[[2]][ , 3, ]
  cond_pos=colSums(site_states)
  
  base_q_mat=multiDiff[[2]][ , 3, ]<=base_q.thresh  
  
  i=1 
  
  for (thresh in 0:100){
      cur_diff_mat=diff_mat>=thresh
      
      cur_call=cur_diff_mat & base_q_mat
      total_tpr=0
      total_fpr=0
      
      for (j in 1:length(cur_terms)){
        term=cur_terms[[j]]
        cur_tpr=sum(site_states[ ,  j ]& cur_call[ , j ])
        cur_fpr=sum((!site_states[ , j ]) & cur_call[ , j ])
        
        total_tpr=total_tpr+cur_tpr
        total_fpr=total_fpr+cur_fpr
        
        cur_tpr=cur_tpr/cond_pos[[j]]
        cur_fpr=cur_fpr/(num_sites-cond_pos[[j]])
        
        roc_dat[ i,  ]=list(term , 'meth', thresh, cur_fpr, cur_tpr)
        
        i=i+1
      }
      
      total_fpr=total_fpr/sum(site_states==0)
      total_tpr=total_tpr/sum(site_states==1)
      
      roc_dat[ i,  ]=list('All Terms' , 'meth', thresh, total_fpr, total_tpr)
      i=i+1
      
  }

  if (add.q.roc){
    for (thresh in seq(-0.1,1, by=0.01)){
      cur_q_mat=q_mat<=thresh
      cur_call=cur_q_mat #& base_meth_mat
      for (j in 1:length(cur_terms)){
        term=cur_terms[[j]]
        cur_tpr=sum(site_states[ ,  j ]& cur_call[ , j ])/cond_pos[[j]]
        cur_fpr=sum((!site_states[ , j ]) & cur_call[ , j ])/(num_sites-cond_pos[[j]])
        roc_dat[ i,  ]=list(term , 'q.value', thresh, cur_fpr, cur_tpr)
        i=i+1
      }
    }
  }
  #base_meth_mat=diff_mat>=base_meth.thresh


  
  
  return(roc_dat)
}

getSpectrumCall<-function( multiDiff, q.thresh=0.01){
  cur_terms=names(multiDiff[[3]])
  num_sites=length(multiDiff[[1]])
  
  diff_mat=abs(getDiffMatrix(multiDiff))
  q_mat=multiDiff[[2]][ , 3, ]<=q.thresh  
  
  threshes=0:100
  spec_dat=data.frame(matrix(nrow = length(threshes),
                             ncol=length(cur_terms)+1) )
  colnames(spec_dat)=c(cur_terms, 'thresh')

  
  for (i in 1:length(threshes)){
    cur_diff_mat=diff_mat>=threshes[[i]]
    
    cur_call=cur_diff_mat & q_mat
    tot_call=colSums(cur_call)
    
    spec_dat[ i,  ]=c(tot_call, threshes[[i]])
    
  }
  
  return(spec_dat)
}


toBinState <-function( state_matrix){
  rowSums(matrix(as.numeric(state_matrix),
               nrow=nrow(state_matrix)) %*% diag(2^c(0:(ncol(state_matrix)-1))))+1
  
}
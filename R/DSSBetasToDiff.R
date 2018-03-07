arcsinLink<-function(x){
  return(asin(2*x-1))
}

sineLink<-function(x){
  if (length(x)>1){
    x[ x<(-pi/2)]=-pi/2
    x[ x>(pi/2)]=pi/2
  }else if (x<(-pi/2)){
    return(0)
  }else if (x>(pi/2)){
    return(1)
  }
  
  return(0.5*(sin(x)+1))
}

#### Get getMaxDiffFromMultiDSS ####
convertDSStoMultiDiff<-function(multiDSS,  cur_design,
                                assembly='mm10',
                                context='CpG',
                                destranded=TRUE,
                                resolution='base'){ 
  
  
  pos.dat=as.data.frame(multiDSS[['gr']])
  colnames(pos.dat)[[1]]='chr'
  pos.dat=pos.dat[ , c('chr', 'start' , 'end', 'strand')]
  
  
  cur_coefs=gsub(' ' , '', strsplit(as.character(multiDSS[[3]])[[2]],
                                    split="+", fix=T)[[1]])
  
  glm.dat=array(rep(0,nrow(pos.dat)*4* length(cur_coefs)),
                dim=c(nrow(pos.dat), 4, length(cur_coefs)))
  
  glm.dat[ , 1, ]=getMaxDiffFromMultiDSS(multiDSS)
  for (i in 1:length(cur_coefs)){
    glm.dat[ , 2:4, i ]=as.matrix(DMLtest.multiFactor(multiDSS,
                                                      coef=i+1)[ ,
                                                                 c('pvals', 'fdrs', 'stat')])
  }
  
  coef.index=1:length(cur_coefs)
  names(coef.index)=cur_coefs
  slot_data=list( sample.ids=rownames(multiDSS$design),
                  assembly=assembly,
                  context=context,
                  treatment=as.integer(rowSums(multiDSS$design)),
                  destranded=destranded,
                  resolution=resolution)
  
  output=list(pos.dat, glm.dat, coef.index, slot_data)
}

addMaxDiffToMultiDSS <-function(multiDSS){
  if (!('max.diff' %in% names(multiDSS))){
    multiDSS[[length(multiDSS)+1]] <-getMaxDiffFromMultiDSS(multiDSS)
    names(multiDSS)[[5]]='max.diff'
  }
  
  
  return(multiDSS)
}


getMaxDiffFromMultiDSS <-function(multiDSS){
  
  cur_design=multiDSS$design
  cur_beta=multiDSS$fit[[ 'beta']]
  cur_diff=convertMultiDSSBetasToDiff(cur_beta, cur_design)
  return(cur_diff)
  #diff_mat=abs(getMaxDiffFromMultiDSS(multiDSS))
}

convertMultiDSSBetasToDiff<-function(beta_matrix, cur_design){
  ##Assume intercept is included
  diff_out=beta_matrix[ , 2:ncol(beta_matrix)]
  coef_numb=ncol(beta_matrix)
  possible_states=laply(0:(2^ncol(cur_design)-1), function (x) { as.integer((intToBits(x))  )  } )[ , 1:ncol(cur_design) ]
  
  for (i in 1:nrow(beta_matrix)){
    cur_coefs=beta_matrix[ i, ]
    if (any(is.na(cur_coefs))){
      diff_out[ i, ]=rep(NA, ncol(diff_out))
    } else {
      for ( j in 2:coef_numb){
        cur_design_col=as.integer(cur_design[ , j-1 ] )
        cur_beta <- cur_coefs[[j]]
        
        cur_rest= cur_coefs [ c(-j) ]
        cur_rest = cur_rest [ c(-1)]
        
        if (length(cur_rest)==1){
          cur_rest=as.matrix(cur_rest)
        }
        
        #cur_se <- summary(cur_obj)$coefficients[ j, 2 ]
        
        #cur_p.value=cur_pvals[[j]]
        
        
        cur_push_meth=sineLink((possible_states[ , c(-(j-1))] %*% cur_rest) + cur_coefs[[1]]+cur_beta)
        cur_base_meth=sineLink((possible_states[ , c(-(j-1)) ] %*% cur_rest) + cur_coefs[[1]])
        
        cur_meth_diff=(cur_push_meth-cur_base_meth)
        
        cur_max_diff=max(cur_meth_diff)
        cur_min_diff=min(cur_meth_diff)
        
        if (abs(cur_max_diff)>=abs(cur_min_diff)){
          cur_meth_diff= cur_max_diff
        } else {
          cur_meth_diff= cur_min_diff
        }
        
        cur_meth_diff=100*cur_meth_diff
        diff_out[ i, j-1]=cur_meth_diff
      }
    }
    
  }
  
  return(diff_out)
}

getMultiDSSqvalues<-function(multiDSS){
  coef_numb=length(strsplit(as.character(multiDSS[[3]])[[2]],
                            split="+", fix=T)[[1]])
  q.out=matrix(nrow=length(multiDSS["gr"]$gr),
               ncol=coef_numb)
  
  for (i in 1:coef_numb){
    q.out[, i]=DMLtest.multiFactor(multiDSS, coef=i+1)[ ,'fdrs']
  }
  return(q.out)
}

getSummaryDiffCallMultiDSS<-function(multiDSS, meth.thresh=15, q.thresh=0.01){
  diff_call=getMultiDSSDiffCall(multiDSS, meth.thresh, q.thresh)
  tmp_mat=diff_call *sign(getMaxDiffFromMultiDSS(multiDSS ))
  
  cur_coefs=gsub(' ' , '', strsplit(as.character(multiDSS[[3]])[[2]],
                                    split="+", fix=T)[[1]])
  
  sum_mat=matrix(nrow=3, ncol=length(cur_coefs))
  colnames(sum_mat)=cur_coefs
  
  rownames(sum_mat)=c('-1', '0', '1')
  
  sum_mat[ 1, ]=colSums(tmp_mat==-1)
  sum_mat[ 3, ]=colSums(tmp_mat==1)
  sum_mat[ 2, ]=rep(nrow(diff_call), length(cur_coefs))-(sum_mat[ 1, ]+sum_mat[ 3, ])
  
  
  return(sum_mat)
}



getMultiDSSDiffCall<-function(multiDSS, meth.thresh=25, q.thresh=0.01){
  cur_terms=gsub(' ' , '', strsplit(as.character(multiDSS[[3]])[[2]],
                                    split="+", fix=T)[[1]])
  
  if ('max.diff' %in% names(multiDSS)){
    diff.matrix=multiDSS['max.diff']
  }else{
    diff.matrix=getMaxDiffFromMultiDSS(multiDSS )
  }
  q.matrix=getMultiDSSqvalues(multiDSS)
  call.matrix=(abs(diff.matrix)>=meth.thresh & q.matrix<q.thresh)
  
  colnames(call.matrix)=c( cur_terms)
  
  return(call.matrix)
}

getDSSRocDat<-function( multiDSS, site_states,base_q.thresh=1,
                        base_meth.thresh=0, add.q.roc=F){
  cur_terms=gsub(' ' , '', strsplit(as.character(multiDSS[[3]])[[2]],
                                    split="+", fix=T)[[1]])
  num_sites=length(multiDSS[[1]])
  
  roc_dat=data.frame(matrix(nrow = 203*( length(cur_terms)+1), ncol=5) )
  colnames(roc_dat)=c('Term', 'type', 'thresh', 'FPR' , 'TPR' )
  diff_mat=abs(getMaxDiffFromMultiDSS(multiDSS))
  q_mat=getMultiDSSqvalues(multiDSS)
  cond_pos=colSums(site_states)
  
  base_q_mat=q_mat<=base_q.thresh  
  
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

getDSSSpectrumCall<-function( multiDSS, q.thresh=0.01){
  cur_terms=gsub(' ' , '', strsplit(as.character(multiDSS[[3]])[[2]],
                                    split="+", fix=T)[[1]])
  num_sites=length(multiDSS[[1]])
  
  threshes=0:100
  spec_dat=data.frame(matrix(nrow = length(threshes),
                             ncol=length(cur_terms)+1) )
  colnames(spec_dat)=c(cur_terms, 'thresh')
  diff_mat=abs(getMaxDiffFromMultiDSS(multiDSS))
  q_mat=getMultiDSSqvalues(multiDSS)
  q_mat=q_mat<=q.thresh
  
  for (i in 1:length(threshes)){
    cur_diff_mat=diff_mat>=threshes[[i]]
    
    cur_call=cur_diff_mat & q_mat
    tot_call=colSums(cur_call, na.rm=T)
    
    spec_dat[ i,  ]=c(tot_call, threshes[[i]])
    
  }
  
  return(spec_dat)
}
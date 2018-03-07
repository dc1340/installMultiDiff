getDiffMatrix<-function(multiDiff){
  pos_dat=multiDiff[[1]]
  #all_sites=with(pos_dat, paste(chr, start, strand,sep='.'))
  
  
  cur_terms=names(multiDiff[[3]])
  
  #meth_diff=matrix(nrow=nrow(pos_dat), ncol=length(cur_terms))
  meth_diff=multiDiff[[2]][ , 1, ]
  colnames(meth_diff)=cur_terms
  rownames(meth_diff)=rownames(pos_dat)
  
  #for (coef in cur_terms){
    
   #meth_diff[ , coef]=multiDiff[[2]][ , 1, multiDiff[[3]][[coef]]]
  #}
  
  return(meth_diff)
}


makeHeatMap <- function(multiDiff, meth.thresh=10,  q.thresh=0.01, make.plot=T, cluster_cols=T,
                        call_col='black', plot_all_sites=F,
                        long_chr=F, relabel=NA, transpose=F, fontsize=18){ 
  
  if (!is.na(relabel)){
    names(multiDiff[[3]])=relabel
  }
  pos_dat=multiDiff[[1]]
  
  chr_string='chr'
  if (long_chr){
   chr_string='chromosome'
  }
  
  plot_all_sites=rep(plot_all_sites, nrow(pos_dat))
  
  cur_terms=names(multiDiff[[3]])
  
  #all_sites=with(pos_dat, paste(chr, start, strand, sep='.'))
  
  diffMatrix=getDiffMatrix(multiDiff)  
  
  callMatrix=getDiffCallMatrix(multiDiff, meth.thresh=meth.thresh, q.thresh = q.thresh)
  
  #cat( length(apply(callMatrix, 1, any)), ' ' , length(plot_all_sites), '\n')
  all_diff_pos=apply(callMatrix, 1, any) | plot_all_sites
  ann_df=data.frame(chr=pos_dat$chr, matrix(as.numeric(callMatrix), ncol=length(cur_terms)))
  colnames(ann_df)=c(chr_string, cur_terms)
  rownames(ann_df)=rownames(diffMatrix)
  
  
  
  cur_term=cur_terms[[1]]

  
  
  #ann_df[ , cur_term]=as.numeric(all_sites %in% diff_sites)
  ann_cols=list()
  ann_cols[[cur_term ]]=c('white', call_col)
  if (length(cur_terms)>1){
      for (cur_term in cur_terms){
      #ann_df[ , cur_term]=as.numeric(all_sites %in% cur_sites)
      ann_cols[[cur_term ]]=c('white', call_col)
      }
  }
  
  
  #diff_sites=unique(diff_sites)
  #cat(length(diff_sites), '\n')
  
  #all_diff_pos=all_sites %in% diff_sites
  ann_df=(ann_df[ all_diff_pos, ] )
  
  if (ncol(diffMatrix)==1){
    diffMatrix=matrix(diffMatrix[ all_diff_pos , ], ncol=1, nrow=sum(all_diff_pos) )
    colnames(diffMatrix)=cur_terms[[1]]
  } else {
    diffMatrix=diffMatrix[ all_diff_pos , ]
  } 
  
  #rownames(ann_df)=all_sites[ all_diff_pos ]
  #rownames(diffMatrix)=all_sites [ all_diff_pos]
  
  
  #  ann_cols[[ 'chr']]=rainbow(length(unique(ann_df$chr)))
  
  #Order column
  ann_df=ann_df[ , (ncol(ann_df):1)]
  rownames(ann_df)=rownames(diffMatrix)
  
  #Order row
  diffMatrix=diffMatrix[ order(ann_df$chr), ]
  ann_df=ann_df[ order(ann_df$chr), ]
  
  #diffMatrix[is.na(diffMatrix)]=0
  if (transpose){
    cur_plot=pheatmap((diffMatrix), cluster_rows=cluster_cols, cluster_cols=F,
                      show_rownames = F, 
                      annotation_legend=F, annotation_row= ann_df,
                      annotation_colors = ann_cols, fontsize=fontsize)+
      theme_set(theme_gray(base_size=fontsize))
  } else {
    cur_plot=pheatmap(t(diffMatrix), cluster_rows=F, cluster_cols=cluster_cols,
                      show_colnames = F, 
                    annotation_legend=F, annotation_col= ann_df,
                    annotation_colors = ann_cols, fontsize=fontsize)+
      theme_set(theme_gray(base_size=fontsize))
}
  if (make.plot){
    print(cur_plot)
    #return(NULL)
  } else{
    #return(cur_plot)
  }
  
  #print(pheatmap(t(diffMatrix), cluster_cols=F, show_colnames = F, annotation=ann_df))
  
  #print(pheatmap(t(diffMatrix[ all_sites %in% diff_sites , ]), cluster_rows=F, show_colnames = F))
  
}
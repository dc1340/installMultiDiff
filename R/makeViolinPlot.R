#' Create a violin and/or barplot for a multiDiff object
#'
#' Allows quick generation of violin and barplots for a multiDiff object.
#' @param multiDiff multiDiff array.
#' @param meth.thresh Effect size threshold for differential calling. Default=0.
#' @param q.thresh Statistical threshold for differential calling. Default=0.01
#' @param make.violin Whether to make the violin plot. Default=T.
#' @param make.bar Whether to make the bar plot. Default=T.
#' @param bw_plot Whether to make the plot black and white. Default=F
#' @param relabel Optional string to relabel covariates. Default=NA.
#' @param bar.ymax Maximum of y-axis in bar plot. Default=NA
#' @param  fontsize What fontsize to use. Default=18.
#' @keywords heatmap
#' @export

makeViolinPlot<- function(multiDiff, meth.thresh=0, q.thresh=0.01,
                          make.violin=T, make.bar=T, make.plot=T, bw_plot=F,
                          relabel=NA, bar.ymax=NA, fontsize=16){
  require(ggplot2)
  require(plyr)

  if (length(relabel)==length(multiDiff[[3]] &  !is.na(relabel)   )){
    names(multiDiff[[3]])=relabel
  }

  if (!make.violin & !make.bar){ return(NULL)}

  cur_terms=names(multiDiff[[3]])


  if (make.violin){
     violin_dat=data.frame(meth.diff=as.vector(multiDiff[[2]][ , 1 , ]),
                        coef=rep(cur_terms, each=nrow(multiDiff[[1]])),
                        q.value=as.vector(multiDiff[[2]][ , 3 , ])
                        )
  violin_dat=subset(violin_dat, abs(meth.diff)>=meth.thresh & q.value<=q.thresh)

    violin_dat$coef=factor(violin_dat$coef, levels=cur_terms)

    cur_violin=qplot(factor(coef), meth.diff, data = violin_dat, geom = "violin", fill=coef)+
      ylim(-100, 100)+ylab('Max. Methylation Difference')+theme(text = element_text(size=fontsize))

 }

  if (make.bar){
    if(nrow(violin_dat)==0){
      bar_dat=data.frame(coef=cur_terms,
                         site_count=0)
      bar_dat$coef=factor(bar_dat$coef, levels=cur_terms)

    }else{

      bar_dat=plyr::count(violin_dat$coef)
      colnames(bar_dat)=c('coef', 'site_count')
      bar_dat$coef=factor(bar_dat$coef, levels=cur_terms)

    }

    #Changed to work with more recent versions of ggplot
    #cur_bar=qplot(x=factor(coef), y=site_count, data=bar_dat, geom='bar', fill=coef, stat='identity')+theme_set(theme_gray(base_size=18))

    cur_bar=qplot(data=violin_dat, x=(coef), geom='bar', fill=coef)+
      scale_y_continuous(label=scientific)+
      ylab('# DMC')+theme(text = element_text(size=fontsize))

    if (!is.na(bar.ymax)){
      cur_bar=cur_bar+ylim(0, bar.ymax)
    }

  }

  if (bw_plot){
    cur_violin=cur_violin+theme_bw()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    cur_bar=cur_bar+theme_bw()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())

  }
  if(!make.bar){
    cur_plot=cur_violin
  } else if (!make.violin){
    cur_plot=cur_bar
  } else{
    cur_plot=multiplot(cur_violin, cur_bar)
  }



  if (make.plot){ print(cur_plot)}




  return(cur_plot)
}

makeBarPlot<- function(multiDiff, meth.thresh=0, q.thresh=0.01,
                       bw_plot=F, fontsize=16) {
  require(ggplot2)
  require(plyr)
  require(scales)

  cur_terms=names(multiDiff[[3]])

  violin_dat=data.frame(meth.diff=as.vector(multiDiff[[2]][ , 1 , ]),
                        coef=rep(cur_terms, each=nrow(multiDiff[[1]])),
                        q.value=as.vector(multiDiff[[2]][ , 3 , ])
  )
  violin_dat=subset(violin_dat, abs(meth.diff)>=meth.thresh & q.value<=q.thresh)

  violin_dat$coef=factor(violin_dat$coef, levels=cur_terms)


  bar_dat=plyr::count(violin_dat$coef)


  colnames(bar_dat)=c('coef', 'site_count')


  bar_dat$coef=factor(bar_dat$coef, levels=cur_terms)

  #cur_bar=qplot(x=factor(coef), y=site_count, data=bar_dat, geom='bar', fill=coef)+theme_set(theme_gray(base_size=18))
  cur_bar=qplot(data=violin_dat, x=(coef), geom='bar', fill=coef)+
    scale_y_continuous(label=scientific)+
    theme(text = element_text(size=fontsize))

  if (bw_plot){
    cur_bar=cur_bar+theme_bw()

  }


  return(cur_bar)
}

scientific_10 <- function(x) {
     parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

scientific <- function(x) {
  parse(text=scientific_format()(x))
}

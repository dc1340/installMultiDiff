#' Filter the sites of a multiDiff array against another
#'
#' Effectively does an inner join between all the sites in the baseMulti,
#' and the differential sites found for "coef" in filtMulti.
#' "inv" can be used to make it a left join instead.
#' on baseMulti.
#' @param baseMulti the base multiDiff object to be filtered
#' @param filtMulti the multiDiff object used for filtering.
#' @param coef Which coefficient in filtMulti to use for filtering.
#' @param meth.thresh Effect size threshold for differential calling. Default=25.
#' @param q.thresh Statistical threshold for differential calling. Default=0.01
#' @param inv The default behavior is to only include sites in the baseMulti
#'     that are called as differential in filtMulti. This flag inverts the filter,
#'     only letting through sites that weren't called. Default=F.
#' @keywords differential,DiffMatrix
#' @export

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

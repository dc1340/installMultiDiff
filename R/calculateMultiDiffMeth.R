#' Calculates the multiDiff Object
#' This function calculates differential methylation for multiple covariates
#' from a methylBase object, and returns a multiDiff array.
#' @param meth methBase object
#' @param cur_design  the design matrix
#' @param cur_formula the design formula, including interactions
#' @param num.cores number of cores if running in parallel.
#' @param .parallel whether to run with multiple cores or not.
#' @param stat.method Which statistical test to use, Wald, or LLT. Defaults to Wald.
#' @param diff.method which method to use for effect size calculation . Defaults to "maximum"
#' @keywords multiDiff
#' @export
#' @examples
#' @calculateMultiDiffMeth()

calculateMultiDiffMeth<- function(meth, cur_design, cur_formula ,
                                  num.cores=2, .parallel=T,
                                  stat.method='Wald',
                                  diff.method='max' ){
  ## meth should be a methylKit methylBase object
  ## cur_design a binary design matrix with colnames named with the model coefficients
  ## diff.method is the method used to calculate the differential methylation - either mean or max

  ##Output is a list with the position data, an array with the output for each term, a list of the coeffiecients
  #, and the metadata for regenerating the methylDiff object

  #Default is max difference method

  require(checkmate)

  assertChoice(diff.method, c('max', 'mean'))
  assertChoice(stat.method, c('Wald','LLT'))


  # meth=setNAstoZero(meth)

  #cat('test','\n')
  pos.dat=getData(meth)[ , c('chr' , 'start', 'end', 'strand')]

  #Add one to allow for intercept
  coef_numb=length(strsplit(cur_formula, '+', fixed=T)[[1]])+1

  cur_terms=strsplit(cur_formula, '+', fixed = T)[[1]]
  if(length(cur_terms)==1){
    cur_design=data.frame(cur_design[ , cur_terms], row.names=rownames(cur_design))
    colnames(cur_design)=cur_terms
  }else {
    cur_design=cur_design[ , cur_terms]
  }

  possible_states=plyr::laply(0:(2^ncol(cur_design)-1), function (x) { as.integer((intToBits(x))  )  } )[ , 1:ncol(cur_design) ]


  #Create auxilury function, which given an index, performs the analysis at the appropriate site
  # Returns an array with p.values, methlylation difference and a slot for q-values, for each coefficient

  all_c_data=getData(meth)[   , meth@numCs.index]
  all_t_data=getData(meth)[   , meth@numTs.index]

  if (stat.method=='LLT'){
    cur_function <- function (i) {
      c_data=all_c_data[ i, ]
      t_data=all_t_data[  i , ]
      samps_with_est=!(is.na(t_data) | is.na(c_data))
      cur_data=data.frame(cur_design, t(c_data), t(t_data))
      colnames(cur_data)=c(colnames(cur_design), 'Cs', 'Ts')
      #cur_obj=glm(as.formula(paste0("cbind(Cs, Ts)~", cur_formula)), data=cur_data, family=binomial(link=logit) )

      #Added a single smoothing count to both Cs and Ts, to allow for better computation of p-values
      cur_obj=glm(as.formula(paste0("cbind(Cs, Ts)~", cur_formula)),
                  data=within(cur_data, { Cs<-Cs+1; Ts<-Ts+1}), family=binomial(link=logit) )
      coef_numb=length(cur_obj$coefficients)

      #Calculate statistics for each coefficient
      cur_coefs=cur_obj$coefficients
      cur_coef_names=names(cur_coefs)
      cur_vars=diag(vcov(cur_obj))
      #cur_pvals= 1-pchisq(cur_coefs^2/cur_vars, df=1)

      #sub_out=array(data=NA, c(1 , 3, ncol(cur_design)), dimnames = c('pos', 'val', 'coef'))
      sub_out=array(data=NA, c(1 , 4, ncol(cur_design))) #, dimnames = list('meth.diff', 'p.value' , 'q.value', 'beta'))
      cur_pvals=summary(cur_obj)$coefficients[ , 'Pr(>|z|)']
      cur_pval_names=names(cur_pvals)
      for ( j in 2:coef_numb){
        cur_design_col=as.integer(cur_design[ , j-1 ] )
        cur_beta <- cur_coefs[[j]]

        cur_rest= cur_coefs [ c(-j) ]
        cur_rest = cur_rest [ c(-1)]
        if (length(cur_rest)==1){
          cur_rest=matrix(cur_rest)
        }
        #cur_se <- summary(cur_obj)$coefficients[ j, 2 ]

        #cur_p.value=cur_pvals[[j]]
        if (!(cur_coef_names[[j]] %in% cur_pval_names)){
          cur_p.value=NA
        } else{
          cur_p.value<- cur_pvals[[j]]
        }

        if (diff.method=='max'){
          cur_push_meth=ilogit((possible_states[ , c(-(j-1))] %*% cur_rest) + cur_coefs[[1]]+cur_beta)
          cur_base_meth=ilogit((possible_states[ , c(-(j-1)) ] %*% cur_rest) + cur_coefs[[1]])

          cur_meth_diff=(cur_push_meth-cur_base_meth)

          cur_max_diff=max(cur_meth_diff)
          cur_min_diff=min(cur_meth_diff)

          ##Handle NA case
          if (is.na(cur_max_diff)){
            cur_meth_diff= cur_min_diff
          }else if (is.na(cur_min_diff)){
            cur_meth_diff= cur_max_diff

          }else if (abs(cur_max_diff)>=abs(cur_min_diff)){
            cur_meth_diff= cur_max_diff
          } else {
            cur_meth_diff= cur_min_diff
          }

          cur_meth_diff=100*cur_meth_diff

        } else if (diff.method=='mean'){
          cur_linear_predictor_est=(cur_obj$linear.predictors+cur_beta*(1-2*cur_design_col)[ samps_with_est])
          cur_meth_diff=100*(2*cur_design_col-1)*(cur_obj$fitted.values-ilogit( cur_linear_predictor_est))
          cur_meth_diff=mean(cur_meth_diff)
        }
        #cur_linear_predictor_est=(cur_obj$linear.predictors+cur_beta*(1-2*cur_design_col)[ samps_with_est])
        #cur_meth_diff=100*(2*cur_design_col-1)*(cur_obj$fitted.values-ilogit( cur_linear_predictor_est))
        #cur_meth_diff=mean(cur_meth_diff)
        #glm.dat[ i ,1:2 , j-1 ] = c(cur_meth_diff,cur_p.value)

        #sub_out[ 1 ,1:2 , j-1 ]= c(cur_meth_diff,cur_p.value)
        sub_out[ 1 , , j-1 ]= c(cur_meth_diff,cur_p.value, 0 , cur_beta)

        #return(c(cur_meth_diff,cur_p.value))


      }
      return(sub_out)
    }
  } else {
    cur_function <- function (i) {
      c_data=all_c_data[ i, ]
      t_data=all_t_data[  i , ]
      samps_with_est=!(is.na(t_data) | is.na(c_data))
      cur_data=data.frame(cur_design, t(c_data), t(t_data))
      colnames(cur_data)=c(colnames(cur_design), 'Cs', 'Ts')
      #cur_obj=glm(as.formula(paste0("cbind(Cs, Ts)~", cur_formula)), data=cur_data, family=binomial(link=logit) )

      #Added a single smoothing count to both Cs and Ts, to allow for better computation of p-values
      cur_obj=glm(as.formula(paste0("cbind(Cs, Ts)~", cur_formula)),
                  data=within(cur_data, { Cs<-Cs+1; Ts<-Ts+1}), family=binomial(link=logit) )
      coef_numb=length(cur_obj$coefficients)

      #Calculate statistics for each coefficient
      cur_coefs=cur_obj$coefficients
      cur_coef_names=names(cur_coefs)
      cur_vars=diag(vcov(cur_obj))
      #cur_pvals= 1-pchisq(cur_coefs^2/cur_vars, df=1)

      #sub_out=array(data=NA, c(1 , 3, ncol(cur_design)), dimnames = c('pos', 'val', 'coef'))
      sub_out=array(data=NA, c(1 , 4, ncol(cur_design))) #, dimnames = list('meth.diff', 'p.value' , 'q.value', 'beta'))
      cur_pvals=summary(cur_obj)$coefficients[ , 'Pr(>|z|)']
      cur_pval_names=names(cur_pvals)
      for ( j in 2:coef_numb){
        cur_design_col=as.integer(cur_design[ , j-1 ] )
        cur_beta <- cur_coefs[[j]]

        cur_rest= cur_coefs [ c(-j) ]
        cur_rest = cur_rest [ c(-1)]
        if (length(cur_rest)==1){
          cur_rest=matrix(cur_rest)
        }
        #cur_se <- summary(cur_obj)$coefficients[ j, 2 ]

        #cur_p.value=cur_pvals[[j]]
        if (!(cur_coef_names[[j]] %in% cur_pval_names)){
          cur_p.value=NA
        } else{
          cur_p.value<- cur_pvals[[j]]
        }
        #cur_p.value<- 1- pchisq((cur_beta/cur_se)^2, 1)
        #cat(cur_beta, cur_se, cur_p.value, 'n')
        #cur_p.value<- 1- pnorm(abs(cur_beta/cur_se))

        #p.value <- 1-pchisq(deviance,df=1)
        #meth_diff=100*(ilogit(sum(obj$coefficients)-beta_last)-ilogit(sum(obj$coefficients)))
        #Samples with NAs won't have estimates, so need to remove them in estimating the methylation difference

        #Calculate differential methylation using selected method
        if (diff.method=='max'){
          cur_push_meth=ilogit((possible_states[ , c(-(j-1))] %*% cur_rest) + cur_coefs[[1]]+cur_beta)
          cur_base_meth=ilogit((possible_states[ , c(-(j-1)) ] %*% cur_rest) + cur_coefs[[1]])

          cur_meth_diff=(cur_push_meth-cur_base_meth)

          cur_max_diff=max(cur_meth_diff)
          cur_min_diff=min(cur_meth_diff)

          ##Handle NA case
          if (is.na(cur_max_diff)){
            cur_meth_diff= cur_min_diff
          }else if (is.na(cur_min_diff)){
            cur_meth_diff= cur_max_diff

          }else if (abs(cur_max_diff)>=abs(cur_min_diff)){
            cur_meth_diff= cur_max_diff
          } else {
            cur_meth_diff= cur_min_diff
          }

          cur_meth_diff=100*cur_meth_diff

        } else if (diff.method=='mean'){
          cur_linear_predictor_est=(cur_obj$linear.predictors+cur_beta*(1-2*cur_design_col)[ samps_with_est])
          cur_meth_diff=100*(2*cur_design_col-1)*(cur_obj$fitted.values-ilogit( cur_linear_predictor_est))
          cur_meth_diff=mean(cur_meth_diff)
        }
        #cur_linear_predictor_est=(cur_obj$linear.predictors+cur_beta*(1-2*cur_design_col)[ samps_with_est])
        #cur_meth_diff=100*(2*cur_design_col-1)*(cur_obj$fitted.values-ilogit( cur_linear_predictor_est))
        #cur_meth_diff=mean(cur_meth_diff)
        #glm.dat[ i ,1:2 , j-1 ] = c(cur_meth_diff,cur_p.value)

        #sub_out[ 1 ,1:2 , j-1 ]= c(cur_meth_diff,cur_p.value)
        sub_out[ 1 , , j-1 ]= c(cur_meth_diff,cur_p.value, 0 , cur_beta)

        #return(c(cur_meth_diff,cur_p.value))


      }
      return(sub_out)
    }
  }


  ## Actually calculate differential ####
  #Use laply function from plyr to apply the auxilury function over each index,

  #and manage the output into a single large array

  if (length(cur_terms)==1){

    glm.dat= array(t(simplify2array(parallel::mclapply(1:nrow(meth),
                                                       cur_function ,
                                                       mc.cores=num.cores))[ 1, , , ]),
                   dim=c(nrow(meth), 3, 1))
    glm.dat[ , 2:3,  1]  = fix.q.values.glm( glm.dat[ , 2, 1],slim=FALSE)

  } else {

    glm.dat=(aperm(simplify2array(parallel::mclapply(1:nrow(meth), cur_function , mc.cores=num.cores))[ 1, , , ], c(3,1,2)))
    #cat(dim(glm.dat), '\n')
    for ( j in 2:coef_numb){

      glm.dat[ , 2:3, j-1]  = fix.q.values.glm( glm.dat[ , 2, j-1],slim=FALSE)
    }
  }
  #glm.dat=array(data=laply(1:nrow(meth), cur_function , .parallel = .parallel), c(nrow(meth), 3, ncol(cur_design)), dimnames = c('pos', 'val', 'coef'))
  #glm.dat=array(data=parallel::mclapply(1:nrow(meth), cur_function , mc.cores=3), c(nrow(meth), 3, ncol(cur_design)))
  #glm.dat[ , , ]=laply(1:nrow(meth), cur_function )


  #Calculate and fill in q-values for each index

  #Create data structure, allowing user to index using coeffieicient names
  coef.index=1:(coef_numb-1)
  names(coef.index)=colnames(cur_design)

  #Need the slot data so can rec-create standard meth.Diff objects
  slot_data=list( sample.ids=meth@sample.ids,assembly=meth@assembly
                  ,context=meth@context,
                  treatment=meth@treatment,destranded=meth@destranded
                  ,resolution=meth@resolution)

  #Final output as list
  output=list(pos.dat, glm.dat, coef.index, slot_data)

}


fix.q.values.glm<-function(pvals,slim=FALSE)
{
  if(slim==FALSE){
    #qvals=qvalue::qvalue(pvals[,3])$qvalues # get qvalues
    qvals=p.adjust(pvals,"BH")
  }else{
    slimObj=SLIMfunc(pvals)
    qvals=QValuesfun(pvals, slimObj$pi0_Est)
  }

  pvals=cbind(pvals,qvalue=qvals) # merge pvals and qvals
  return(pvals)
}

QValuesfun <-function(rawp,pi0)
{
  order_rawp=sort(rawp);
  #order_rawp=sort(rawp,method="qu");
  qvalues=pi0*length(order_rawp)*order_rawp/c(1:length(order_rawp));
  temp=cummin(qvalues[seq(length(qvalues),1,-1)])
  qvalues=temp[seq(length(temp),1,-1)];
  qvalues=qvalues[order(order(rawp))]
}

SLIMfunc<-function(rawp,STA=.1,Divi=10,Pz=0.05,B=100,Bplot=FALSE)
{


  ####################
  m <- length(rawp)

  ########################
  alpha_mtx=NULL;#
  pi0s_est_COM=NULL;
  SzCluster_mtx=NULL;
  P_pi1_mtx=NULL;
  pi1_act_mtx=NULL;
  Dist_group=NULL;
  Num_mtx=NULL;
  Gamma=NULL;
  PI0=NULL;

  #############
  ##observed points
  lambda_ga=seq(0,1,0.001);
  gamma_ga=sapply(lambda_ga,f1,rawp=rawp);
  Gamma=c(Gamma,gamma_ga);
  alpha_mtx=c(alpha_mtx,gamma_ga[which(lambda_ga==0.05)]);


  ###esimation
  pi0_mtx=NULL;
  x.axis=NULL;
  y.axis=NULL;
  itv=(1-STA)/Divi;
  for (i in 1:Divi)##10 for uniform data
  {
    cutoff=STA+(i/Divi)*(1-STA);##10 for uniform data
    lambda=seq(cutoff-itv,cutoff,itv/10);
    gamma_mtx=sapply(lambda,f1,rawp=rawp);
    LModel=lm(gamma_mtx~lambda);
    pi0_mtx=c(pi0_mtx,coefficients(LModel)[2]);
  }


  ##################################
  ########searching
  N_COM=NULL;
  N_rawp=NULL;
  maxFDR_mtx=NULL;
  quapoint_mtx=NULL;
  if (B<=1) B=100;
  quapoint_mtx=seq(0.01,0.99,1/B);
  for (k in 1:length(quapoint_mtx))
  {
    qua_point=quapoint_mtx[k];

    pi0_combLR=min(quantile(pi0_mtx,qua_point),1);#mean(pi0_mtx);#median();# qua_point=0.78 for desreasing distribution;
    ##0.4 for uniform or normal distribution;
    pi0_est=pi0_combLR;

    ###########Calculate independent index of raw p vlaues
    PI0=rbind(PI0,pi0_mtx);

    pi0s_est_COM=c(pi0s_est_COM,pi0_est);
    ##Condition1
    P_pi1=sort(rawp)[max(length(rawp)*(1-pi0_est),1)];##
    P_pi1_mtx=c(P_pi1_mtx,P_pi1);

    pi0=pi0_est;
    if (is.null(Pz)) Pz=0.05;
    maxFDR=Pz*pi0/(1-(1-Pz)*pi0);
    maxFDR_mtx=c(maxFDR_mtx,maxFDR);

    qvalues_combLR=QValuesfun(rawp,pi0);
    qvalue_cf=maxFDR;
    selected=which(qvalues_combLR<qvalue_cf);
    Sel_qvalues_combLR=selected;

    pi1_act_mtx=c(pi1_act_mtx,length(Sel_qvalues_combLR)/length(rawp));
    N_COM=c(N_COM,list(Sel_qvalues_combLR));
    Num_mtx=c(Num_mtx,length(Sel_qvalues_combLR));
  }
  length(N_COM)
  length(quapoint_mtx)

  ####doing judging
  ##by max FDR
  pi1s_est_COM=1-pi0s_est_COM;
  Diff=sum(rawp<=Pz)/length(rawp)-pi1_act_mtx;

  ###
  loc=which.min(abs(Diff));
  Diff.loc=Diff[loc];
  selQuantile=quapoint_mtx[loc];
  pi0_Est=min(1,pi0s_est_COM[loc]);
  maxFDR.Pz=Pz*pi0_Est/(1-(1-Pz)*pi0_Est);

  if(Bplot)
  {
    #windows();
    par(mfrow=c(1,2));
    hist(rawp,main="Histogram of p-value");
    gamma_ga=sapply(lambda_ga,f1,rawp=rawp);
    plot(lambda_ga,gamma_ga,type="l",main="Relationship of p- and q value",xlab=expression(lambda),ylab=expression(gamma),cex.lab=1.45,cex.axis=1.42)
    #par(xaxp=c(0,1,10));
    #axis(1);
    ##qvalues
    qValues=QValuesfun(rawp,pi0=pi0_Est);
    gammaq_ga=sapply(lambda_ga,f1,rawp=qValues);
    lines(lambda_ga,gammaq_ga,col="blue",lwd=2,lty="dashed")
    abline(v=Pz,col="black",lwd=2,lty="dotdash")
    abline(v=maxFDR.Pz,col="blue",lwd=2,lty="dotdash")
    text(0.75,0.6,labels=paste("L=",round(abs(Diff.loc),4),sep=""));
    leg=list(bquote("CPD of p-value"),bquote("CPD of q-value"),bquote("Pmax"==.(Pz)),bquote("FDRmax"==.(round(maxFDR.Pz,2))));
    legend("bottomright",legend=as.expression(leg),lwd=2,lty=c("solid","dashed","dotdash","dotdash"),col=c("black","blue","black","blue"));
  }

  return(list(pi0_Est=pi0_Est,selQuantile=selQuantile));
}

f1<-function(cutoff,rawp){sum(rawp<cutoff)/length(rawp)};

ilogit <- function(x) { return (1/(1+exp(-x))) }
logic <-function(x){ return (log(x/(1-x)))}

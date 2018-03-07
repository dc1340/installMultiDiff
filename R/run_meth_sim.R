#' Run methylation simulation
#'
#' Runs multivariate simulation of methylation for methylKit and DSS
#' @param sim_num_sites_per_cond Number of sites per condition. Default=1000
#' @param sim_read_depth Read depth- simulated uniformly. Default=10
#' @param seed Random seed. Default=1
#' @param sim_control_region_size_multiplier Allows the control region to be larger to
#'     simulate that most of the genome should not be differential. Default=1
#' @param sim_use_sawtooth_base Instead of a flat baseline, simulates a sawtooth patten. Default=F.
#' @param sim_mean_beta Mean linear effect of the covariates. Default=3
#' @param sim_beta_multiplier=c(1,-1,-1) Vector of 3 numbers allowing each covariate
#'     and their interaction term to have different strengths. Default =c(1,-1,-1).
#' @param sim_condition_noise Controls how much variance there is in the effect size across each region. Default=0.01
#' @param sim_biological_noise Simulates the biological noise across replicates. Default=0.01
#' @param run_DSS Whether to run DSS or not. Default =T.
#' @param link_func The link function. Expects inputs in [0,1]. Defaults to logit.
#' @param inv_link_func The inverse link function. Expected to output in [0,1], should be inverse of the link.
#' @param sim_cores Number of cores to use in the simulation
run_meth_sim <-function ( sim_num_sites_per_cond=1000, sim_read_depth=10,
              seed=1,
              sim_control_region_size_multiplier=1,
              sim_use_sawtooth_base=F,
              sim_mean_beta=3, sim_beta_multiplier=c(1,-1,-1),
              sim_condition_noise=0.01,
              sim_biological_noise=0.01,
              run_DSS=T, make_confusion_matrices=T,
              link_func=logit, inv_link_func=inv_link_func,
              sim_cores=4){


  set.seed(seed)

  sim_bin_names=laply (0:7 , function(x) { paste0(as.character(as.integer(intToBits(x)))[ 1:3 ], collapse='') } )

  sim_simp_design=data.frame(Cov1=c(rep(1,8), rep(0,8)) , Cov2=c(rep(c(rep(1,4), rep(0,4)), 2)))

  sim_simp_design$IsBoth=with(sim_simp_design, Cov1* Cov2)

  sim_rev_design=-sim_simp_design+1
  colnames(sim_rev_design)=c('NotCov1', 'NotCov2', 'Neither')
  sim_rev_design$Neither=with(sim_rev_design, NotCov1 * NotCov2)

  sim_mixed_design=cbind(sim_rev_design[ , c(1,2)], sim_simp_design[  , 3])

  sim_possible_site_states=unique(sim_simp_design)
  #sim_possible_site_states
  sim_possible_site_states=rbind(sim_possible_site_states,
                                 c(1,1,0), c(1,0,1), c(0,1,1), c(0,0,1))

  #Choose generative deisgn, and the design used for differential
  sim_design=sim_simp_design
  sim_diff_design=sim_simp_design


  sim_num_sites_for_control_region=sim_num_sites_per_cond*sim_control_region_size_multiplier

  sim_num_sites=(nrow(sim_possible_site_states)-1)*sim_num_sites_per_cond+
    sim_num_sites_for_control_region



  sim_control_region_state_index=which(rowSums(sim_possible_site_states)==0)[[1]]


  if (sim_use_sawtooth_base){

    sim_base_linear_meth=link_func(
      c(rep.int(seq(0.01, 0.99, length.out = sim_num_sites_per_cond),
                sim_control_region_state_index-1),
        seq(0.01, 0.99, length.out = sim_num_sites_for_control_region),
        rep.int(seq(0.01, 0.99, length.out = sim_num_sites_per_cond),
                nrow(sim_possible_site_states)-sim_control_region_state_index)
      ))

    sim_base_linear_meth=matrix(rep(sim_base_linear_meth, sim_num_samps),
                                nrow=sim_num_sites,
                                ncol=sim_num_samps)
  } else{
    sim_base_meth=0.5
    sim_base_linear_meth=link_func(sim_base_meth)

  }

  sim_num_samps=nrow(sim_design)
  #Prepare data matrix
  sim_dat=data.frame(matrix(nrow=sim_num_sites, ncol=4+3*sim_num_samps))
  sim_cov_cols=4+seq(1, by=3, length.out = sim_num_samps)
  sim_c_cols=4+seq(2, by=3, length.out = sim_num_samps)
  sim_t_cols=4+seq(3, by=3, length.out = sim_num_samps)

  base_cols=c('chr', 'start', 'end', 'strand')
  colnames(sim_dat)[ 1:4]=base_cols
  sim_dat$chr=c(rep(paste0('chr1'), sim_num_sites_for_control_region),
                rep(paste0('chr', 2:nrow(sim_possible_site_states)), each=sim_num_sites_per_cond))
  sim_dat$start=1:sim_num_sites
  sim_dat$end=1:sim_num_sites
  sim_dat$strand='+'

  sim_dat[ , sim_cov_cols]=sim_read_depth

  colnames(sim_dat)[ sim_cov_cols]=paste0('coverage' , 1:sim_num_samps)
  colnames(sim_dat)[ sim_c_cols]=paste0('numCs' , 1:sim_num_samps)
  colnames(sim_dat)[ sim_t_cols]=paste0('numTs' , 1:sim_num_samps)

  sim_site_states=Matrix( nrow=sim_num_sites, ncol=ncol(sim_design))

  #Handle Control Region
  i=sim_control_region_state_index
  sim_start=(i-1)*sim_num_sites_per_cond+1
  sim_end=(i-1)*sim_num_sites_per_cond+sim_num_sites_for_control_region
  sim_cur_dat=Matrix(unlist(rep(sim_possible_site_states[ i , ],
                                sim_num_sites_for_control_region)),
                     nrow=sim_num_sites_for_control_region,
                     ncol=ncol(sim_design),
                     byrow = T)

  sim_site_states[ sim_start:sim_end , ]=sim_cur_dat

  #Handle other regions
  for ( i in 1:(sim_control_region_state_index-1)){
    sim_start=(i-1)*sim_num_sites_per_cond+1
    sim_end=(i-1)*sim_num_sites_per_cond+sim_num_sites_per_cond
    sim_cur_dat=Matrix(unlist(rep(sim_possible_site_states[ i , ],
                                  sim_num_sites_per_cond)),
                       nrow=sim_num_sites_per_cond, ncol=ncol(sim_design), byrow = T)

    sim_site_states[ sim_start:sim_end , ]=sim_cur_dat
  }

  for ( i in (sim_control_region_state_index+1):nrow(sim_possible_site_states)){
    sim_start=sim_num_sites_for_control_region+((i-2)*sim_num_sites_per_cond)+1
    sim_end=sim_num_sites_for_control_region+((i-2)*sim_num_sites_per_cond)+
      sim_num_sites_per_cond
    sim_cur_dat=Matrix(unlist(rep(sim_possible_site_states[ i , ],
                                  sim_num_sites_per_cond)),
                       nrow=sim_num_sites_per_cond, ncol=ncol(sim_design), byrow = T)

    sim_site_states[ sim_start:sim_end , ]=sim_cur_dat
  }
  sim_bin_site_state=rowSums(matrix(as.numeric(sim_site_states),
                                    nrow=sim_num_sites) %*% diag(2^c(0:(ncol(sim_design)-1))))+1

  sim_beta_multipler_matrix=matrix((rep(sim_beta_multiplier, sim_num_sites)),
                                   nrow=sim_num_sites,
                                   ncol=ncol(sim_design),
                                   byrow=T)
  sim_raw_betas=Matrix(rnorm(sim_num_sites,
                             sim_mean_beta, sim_condition_noise ),
                       nrow=sim_num_sites, ncol=ncol(sim_design))



  sim_betas=sim_raw_betas *sim_beta_multipler_matrix* sim_site_states
  sim_linear_predictors=sim_base_linear_meth+(sim_betas %*% t(sim_design))+
    matrix(rnorm(sim_num_sites*sim_num_samps,
                 0,
                 sim_biological_noise),
           nrow=sim_num_sites,
           ncol=sim_num_samps)

  for (i in 1:sim_num_samps){

    sim_cur_c_count=aaply(inv_link_func(sim_linear_predictors[ , i]),
                          1,
                          function (x) {
                            sum(rbinom(n=sim_read_depth, prob = x, size=1))
                          }
    )
    sim_cur_t_count=sim_read_depth-sim_cur_c_count
    sim_dat[ , sim_c_cols[[i]]]=sim_cur_c_count
    sim_dat[ , sim_t_cols[[i]]]=sim_cur_t_count
  }

  sim_sample.ids=paste0('sample',1:sim_num_samps)

  sim_meth_base=new("methylBase",sim_dat,
                    sample.ids=sim_sample.ids,assembly='mm10'
                    ,context='CpG',
                    treatment=rep(c(0,1), sim_num_samps/2),
                    coverage.index=sim_cov_cols,
                    numCs.index=sim_c_cols,
                    numTs.index=sim_t_cols,
                    destranded=T,
                    resolution='base')


  ### Calculate differential methylation ####
  # sim_diff=cem.calculateMultiDiffMeth(sim_meth_base, sim_diff_design, paste(colnames(sim_diff_design), collapse = '+') , num.cores = sim_cores);
  # sim_dc_diff=dc.calculateMultiDiffMeth(sim_meth_base, sim_diff_design, paste(colnames(sim_diff_design), collapse = '+') , num.cores = sim_cores);
  sim_jm_diff=calculateMultiDiffMeth(sim_meth_base,
                                     sim_diff_design,
                                     paste(colnames(sim_diff_design),
                                           collapse = '+') ,
                                     num.cores = sim_cores);

  if (run_DSS){
    sim_multiDSS = DMLfit.multiFactor(convertMethBaseToBSeqData(sim_meth_base),
                                      design=sim_diff_design,
                                      formula=~Cov1+Cov2+IsBoth)
  }
  rm(sim_base_meth)

  # sim_sig_diff=apply(sim_diff[[2]][ , 1, ], 2, mean)
  # sim_sig_diff=sim_sig_diff/(sim_sig_diff[[1]])

  ### Generate confusion matrices ####
  if (make_confusion_matrices){
    ##Ideal Confusion
    sim_ideal_confusion_matrix=matrix(rep(0, 64), nrow=8, ncol=8)
    diag(sim_ideal_confusion_matrix)=rep(sim_num_sites_per_cond, 8)
    colnames(sim_ideal_confusion_matrix)=sim_bin_names
    sim_ideal_confusion_matrix[ which(sim_bin_names=="000"),
                                which(sim_bin_names=="000")]=sim_num_sites_for_control_region
    ##Confusion for multiDiff
    sim_jm_diff_matrix=getDiffMatrix(sim_jm_diff)
    sim_jm_diff_call=getDiffCallMatrix(sim_jm_diff)

    sim_jm_call_state=rowSums(matrix(as.numeric(sim_jm_diff_call),
                                     nrow=sim_num_sites) %*% diag(2^c(0:(ncol(sim_design)-1))))+1

    sim_jm_confusion_matrix=matrix(rep(0, 64), nrow=8, ncol=8)
    for (i in 1:length(sim_bin_site_state)){
      sim_jm_confusion_matrix[ sim_bin_site_state[[i]], sim_jm_call_state[[i]]]=
        sim_jm_confusion_matrix[ sim_bin_site_state[[i]], sim_jm_call_state[[i]]]+1
    }

    ##Confusion for multiDSS
    sim_multiDSS_diff_matrix=getMaxDiffFromMultiDSS(sim_multiDSS)
    sim_multiDSS_diff_call=getMultiDSSDiffCall(sim_multiDSS)

    sim_multiDSS_call_state=rowSums(matrix(as.numeric(sim_multiDSS_diff_call),
                                           nrow=sim_num_sites) %*% diag(2^c(0:(ncol(sim_design)-1))))+1

    sim_multiDSS_confusion_matrix=matrix(rep(0, 64), nrow=8, ncol=8)
    for (i in 1:length(sim_bin_site_state)){
      sim_multiDSS_confusion_matrix[ sim_bin_site_state[[i]], sim_multiDSS_call_state[[i]]]=
        sim_multiDSS_confusion_matrix[ sim_bin_site_state[[i]], sim_multiDSS_call_state[[i]]]+1
    }

    ##Confusion by q-value for multiDiff

    sim_jm_diff_call=getDiffCallMatrix(sim_jm_diff, meth.thresh = 0)

    sim_jm_q_confusion_matrix=matrix(rep(0, 64), nrow=8, ncol=8)
    for (i in 1:length(sim_bin_site_state)){
      sim_jm_q_confusion_matrix[ sim_bin_site_state[[i]], sim_jm_call_state[[i]]]=
        sim_jm_q_confusion_matrix[ sim_bin_site_state[[i]], sim_jm_call_state[[i]]]+1
    }


    ##Confusion by q-value for DSS


    sim_multiDSS_diff_call=getMultiDSSDiffCall(sim_multiDSS, meth.thresh =0)


    sim_multiDSS_q_confusion_matrix=matrix(rep(0, 64), nrow=8, ncol=8)
    for (i in 1:length(sim_bin_site_state)){
      sim_multiDSS_q_confusion_matrix[ sim_bin_site_state[[i]], sim_multiDSS_call_state[[i]]]=
        sim_multiDSS_q_confusion_matrix[ sim_bin_site_state[[i]], sim_multiDSS_call_state[[i]]]+1
    }


    # rownames(sim_dc_confusion_matrix)=sim_bin_names
    # rownames(sim_jm_confusion_matrix)=sim_bin_names
    colnames(sim_jm_confusion_matrix)=sim_bin_names
    colnames(sim_multiDSS_confusion_matrix)=sim_bin_names
    colnames(sim_multiDSS_q_confusion_matrix)=sim_bin_names

    output=
  }

 return (output)
}

#cat('Input: ', sim_beta_multiplier, '\n', 'Output 1: ', sim_sig_diff, '\n', 'Output 2: ', sim_sig_dc, '\n')

#makeHeatMap(sim_dc_diff, cluster_cols = F, plot=T)

#sim_rev_diff=cem.calculateMultiDiffMeth(sim_meth_base, sim_mixed_design, paste(colnames(sim_mixed_design), collapse = '+') , num.cores = sim_cores);
#makeHeatMap(sim_rev_diff, cluster_cols = F)

#
# sim_diff_matrix=getDiffMatrix(sim_diff)
# sim_diff_call=getDiffCallMatrix(sim_diff)
#
# sim_call_state=rowSums(matrix(as.numeric(sim_diff_call),
#                               nrow=sim_num_sites) %*% diag(2^c(0:(ncol(sim_design)-1))))+1


#
# sim_confusion_matrix=matrix(rep(0, 64), nrow=8, ncol=8)
# for (i in 1:length(sim_bin_site_state)){
#   sim_confusion_matrix[ sim_bin_site_state[[i]], sim_call_state[[i]]]=
#     sim_confusion_matrix[ sim_bin_site_state[[i]], sim_call_state[[i]]]+1
# }

##Confusion for dc
#
# sim_dc_diff_matrix=getDiffMatrix(sim_dc_diff)
# sim_dc_diff_call=getDiffCallMatrix(sim_dc_diff)
#
# sim_dc_call_state=rowSums(matrix(as.numeric(sim_dc_diff_call),
#                                  nrow=sim_num_sites) %*% diag(2^c(0:(ncol(sim_design)-1))))+1
#
# sim_dc_confusion_matrix=matrix(rep(0, 64), nrow=8, ncol=8)
# for (i in 1:length(sim_bin_site_state)){
#   sim_dc_confusion_matrix[ sim_bin_site_state[[i]], sim_dc_call_state[[i]]]=
#     sim_dc_confusion_matrix[ sim_bin_site_state[[i]], sim_dc_call_state[[i]]]+1
# }


#colnames(sim_dc_confusion_matrix)=sim_bin_names
# pheatmap(sim_multiDSS_confusion_matrix, cluster_rows = F, cluster_cols = F, display_numbers = T)
# pheatmap(sim_jm_confusion_matrix, cluster_rows = F, cluster_cols = F, display_numbers = T)

# sum(diag(sim_multiDSS_confusion_matrix))/sim_num_sites
# sum(diag(sim_jm_confusion_matrix))/sim_num_sites

#sim_jm_confusion_matrix- sim_multiDSS_confusion_matrix
#sim_roc=getRocDat(sim_dc_diff, sim_site_states)
#png('sim_meth_heatmap.png' , width = 480, height=480)
# pheatmap(t(getMethMatrix(sim_meth_base)),
#          cluster_rows=F, culster_cols=F,
#          show_rownames = F, show_colnames = F
#         )
#
# dev.off()

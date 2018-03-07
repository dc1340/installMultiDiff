# installMultiDiff
Quick multiDiff installation. multiDiff is an R package designed to extend methyKit for use with multiple covariates. It implements the maximum difference estimate to allow the ability to assign biologically meaningful effect sizes in complex designs. 


methyKit, created by Dr. Altuna Akalin and Dr. Sheng Li, is on github [here](https://github.com/al2na/methylKit)

# Features
* Calling and filtering sites
* Visualizing and summarizing results
* Running simulations of two covariates and their interaction.
* Convert back to methylKit objects for use in established workflows



# Installation

You can install multiDiff from github:

```
library(devtools)
install_github("dc1340/installMultiDiff", build_vignettes=FALSE,
  dependencies=TRUE)
```

# Test

You can run a quick test of the simulation and visualization with:

```
library(installMultiDiff)
testMultiDiff=run_meth_sim(sim_num_sites_per_cond = 50, run_DSS = F, make_confusion_matrices = F)
makeHeatMap(testMultiDiff, cluster_cols = F)
makeViolinPlot(testMultiDiff)
```

# installMultiDiff
Quick multiDiff installation. multiDiff is an R package designed to extend methyKit for use with multiple covariates, with the ability to assign biological effect sizes in complex designs. Using it you can:

* Call and filter sites
* Visualize and summarize results
* Run simulations
* Convert back to methylKit objects for use in established workflows

methyKit, created by Dr. Altuna Akalin is on github [here](https://github.com/al2na/methylKit)

# Installation

You can install multiDiff from github:

```
library(devtools)
install_github(dc1340/installMultiDiff", build_vignettes=FALSE,
  dependencies=TRUE)
```

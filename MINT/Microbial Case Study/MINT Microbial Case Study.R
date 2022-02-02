## ----setup, include=FALSE----------------------------------------------------------------------
knitr::opts_chunk$set(dpi = 100, 
                      echo= TRUE, 
                      warning=FALSE, 
                      message=FALSE, 
                      fig.show=TRUE, 
                      fig.keep = 'all',
                      out.width = "90%")


## ----------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library


## ----------------------------------------------------------------------------------------------
load("Microbial Data/mint_phenol_ammonia.RData") # load the data

treatment <- metadata$inhib_inoc # extract the Y vector, the type of inoculant
study = metadata$experiment # extract the study each sample is from


## ----------------------------------------------------------------------------------------------
TSS.divide = function(x){ # function to apply a TSS transformation to the data
 x/sum(x)
}

# convert from raw to compositional data as well as transpose abundance so
# that OTUs are columns 
t.abundance <- apply(t(abundance), 1, TSS.divide)


## ----------------------------------------------------------------------------------------------
abundance.offset <- t.abundance + 0.01 # apply offset

# can see there is a lot of zeroes
cat("Number of zeroes prior to offset: ", length(which(t.abundance==0))) 

cat("Number of zeroes after offset: ", length(which(abundance.offset==0)))


## ----------------------------------------------------------------------------------------------
# function to remove any OTUs which has a count lower than cutoff % of the total
low.count.removal <- function(data, cutoff=0.01) {
  
    # find which OTUs (rows) have counts above cutoff% of the total
    keep.otu = which(rowSums(data)*100/(
      sum(rowSums(data))) > 0.01)
    
    data.filter = data[keep.otu,] # filter out those with lower than cutoff %
    return(list(data.filter = data.filter, keep.otu = keep.otu))
}

result.filter <- low.count.removal(abundance.offset, cutoff=0.01)
abundance.filter <- result.filter$data.filter

cat("Number of OTU's prior to filtering: ", nrow(abundance.offset))
cat("Number of OTU's remaining after filtering: ", nrow(abundance.filter))


## ----------------------------------------------------------------------------------------------
# apply the clr transformation, offset = 0 as this was done above
# note that the abundance dataframe was transposed for the last time, leaving 
# the OTUs in the columns
abundance.processed <- logratio.transfo(t(abundance.filter),
                                  logratio = 'CLR', offset = 0)

class(abundance.processed) <- "matrix" # convert from clr object to matrix
abundance.processed <- data.frame(abundance.processed) # then convert to df


## ---- fig.align = "center"---------------------------------------------------------------------
# undergo PCA
ab.pca <- pca(abundance.processed, scale = TRUE, center = TRUE, ncomp = 5)

# plot projection of samples in PC space
plotIndiv(ab.pca, group=treatment, 
          ind.names = F,legend=T,
          pch = as.numeric(factor(study))+15,
          pch.levels=(study), ellipse = TRUE, 
          title="PCA",legend.title = "Inhibitor", 
          legend.title.pch = "Experiment", size.legend = rel(2.4),
          size.legend.title = rel(2.5))


## ---- eval = FALSE-----------------------------------------------------------------------------
## # undergo repeated CV tuning, using the maximum distance metric
## ab.splsda.tuning <- tune.splsda(abundance.processed, treatment,
##                                 ncomp = 5, test.keepX = seq(5,80, 5),
##                                 validation = c('Mfold'), dist = 'max.dist',
##                                 folds = 5, nrepeat = 10)
## 
## # To reduce run time, the tuning object is saved externally and loaded in
## #save(ab.splsda.tuning, file="Microbial Data/ab_splsda_tuning.RData")


## ---- echo = FALSE-----------------------------------------------------------------------------
load("Microbial Data/ab_splsda_tuning.RData")


## ---- fig.align = "center"---------------------------------------------------------------------
plot(ab.splsda.tuning) # plot the tuning


## ----------------------------------------------------------------------------------------------
# extract the tuned ncomp and keepX parameters
optimal.ncomp <- ab.splsda.tuning$choice.ncomp$ncomp
optimal.keepX <- ab.splsda.tuning$choice.keepX[1:optimal.ncomp]


## ----------------------------------------------------------------------------------------------
optimal.ncomp 
optimal.keepX


## ---- fig.align = "center"---------------------------------------------------------------------
# generate sPLS-DA model using tuned parameters
ab.splsda.tuned <- splsda(abundance.processed, treatment, scale = TRUE, 
                    ncomp = optimal.ncomp, keepX = optimal.keepX)

# plot projection of samples onto the sPLS-DA components
plotIndiv(ab.splsda.tuned,  
          ind.names = F,legend=T,
          pch = as.numeric(factor(study))+15,
          pch.levels=(study), ellipse = TRUE,
          title="sPLS-DA",legend.title = "Inhibitor", 
          legend.title.pch = "Experiment", 
          size.legend = rel(2.4), size.legend.title = rel(2.5))


## ----------------------------------------------------------------------------------------------
# asssess the performance of the model
splsda_perf = perf(ab.splsda.tuned, validation = 'Mfold', 
                   folds = 5, nrepeat = 20,
                   progressBar = FALSE)

splsda_perf$error.rate # print the error rate


## ---- eval = FALSE-----------------------------------------------------------------------------
## # tune the ncomp and keepX parameters for the MINT sPLS-DA model
## ab.mint.splsda.tuning <- tune(method = "mint.splsda",
##                               X=abundance.processed,
##                               treatment,
##                               study = study,
##                               ncomp = 5,
##                               test.keepX = seq(5,80, 5),
##                               measure = 'BER', # balanced error rate
##                               dist = "centroids.dist")
## 
## # To reduce run time, the tuning object is saved externally and loaded in
## #save(ab.mint.splsda.tuning, file="Microbial Data/ab_mint_splsda_tuning.RData")


## ---- echo = FALSE-----------------------------------------------------------------------------
load("Microbial Data/ab_mint_splsda_tuning.RData")


## ----------------------------------------------------------------------------------------------
# this determines the indices of the lowest error rate and extracts the 
# second dimension which corresponds to the optimal component number 
optimal.ncomp <- which(ab.mint.splsda.tuning$error.rate == 
                         min(ab.mint.splsda.tuning$error.rate),
                       arr.ind = TRUE)[2]

# extract the optimal keepX parameter
optimal.keepX <- ab.mint.splsda.tuning$choice.keepX[1:optimal.ncomp]

optimal.ncomp
optimal.keepX


## ----------------------------------------------------------------------------------------------
# form tuned MINT sPLS-DA model
ab.mint.splsda <- mint.splsda(X = abundance.processed, Y = treatment,
                              ncomp = optimal.ncomp, keepX = optimal.keepX,
                              study = study)


## ---- fig.show = "hold", out.width = "49%"-----------------------------------------------------
plotIndiv(ab.mint.splsda, ind.names = F,legend=T,
          pch = as.numeric(factor(study))+15,
          pch.levels=study,
          subtitle="Figure 4(a)",legend.title = "Inhibitor",
          legend.title.pch = "Experiment", 
          size.legend = rel(0.8))

plotIndiv(ab.mint.splsda, study = "all.partial",
          ind.names = F,legend=T,
          pch = as.numeric(factor(study))+15, pch.levels=study,
          title = "Figure 4(b)", subtitle = paste("Study",1:2), 
          size.legend = rel(0.8),legend.title = "Inhibitor")


## ---- fig.show = "hold", out.width = "49%"-----------------------------------------------------
plotLoadings(ab.mint.splsda, method = "median", comp = 1, 
             legend = T, study = "global", contrib = "max",
             title = "(a) All Studies, Comp 1")

plotLoadings(ab.mint.splsda, method = "median", comp = 1, 
             legend = F, study = "all.partial", contrib = "max",
             title = "(b) Individual Studies, Comp 1", 
             subtitle = c("Study 1", "Study 2"))


## ---- fig.align = "center"---------------------------------------------------------------------
plotLoadings(ab.mint.splsda, method = "median", comp = 2, 
             legend = T, study = "global", contrib = "max",
             title = "All Studies, Comp 2")


## ---- fig.align = "center"---------------------------------------------------------------------
plotVar(ab.mint.splsda, var.names = FALSE,
        pch = 16)


## ---- echo = FALSE-----------------------------------------------------------------------------
cim(ab.mint.splsda, 
    row.sideColors = cbind(color.mixo(as.numeric(treatment)),
                           color.mixo(as.numeric(study)+4)),
    
    legend = list(legend = cbind(c(levels(treatment)), c(levels(study))), 
                col = cbind(c(color.mixo(1:3)), c(color.mixo(5:6))),
                title = "Treatment and Study", cex = 0.8),
    
    save = 'jpeg', name.save = 'MINT Microbial CIM')


## ---- eval = FALSE-----------------------------------------------------------------------------
## cim(ab.mint.splsda,
##     row.sideColors = cbind(color.mixo(as.numeric(treatment)), color.mixo(as.numeric(study)+4)),
## 
##     legend = list(legend = cbind(c(levels(treatment)), c(levels(study))),
##                 col = cbind(c(color.mixo(1:3)), c(color.mixo(5:6))),
##                 title = "Treatment and Study", cex = 0.8)
##     )


## ---- echo = FALSE-----------------------------------------------------------------------------
network(ab.mint.splsda, cutoff = 0.6,
        color.node = c(color.mixo(1), color.mixo(2)), 
        shape.node = c("rectangle", "circle"),
        save = 'jpeg', name.save = 'MINT Microbial Network')


## ---- eval = FALSE-----------------------------------------------------------------------------
## network(ab.mint.splsda, cutoff = 0.6,
##         color.node = c(color.mixo(1), color.mixo(2)),
##         shape.node = c("rectangle", "circle"))


## ---- fig.align = "center"---------------------------------------------------------------------
ab.mint.splsda.perf <- perf(ab.mint.splsda, folds = 5, nrepeat = 10)
plot(ab.mint.splsda.perf)


## ---- eval = FALSE-----------------------------------------------------------------------------
## auroc(ab.mint.splsda, roc.comp = 1, print = FALSE)
## auroc(ab.mint.splsda, roc.comp = 2, print = FALSE)
## auroc(ab.mint.splsda, roc.comp = 3, print = FALSE)


## ---- fig.show = "hold", out.width = "49%"-----------------------------------------------------
auroc(ab.mint.splsda, roc.comp = 1, roc.study = "Study 1", 
      print = FALSE, title = "(a) Study 1 AUROC")
auroc(ab.mint.splsda, roc.comp = 1, roc.study = "Study 2", 
      print = FALSE, title = "(b) Study 2 AUROC")


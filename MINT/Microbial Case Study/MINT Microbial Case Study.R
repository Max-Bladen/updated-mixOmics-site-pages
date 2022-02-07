## ----setup, include=FALSE--------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(dpi = 100, 
                      echo= TRUE, 
                      warning=FALSE, 
                      message=FALSE, 
                      fig.show=TRUE, 
                      fig.keep = 'all',
                      out.width = "90%")


## --------------------------------------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library


## --------------------------------------------------------------------------------------------------------------------------
load("Microbial Data/data_prediction_data.RData") # load the data

# extract the Y vectors
treatment.1.2 <- metadata_studies_1_2$inhib_inoc 
treatment.3.4 <- metadata_studies_3_4$inhib_inoc 
treatment <- c(treatment.1.2, treatment.3.4)

# extract the study each sample is from
study.1.2 <- metadata_studies_1_2$experiment
study.3.4 <- metadata_studies_3_4$experiment
study <- c(study.1.2,study.3.4)

# combine abundance datasets into a single dataframe
abundance <- cbind(abundance_studies_1_2, abundance_studies_3_4)


## ---- echo = FALSE---------------------------------------------------------------------------------------------------------
# function to remove any OTUs which has a count lower than cutoff % of the total
low.count.removal <- function(data, cutoff=0.01) {
  
    # find which OTUs (rows) have counts above cutoff% of the total
    keep.otu = which(rowSums(data)*100/(
      sum(rowSums(data))) > 0.01)
    
    data.filter = data[keep.otu,] # filter out those with lower than cutoff %
    return(list(data.filter = data.filter, keep.otu = keep.otu))
}

# -----------------------------------------------------------------------------#

abundance.offset <- abundance + 1 # apply offset

result.filter <- low.count.removal(abundance.offset, cutoff=0.01)
abundance.filter <- result.filter$data.filter

# apply the clr transformation, offset = 0 as this was done above
# note that the abundance dataframe was transposed for the last time, leaving 
# the OTUs in the columns
abundance.processed <- logratio.transfo(t(abundance.filter),
                                  logratio = 'CLR', offset = 0)

class(abundance.processed) <- "matrix" # convert from clr object to matrix
abundance.processed <- data.frame(abundance.processed) # then convert to df


## ---- fig.align = "center"-------------------------------------------------------------------------------------------------
# undergo PCA
ab.pca <- pca(abundance.processed, scale = TRUE, center = TRUE, ncomp = 5)

# plot projection of samples in PC space
plotIndiv(ab.pca, group=treatment, 
          ind.names = F,legend=T,
          pch = as.numeric(factor(study))+14,
          pch.levels=(study), ellipse = TRUE, 
          title="PCA",legend.title = "Inhibitor", 
          legend.title.pch = "Experiment", size.legend = rel(2.4),
          size.legend.title = rel(2.5))


## --------------------------------------------------------------------------------------------------------------------------
# determine indices of training and testing sets
train.idx <- which(study %in% c("Study 1", "Study 2"))
test.idx <- which(study %in% c("Study 3", "Study 4"))

# separate the predictors
X.train <- abundance.processed[train.idx, ]
X.test <- abundance.processed[test.idx, ]

# separate the class labels 
Y.train <- treatment[train.idx]
Y.test <- treatment[test.idx]

# separate the study associations
study.train <- as.factor(as.character(study[train.idx]))
study.test <- as.factor(as.character(study[test.idx]))


## --------------------------------------------------------------------------------------------------------------------------
# tune the ncomp and keepX parameters for the MINT sPLS-DA model
ab.mint.splsda.tuning <- tune(method = "mint.splsda",
                              X = X.train, 
                              Y = Y.train,
                              study = study.train,
                              ncomp = 5,
                              test.keepX = seq(5,50, 3),
                              measure = 'BER', # balanced error rate
                              dist = "centroids.dist")


## --------------------------------------------------------------------------------------------------------------------------
ab.mint.splsda.tuning$error.rate

optimal.ncomp <- 2

# extract the optimal keepX parameter
optimal.keepX <- ab.mint.splsda.tuning$choice.keepX[1:optimal.ncomp]

optimal.keepX


## --------------------------------------------------------------------------------------------------------------------------
# form tuned, trained MINT sPLS-DA model
ab.mint.splsda <- mint.splsda(X = X.train, Y = Y.train,
                              ncomp = optimal.ncomp, keepX = optimal.keepX,
                              study = study.train)


## ---- fig.align = "center"-------------------------------------------------------------------------------------------------
plotIndiv(ab.mint.splsda, ind.names = F,legend=T,
          pch = as.numeric(factor(study.train))+14,
          pch.levels=study.train,
          ellipse = T,
          subtitle="sPLS-DA Sample Projection",legend.title = "Inhibitor",
          legend.title.pch = "Experiment", 
          size.legend = rel(0.8))


## ---- fig.show = "hold", out.width = "49%"---------------------------------------------------------------------------------
plotLoadings(ab.mint.splsda, method = "median", comp = 1, 
             legend = T, study = "global", contrib = "max",
             title = "(a) All Studies, Comp 1")

plotLoadings(ab.mint.splsda, method = "median", comp = 2, 
             legend = T, study = "global", contrib = "max",
             title = "(b) All Studies, Comp 2")


## ---- fig.align = "center"-------------------------------------------------------------------------------------------------
plotVar(ab.mint.splsda, var.names = FALSE,
        pch = 16, cutoff = 0.5)


## ---- echo = FALSE, eval = FALSE-------------------------------------------------------------------------------------------
cim(ab.mint.splsda,
    row.sideColors = cbind(color.mixo(as.numeric(Y.train)),
                           color.mixo(as.numeric(study.train)+4)),

    legend = list(legend = cbind(c(levels(Y.train)), c(levels(study.train))),
                col = cbind(c(color.mixo(1:2)), c(color.mixo(5:6))),
                title = "Treatment and Study", cex = 0.8),

    save = 'jpeg', name.save = 'MINT Microbial CIM')


## ---- eval = FALSE---------------------------------------------------------------------------------------------------------
## cim(ab.mint.splsda,
##     row.sideColors = cbind(color.mixo(as.numeric(Y.train)),
##                            color.mixo(as.numeric(study.train)+4)),
## 
##     legend = list(legend = cbind(c(levels(Y.train)), c(levels(study.train))),
##                 col = cbind(c(color.mixo(1:2)), c(color.mixo(5:6))),
##                 title = "Treatment and Study", cex = 0.8)
##     )


## ---- echo = FALSE, eval = FALSE-------------------------------------------------------------------------------------------
names <- gsub("[^0-9.-]", "", colnames(X.train))
names[length(names)] <- "C3"

network(ab.mint.splsda, cutoff = 0.7, comp = 1,
        color.node = c(color.mixo(1), color.mixo(2)),
        shape.node = c("circle", "rectangle"),
        lty.edge = c("dotted", "solid"),
        cex.node.name = 0.7,
        alpha.node = 0.5,
        row.names = names,
        save = 'jpeg', name.save = 'MINT Microbial Network')


## ---- eval = FALSE---------------------------------------------------------------------------------------------------------
## network(ab.mint.splsda, cutoff = 0.7, comp = 1,
##         color.node = c(color.mixo(1), color.mixo(2)),
##         shape.node = c("circle", "rectangle"),
##         lty.edge = c("dotted", "solid"),
##         cex.node.name = 0.7,
##         alpha.node = 0.5,
##         row.names = names)


## ---- fig.align = "center"-------------------------------------------------------------------------------------------------
ab.mint.splsda.perf <- perf(ab.mint.splsda, folds = 5, nrepeat = 10)
plot(ab.mint.splsda.perf)


## --------------------------------------------------------------------------------------------------------------------------
# make predictions of stem cell type of study 3 samples
predict.splsda <- predict(ab.mint.splsda, newdata = X.test, 
                             dist = "max.dist", 
                             study.test = study.test)


## ---- echo = FALSE---------------------------------------------------------------------------------------------------------
YieldErrorRates <- function(comp) {
  # extract the predictions
  test.prediction <- predict.splsda$class$max.dist[, comp]

  # generate the classification confusion matrix
  conf.mat <- get.confusion_matrix(truth = Y.test, 
                                   predicted = test.prediction)
  
  cat("Metrics for model with", comp, "component(s): ", "\n")
  cat("Error rate: ", (sum(conf.mat) - sum(diag(conf.mat)))/sum(conf.mat), "\n")
  cat("Balanced error rate: ", get.BER(conf.mat), "\n")
}


## --------------------------------------------------------------------------------------------------------------------------
YieldErrorRates(1)
YieldErrorRates(2)


## ---- fig.align = "center"-------------------------------------------------------------------------------------------------
auroc(ab.mint.splsda, roc.comp = 2, print = FALSE)


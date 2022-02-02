## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', fig.height= 6, fig.width=7,
                      knitr.duplicate.label = 'allow')




## -------------------------------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library


## -------------------------------------------------------------------------------------------------------------------
data(stemcells) # extract stem cells data

X <- stemcells$gene # use genetic expression levels as predictors
Y <- stemcells$celltype # use stem cell type as response
study <- stemcells$study # extract study allocation of samples

dim(X)
summary(Y)
table(Y,study)


## ---- fig.cap = "FIGURE 1:  Sample plot from the MINT PLS-DA performed on the `stemcells` gene expression data. Samples are projected into the space spanned by the first two components. Samples are coloured by their cell types and symbols indicate the study membership."----
basic.plsda.model <- mint.plsda(X, Y, study = study, ncomp = 2) # generate basic MINT pls-da model

plotIndiv(basic.plsda.model, legend = TRUE, title = ' ', subtitle = ' ', ellipse = TRUE) # plot the samples


## -------------------------------------------------------------------------------------------------------------------
basic.splsda.model <- mint.plsda(X, Y, study = study, ncomp = 5)


## ---- fig.cap = "FIGURE 2: Choosing the number of components in `mint.splsda` using `perf()` with LOGOCV in the `stemcells` study. Classification error rates are represented on the y-axis with respect to the number of components on the x-axis for each prediction distance"----
set.seed(5249) # For reproducible results here, remove for your own analyses

splsda.perf <- perf(basic.splsda.model) # undergo performance optimisation
plot(splsda.perf)


## -------------------------------------------------------------------------------------------------------------------
splsda.perf$choice.ncomp

optimal.ncomp <- 2 

# note that usually a call in the below line would be more appropriate:
# optimal.ncomp <- splsda.perf$choice.ncomp["overall", "centroids.dist"]


## ---- fig.cap = "FIGURE 3: Tuning keepX in MINT sPLS-DA performed on the `stemcells` gene expression data. Each coloured line represents the balanced error rate (y-axis) per component across all tested keepX values (x-axis). The diamond indicates the optimal keepX value on a particular component which achieves the lowest classification error rate as determined with a one-sided t−test across the studies"----
set.seed(5249) # For a reproducible result here, remove for your own analyses

splsda.tune <- tune(X = X, Y = Y, study = study,  # tune the number of features
                              ncomp = optimal.ncomp,# using optimal comp number
                              test.keepX = seq(1, 100, 1), 
                              method = 'mint.splsda', 
                              measure = 'BER', # balanced error rate
                              dist = "centroids.dist")

plot(splsda.tune)
 
optimal.keepX <- splsda.tune$choice.keepX # extract optimal values

splsda.tune$choice.keepX


## -------------------------------------------------------------------------------------------------------------------
#generate optimal model using tuned parameters
final.splsda.model <- mint.splsda(X = X, Y = Y, study = study, 
                                  ncomp = optimal.ncomp, 
                                  keepX = optimal.keepX)


## ---- fig.show = "hold", out.width = "49%", fig.cap = "FIGURE 4: : Sample plots from the MINT sPLS-DA performed on the `stemcells` gene expression data. Samples are projected into the space spanned by the first two components. Samples are coloured by their cell types and symbols indicate study membership. (a) Global components from the model with 95% ellipse confidence intervals around each sample class. (b) Partial components per study show a good agreement across studies"----
plotIndiv(final.splsda.model, study = 'global', 
          legend = TRUE,
          title = '(a) Stem cells, MINT sPLS-DA',
          subtitle = 'Global', ellipse = T)

plotIndiv(final.splsda.model, study = 'all.partial', 
          legend = TRUE,
          title = '(b) Stem cells, MINT sPLS-DA',
          subtitle = paste("Study",1:4))


## ---- fig.cap = "FIGURE 5: Correlation circle plot representing the genes selected by MINT sPLS-DA performed on the `stemcells` gene expression data to examine the association of the genes selected on the first two components"----
# all gene names have the same first 10 characters, shorten all the names to reduce visual clutter
shortenedNames <- list(unlist(lapply(final.splsda.model$names$colnames$X, substr, start = 10, stop = 16)))

plotVar(final.splsda.model,
        cutoff = 0.5,
        var.names = shortenedNames)


## ---- fig.cap = "FIGURE 6: Clustered Image Map of the genes selected by MINT sPLS-DA on the `stemcells` gene expression data for component 1 only. A hierarchical clustering based on the gene expression levels of the selected genes on component 1, with samples in rows coloured according to cell type"----
cim(final.splsda.model, comp = 1, margins=c(10,5), 
    row.sideColors = color.mixo(as.numeric(Y)), 
    row.names = FALSE, title = "MINT sPLS-DA, component 1")


## ---- fig.cap = "FIGURE 7: Relevance network of the genes selected by MINT sPLS-DA performed on the `stemcells` gene expression data for component 1 only."----
network(final.splsda.model, comp = 1, 
        color.node = c(color.mixo(1), color.mixo(2)), 
        shape.node = c("rectangle", "circle"))


## ---- fig.cap = "FIGURE 8: ROC curve and AUC from the MINT sPLS-DA performed on the `stemcells` gene expression data for global component 1, averaged across one-vs-all comparisons. Numerical outputs include the AUC and a Wilcoxon test p−value for each ‘one vs. other’ class comparison that are performed per component."----
auroc(final.splsda.model, roc.comp = 1)


## -------------------------------------------------------------------------------------------------------------------
# determine what rows correpsond to the third study
indiv.test <- which(study == "3")

# train a model on studies 1,2 and 4
perf.splsda.model <- mint.splsda(X = X[-c(indiv.test), ], 
                                 Y = Y[-c(indiv.test)], 
                                 study = droplevels(study[-c(indiv.test)]),
                                 ncomp = optimal.ncomp,
                                 keepX = optimal.keepX)

# make predictions of stem cell type of study 3 samples
predict.splsda <- predict(perf.splsda.model, newdata = X[indiv.test, ], 
                             dist = "centroids.dist", 
                             study.test = factor(study[indiv.test]))

# store class prediction with a model of 1 component
indiv.prediction <- predict.splsda$class$centroids.dist[, 1]

# generate the classification confusion matrix
conf.mat <- get.confusion_matrix(truth = Y[indiv.test], predicted = indiv.prediction)

conf.mat

# calculate the prediction error rate
(sum(conf.mat) - sum(diag(conf.mat)))/sum(conf.mat)



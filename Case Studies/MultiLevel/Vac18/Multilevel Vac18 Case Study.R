## ----global_options, include=FALSE------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', out.width = '90%')


## ---------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library

set.seed(5249) # for reproducibility, remove for normal use


## ---------------------------------------------------------------------------------------------
data(vac18) # extract the vac18 data
X <- vac18$genes # use the genetic expression data as the X (predictor) dataframe
dim(X) # check dimensions of gene expression data

summary(vac18$stimulation) # observe distribution of response variable
summary(as.factor(vac18$sample)) # observe distribution of subjects


## ---- fig.show = "hold", out.width = "49%", fig.cap = "FIGURE 1: PCA and multilevel PCA sample plot on the gene expression data from the `vac18` study. (a) The sample plot shows that samples within the same individual tend to cluster, as indicated by the individual ID, but we observe no clustering according to treatment. (b) After multilevel PCA, we observe some clustering according to treatment type"----

pca.vac18 <- pca(X, scale = TRUE, center = TRUE) # undergo normal PCA after scaling/centering

design <- data.frame(sample = vac18$sample) # set the multilevel design

pca.multilevel.vac18 <- pca(X, scale = TRUE, center = TRUE, # undergo multilevel PCA after scaling/centering
                            multilevel = design)

plotIndiv(pca.vac18, group = vac18$stimulation, # plot the samples on normal PCs
          ind.names = vac18$sample,
          legend = TRUE, legend.title = 'Stimulation',
          title = '(a) PCA on VAC18 data')


plotIndiv(pca.multilevel.vac18, group = vac18$stimulation, # plot the samples on multilevel PCs
          ind.names = vac18$sample,
          legend = TRUE, legend.title = 'Stimulation',
          title = '(b) Multilevel PCA on VAC18 data')


## ---- echo = FALSE, eval = FALSE--------------------------------------------------------------
## Y <- vac18$stimulation
## 
## splsda.multilevel.vac18 <- splsda(X, Y, ncomp = 10)
## 
## list.keepX <- c(1:10,  seq(20, 300, 10))
## 
## # undergo the tuning process to determine the optimal number of variables
## tune.splsda.vac18 <- tune.splsda(X, Y, ncomp = 10, # calculate for first 10 components
##                                  validation = 'Mfold',
##                                  folds = 5, nrepeat = 10, # use repeated cross-validation
##                                  dist = 'max.dist', # use max.dist measure
##                                  measure = "BER", # use balanced error rate of dist measure
##                                  test.keepX = list.keepX,
##                                  cpus = 2)
## 
## splsdaTuning <- list(ncomp = tune.splsda.vac18$choice.ncomp$ncomp,
##                      keepX = tune.splsda.vac18$choice.keepX[1:tune.splsda.vac18$choice.ncomp$ncomp])
## 
## save(splsdaTuning, file = "splsdaTuning.RData")
## 


## ---- echo = FALSE----------------------------------------------------------------------------
load("splsdaTuning.RData")
optimal.ncomp <- splsdaTuning$ncomp
optimal.keepX <- splsdaTuning$keepX


## ---------------------------------------------------------------------------------------------
Y <- vac18$stimulation # use the stimulation method as the response variable vector

final.splsda.multilevel.vac18 <- splsda(X, Y, ncomp = optimal.ncomp, # undergo sPLS-DA after parameter tuning
                                        keepX = optimal.keepX,
                                        multilevel = design)


## ---------------------------------------------------------------------------------------------
optimal.ncomp # tuned number of components
optimal.keepX # tuned number of features per component


## ---- fig.cap = "FIGURE 2: Sample plot for sPLS-DA performed on the `vac18` data. Samples are projected into the space spanned by the first two components yielded by this method."----
plotIndiv(final.splsda.multilevel.vac18, group = vac18$stimulation, 
          ind.names = vac18$sample,
          legend = TRUE, legend.title = 'Treatment',
          title = 'Sample Plot of sPLS-DA on Vac18 data')


## ---- eval = FALSE----------------------------------------------------------------------------
## plotIndiv(final.splsda.multilevel.vac18,
##           group = vac18$stimulation,
##           ind.names = vac18$stimulation,
##           style = '3d')


## ---- fig.cap = "FIGURE 3: Sample plot for sPLS-DA performed on the `vac18` data. Samples are projected into the space spanned by the first three components yielded by this method.", echo = FALSE----
knitr::include_graphics("Figures/3D Sample Plot.png")


## ---- eval = FALSE----------------------------------------------------------------------------
## # set the colours used for the subject assocaited with each sample (left-most column)
## col.ID <- c("lightgreen", "red", "lightblue", "darkorange",
##               "purple", "maroon", "blue", "chocolate", "turquoise",
##               "tomato1", "pink2", "aquamarine")[vac18$sample]
## 
## cim(final.splsda.multilevel.vac18,
##     row.sideColors = cbind(color.mixo(c(vac18$stimulation)), col.ID),
##     row.names = paste(vac18$stimulation, vac18$sample, sep = "_"),
##     col.names = FALSE, legend=list(legend = c(levels(vac18$stimulation)),
##     col = c(color.mixo(1:4)),
##     title = "Stimulation", cex = 0.8))


## ---- fig.cap = "FIGURE 4: Clustered Image Map from the sPLS-DA performed on the `vac18` data. The plot displays the genetic expression levels of each measured feature for each sample. Hierarchical clustering was done with a complete Euclidean distance method. The left-most column of colours denotes which subject a given sample belongs to. The colour column to the right of this shows the treatment associated with each sample.", echo = FALSE----
knitr::include_graphics("Figures/CIM.png")


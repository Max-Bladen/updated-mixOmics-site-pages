## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', out.width = '90%')


## -------------------------------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library

set.seed(5249) # for reproducibility, remove for normal use


## -------------------------------------------------------------------------------------------------------------------
data(liver.toxicity) # extract the liver toxicity data
X <- liver.toxicity$gene # use the gene expression data as the X matrix
Y <- liver.toxicity$clinic # use the clinical data as the Y matrix

dim(X) # check the dimensions of the X dataframe
dim(Y) # check the dimensions of the Y dataframe


## ---- out.width = "49%", fig.show = "hold", fig.cap = "FIGURE 1: Barplot of the variance each principal component explains of the liver toxicity data."----
pca.gene <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
pca.clinic <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)

plot(pca.gene)
plot(pca.clinic)


## ---- out.width = "49%", fig.show = "hold", fig.cap = "FIGURE 2: Preliminary (unsupervised) analysis with PCA on the liver toxicity data"----
plotIndiv(pca.gene, comp = c(1, 2), 
          group = liver.toxicity$treatment[, 4], 
          ind.names = liver.toxicity$treatment[, 3], 
          legend = TRUE, title = 'Liver gene, PCA comp 1 - 2')

plotIndiv(pca.clinic, comp = c(1, 2), 
          group = liver.toxicity$treatment[, 4], 
          ind.names = liver.toxicity$treatment[, 3], 
          legend = TRUE, title = 'Liver clinic, PCA comp 1 - 2')


## -------------------------------------------------------------------------------------------------------------------
spls.liver <- spls(X = X, Y = Y, ncomp = 5, mode = 'regression') # form basic model


## ---- fig.cap = "FIGURE 3: Tuning the number of components in PLS on the liver toxicity data. For each component, the repeated cross-validation (5 × 10−fold CV) $Q^2$ score is shown. Horizontal line depicts $Q^2$ = 0.0975. The bars represent the variation of these values across the repeated folds."----
perf.spls.liver <- perf(spls.liver, validation = 'Mfold',
                         folds = 10, nrepeat = 5) # repeated CV tuning of component count
 
plot(perf.spls.liver, criterion = 'Q2.total')


## ---- fig.cap = "FIGURE 4: Tuning plot for sPLS2."------------------------------------------------------------------
list.keepX <- c(seq(20, 50, 5)) # set range of test values for number of variables to use from X dataframe
list.keepY <- c(3:10) # set range of test values for number of variables to use from Y dataframe


tune.spls.liver <- tune.spls(X, Y, ncomp = 2,
                              test.keepX = list.keepX,
                              test.keepY = list.keepY,
                              nrepeat = 1, folds = 10, # use 10 folds
                              mode = 'regression', measure = 'cor') # use the correlation measure for tuning
plot(tune.spls.liver)


## -------------------------------------------------------------------------------------------------------------------
tune.spls.liver$choice.keepX
tune.spls.liver$choice.keepY


## -------------------------------------------------------------------------------------------------------------------
optimal.keepX <- tune.spls.liver$choice.keepX # extract optimal number of variables for X dataframe
optimal.keepY <- tune.spls.liver$choice.keepY # extract optimal number of variables for Y dataframe
optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components



## -------------------------------------------------------------------------------------------------------------------
final.spls.liver <- spls(X, Y, ncomp = optimal.ncomp, # use all tuned values from above
                    keepX = optimal.keepX,
                    keepY = optimal.keepY,
                    mode = "regression") # explanitory approach being used, hence use regression mode


## ---- out.width = "49%", fig.show = "hold", fig.cap = "FIGURE 5: Sample plot for sPLS2 performed on the liver.toxicity data. Samples are projected into the space spanned by the components associated to each data set (or block)."----
plotIndiv(final.spls.liver, ind.names = FALSE, 
         rep.space = "X-variate", # plot in X-variate subspace
         group = liver.toxicity$treatment$Time.Group, # colour by time group
         pch = as.factor(liver.toxicity$treatment$Dose.Group), # select symbol by dose group
         col.per.group = color.mixo(1:4), 
         legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')

plotIndiv(final.spls.liver, ind.names = FALSE,
         rep.space = "Y-variate", # plot in Y-variate subspace
         group = liver.toxicity$treatment$Time.Group, # colour by time group
         pch = as.factor(liver.toxicity$treatment$Dose.Group), # select symbol by dose group
         col.per.group = color.mixo(1:4), 
         legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')


## ---- fig.cap = "FIGURE 6: Sample plot for sPLS2 performed on the liver.toxicity data. Samples are projected into the space spanned by the averaged components of both datasets."----
plotIndiv(final.spls.liver, ind.names = FALSE, 
         rep.space = "XY-variate", # plot in averaged subspace
         group = liver.toxicity$treatment$Time.Group, # colour by time group
         pch = as.factor(liver.toxicity$treatment$Dose.Group), # select symbol by dose group
         col.per.group = color.mixo(1:4), 
         legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')


## ---- eval = FALSE--------------------------------------------------------------------------------------------------
## col.tox <- color.mixo(as.numeric(as.factor(liver.toxicity$treatment[, 4]))) # create set of colours
## plotIndiv(final.spls.liver, ind.names = FALSE,
##           rep.space = "XY-variate", # plot in averaged subspace
##           axes.box = "both", col = col.tox, style = '3d')


## ---- fig.cap = "FIGURE 8:  Arrow plot from the sPLS2 performed on the liver.toxicity data. The start of the arrow indicates the location of a given sample in the space spanned by the components associated to the gene data set, and the tip of the arrow the location of that same sample in the space spanned by the components associated to the clinical data set."----
plotArrow(final.spls.liver, ind.names = FALSE,
          group = liver.toxicity$treatment$Time.Group, # colour by time group
          col.per.group = color.mixo(1:4),
          legend.title = 'Time.Group')


## ---- fig.cap = "FIGURE 9:  Stability of variable selection from the sPLS on the Liver Toxicity gene expression data. The barplot represents the frequency of selection across repeated CV folds for each selected gene for component 1 (a) and 2 (b)."----
# form new perf() object which utilises the final model
perf.spls.liver <- perf(final.spls.liver, 
                          folds = 5, nrepeat = 10, # use repeated cross-validation
                          validation = "Mfold", dist = "max.dist",  # use max.dist measure
                          progressBar = FALSE)

# plot the stability of each feature for the first two components, 'h' type refers to histogram
par(mfrow=c(1,2)) 
plot(perf.spls.liver$features$stability.X[[1]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(a) Comp 1', las =2,
     xlim = c(0, 150))
plot(perf.spls.liver$features$stability.X$comp2, type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(b) Comp 2', las =2,
     xlim = c(0, 300))


## ---- fig.cap = "FIGURE 10:  Correlation circle plot from the sPLS2 performed on the liver.toxicity data. This plot should be interpreted in relation to Figure 5 to better understand how the expression levels of these molecules may characterise specific sample groups."----
plotVar(final.spls.liver, cex = c(3,4), var.names = c(FALSE, TRUE))


## ---- eval = FALSE--------------------------------------------------------------------------------------------------
## color.edge <- color.GreenRed(50)  # set the colours of the connecting lines
## 
## # X11() # To open a new window for Rstudio
## network(final.spls.liver, comp = 1:2,
##         cutoff = 0.7, # only show connections with a correlation above 0.7
##         shape.node = c("rectangle", "circle"),
##         color.node = c("cyan", "pink"),
##         color.edge = color.edge,
##         save = 'png', # save as a png to the current working directory
##         name.save = 'sPLS Liver Toxicity Case Study Network Plot')


## ---- eval = FALSE--------------------------------------------------------------------------------------------------
## cim(final.spls.liver, comp = 1:2, xlab = "clinic", ylab = "genes")


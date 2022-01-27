## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', out.width = '90%')


## -------------------------------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library

set.seed(5249) # for reproducibility, remove for normal use


## -------------------------------------------------------------------------------------------------------------------
data(nutrimouse)
X <- nutrimouse$lipid # extract all lipid concentration variables
Y <- nutrimouse$gene # extract all gene expression variables

dim(X) # check the dimensions of the X dataframe
dim(Y) # check the dimensions of the Y dataframe


## ---- eval = FALSE--------------------------------------------------------------------------------------------------
## imgCor(X, Y, sideColors = c("purple", "green")) # produce a heat map of the cross correlation matrix


## ---- out.width='80%', fig.cap= "FIGURE 2: Heatmap of lambda1 and lambda2 values coloured by the resulting cross-validation score"----
grid1 <- seq(0.001, 0.2, length = 10) # set grid search values for each regularisation parameter
grid2 <- seq(0.001, 0.2, length = 10)

cv.tune.rcc.nutrimouse <- tune.rcc(X, Y, grid1 = grid1, grid2 = grid2, validation = "loo") # optimise the regularisation parameter values


## -------------------------------------------------------------------------------------------------------------------
cv.tune.rcc.nutrimouse # examine the results of CV tuning

opt.l1 <- cv.tune.rcc.nutrimouse$opt.lambda1 # extract the optimal lambda values
opt.l2 <- cv.tune.rcc.nutrimouse$opt.lambda2

CV.rcc.nutrimouse <- rcc(X, Y, method = "ridge", lambda1 = opt.l1, lambda2 = opt.l2) # formed optimised CV rCCA


## -------------------------------------------------------------------------------------------------------------------
shrink.rcc.nutrimouse <- rcc(X,Y, method = 'shrinkage') # run the rCCA method using shrinkage
shrink.rcc.nutrimouse$lambda # examine the optimal lambda values after shrinkage 


## ---- fig.cap= "FIGURE 3: Barplots showing canonical correlation values for each novel dimension. Left diagram depicts values using the cross-valiation method, right diagram depicts values using shrinkage method", fig.show='hold', out.width = '49%'----
plot(CV.rcc.nutrimouse, type = "barplot", main = "Cross Validation") # barplot of cross validation method rCCA canonical correlations
plot(shrink.rcc.nutrimouse, type = "barplot", main = "Shrinkage") # barplot of shrinkage method rCCA canonical correlations


## ---- fig.show = "hold", out.width = "49%", fig.cap = "FIGURE 4: rCCA sample plots from the CV (a) or shrinkage (b) method. Canonical variates corresponding to each data set are first averaged using the argument `rep.space = 'XY-variate'`. Samples are projected into the space spanned by the averaged canonical variates and coloured according to genotype information."----
plotIndiv(CV.rcc.nutrimouse, comp = 1:2, # plot the projection of samples for CV rCCA data
          ind.names = nutrimouse$genotype,
          group = nutrimouse$diet, rep.space = "XY-variate", # used averaged variate subspace
          legend = TRUE, title = '(a) Nutrimouse, rCCA CV XY-space')

plotIndiv(shrink.rcc.nutrimouse, comp = 1:2, # plot the projection of samples for shrinkage rCCA data
          ind.names = nutrimouse$genotype,
          group = nutrimouse$diet, rep.space = "XY-variate", # used averaged variate subspace
          legend = TRUE, title = '(b) Nutrimouse, rCCA shrinkage XY-space')


## ---- fig.show = "hold", out.width = "49%", fig.cap = "FIGURE 5: Arrow sample plots from the rCCA performed on the nutrimouse data to represent the samples projected onto the first two canonical variates. Each arrow represents one sample."----
plotArrow(CV.rcc.nutrimouse, group = nutrimouse$diet, # plot the arrow plot of samples for CV rCCA data
          col.per.group = color.mixo(1:5),
          title = '(a) Nutrimouse, CV method')

plotArrow(shrink.rcc.nutrimouse, group = nutrimouse$diet, # plot the arrow plot of samples for shrinkage rCCA data
          col.per.group = color.mixo(1:5),
          title = '(b) Nutrimouse, shrinkage method')


## ---- fig.show = "hold", out.width = "49%", fig.cap = "FIGURE 6:  Correlation circle plots from the rCCA performed on the nutrimouse data showing the correlation between the gene expression and lipid concentration data."----
plotVar(CV.rcc.nutrimouse, var.names = c(TRUE, TRUE),
        cex = c(4, 4), cutoff = 0.5,
        title = '(a) Nutrimouse, rCCA CV comp 1 - 2')

plotVar(shrink.rcc.nutrimouse, var.names = c(TRUE, TRUE),
        cex = c(4, 4), cutoff = 0.5,
        title = '(b) Nutrimouse, rCCA shrinkage comp 1 - 2')


## ---- eval = FALSE--------------------------------------------------------------------------------------------------
## network(CV.rcc.nutrimouse, comp = 1:2, interactive = FALSE,
##         lwd.edge = 2,
##         cutoff = 0.5)


## ---- eval = FALSE--------------------------------------------------------------------------------------------------
## cim(CV.rcc.nutrimouse, comp = 1:2, xlab = "genes", ylab = "lipids")


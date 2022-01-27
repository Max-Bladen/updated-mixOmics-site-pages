## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', out.width = '90%')


## -------------------------------------------------------------------------------------------------------------------
library(mixOmics)


## -------------------------------------------------------------------------------------------------------------------
data(multidrug) # call multidrug dataset
X <- multidrug$ABC.trans # extract ABC transporter data
dim(X) # confirm the dimension of data


## ---- fig.cap= "FIGURE 1: Explained variance of Principal Components on the Multidrug ABC Transporter data"---------
trans.spca <- spca(X, ncomp = 10, center = TRUE, scale = TRUE) # run preliminary model
plot(trans.spca) # plot the explained variance per component


## ---- eval = FALSE--------------------------------------------------------------------------------------------------
## explainedVariance <- tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE) # extract the percentages of explained variance
## plot(explainedVariance) # plot these values


## ---- fig.cap= "FIGURE 2: Tuning the number of variables to select with sPLCA on the ABC Transporter data."---------
set.seed(8589) # for reproducibility with this case study, remove otherwise
test.keepX <- c(seq(5, 30, 5)) # set the number of variable values that will be tested

tune.spca.res <- tune.spca(X, ncomp = 3, # generate the first three components
                           nrepeat = 5, # repeat the cross-validation process five times
                           folds = 3, # use three folds for the cross-validation
                           test.keepX = test.keepX)
plot(tune.spca.res) # plot the optimisation output


## -------------------------------------------------------------------------------------------------------------------
tune.spca.res$choice.keepX # how many variables per component is optimal


## -------------------------------------------------------------------------------------------------------------------
final.spca <- spca(X, ncomp = 3, # based off figure 1, three components is best
                   keepX = tune.spca.res$choice.keepX) # based off figure 2, these variables are best


## ---- fig.cap = "FIGURE 3: Sample plot from the sPCA performed on the ABC Transporter data. Samples are coloured by cell line type and numbers indicate the sample IDs"----
plotIndiv(final.spca, comp = c(1, 2), ind.names = TRUE, # plot final sPCA samples for first two components
          group = multidrug$cell.line$Class,  # use the class to colour each sample
          legend = TRUE, title = 'Multidrug transporter, sPCA comp 1 - 2')


## ---- fig.cap = "FIGURE 4: Correlation circle plot from the sPCA performed on the ABC Transporter data. Only the transporters selected by the sPCA are shown on this plot."----
plotVar(final.spca, comp = c(1, 2), var.names = TRUE,  # plot variables against the sPCA components
        title = 'Multidrug transporter, sPCA comp 1 - 2')


## ---- fig.cap = "FIGURE 5: Biplot from the sPCA performed on the ABS.trans data after variable selection. The plot highlights which transporter expression levels may be related to specific cell lines, such as melanoma."----
biplot(final.spca, cex = 0.7, # plot samples and variables against the sPCA components
       xlabs = paste(multidrug$cell.line$Class, 1:nrow(X)), # make sample names easier to discern
       group = multidrug$cell.line$Class,  # colour by sample class
       title = 'Multidrug transporter, sPCA comp 1 - 2')


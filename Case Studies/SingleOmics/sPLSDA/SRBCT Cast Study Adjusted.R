library(mixOmics) # import the mixOmics library

## -------------------------------------------------------------------------------------------------------------------

data(srbct) # extract the small round bull cell tumour data
X <- srbct$gene # use the gene expression data as the X matrix
Y <- srbct$class # use the class data as the Y matrix

dim(X) # check the dimensions of the X dataframe
Y <- unique(Y)[sample(1:4, 63, replace = TRUE)]
summary(Y) # check the distribution of class labels

## -------------------------------------------------------------------------------------------------------------------

pca.srbct = pca(X, ncomp = 10, center = TRUE, scale = TRUE) # run pca method on data
plot(pca.srbct)  # barplot of the eigenvalues (explained variance per component)

plotIndiv(pca.srbct, group = Y, ind.names = FALSE, # plot the samples projected
          legend = TRUE, title = 'PCA on SRBCT, comp 1 - 2') # onto the PCA subspace

## -------------------------------------------------------------------------------------------------------------------

srbct.splsda <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance assessment later

## -------------------------------------------------------------------------------------------------------------------
# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda , comp = 1:2, 
          group = srbct$class, ind.names = FALSE,  # colour points by class
          ellipse = TRUE, # include 95% confidence ellipse for each class
          legend = TRUE, title = '(a) PLSDA with confidence ellipses')

# use the max.dist measure to form decision boundaries between classes based on PLS-DA data
background = background.predict(srbct.splsda, comp.predicted=2, dist = "max.dist")

# plot the samples projected onto the first two components of the PLS-DA subspace
plotIndiv(srbct.splsda, comp = 1:2,
          group = srbct$class, ind.names = FALSE, # colour points by class
          background = background, # include prediction background for each class
          legend = TRUE, title = " (b) PLSDA with prediction background")

## -------------------------------------------------------------------------------------------------------------------


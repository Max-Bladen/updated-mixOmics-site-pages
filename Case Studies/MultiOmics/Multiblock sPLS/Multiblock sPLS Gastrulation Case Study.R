## ----global_options, include=FALSE------------------------------------------------------------------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, 
                      echo= TRUE, 
                      warning=FALSE, 
                      message=FALSE, 
                      fig.show=TRUE, 
                      fig.keep = 'all', 
                      fig.align = "center",
                      out.width = "70%")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library

set.seed(123) # for reproducibility, remove for normal use


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
load('Data/nmt_data_processed.RData') # load the gastrulation data

X1 <- data$rna # select three of the five dataframes to explore
X2 <- data$met_genebody
X3 <- data$acc_genebody
X <- list(rna = X1, methylation = X2, accessibility = X3) # compile these into a single X object

lapply(X, dim) # check dimensions


## ---- fig.show = "hold", out.width = "33%", fig.height = 6, fig.cap = "FIGURE 1: Circle Correlation Plots for pairwise PLS models on the gastrulation data. Only displays the top 25 features for each dimension, subsetting by those with a correlation above 0.5."----
# select arbitrary values of features to keep
list.keepX = c(25, 25)
list.keepY = c(25, 25)

# generate three pairwise PLS models
pls1 <- spls(X[["rna"]], X[["methylation"]], keepX = list.keepX, keepY = list.keepY)
pls2 <- spls(X[["rna"]], X[["accessibility"]], keepX = list.keepX, keepY = list.keepY)
pls3 <- spls(X[["methylation"]], X[["accessibility"]], keepX = list.keepX, keepY = list.keepY)

# plot features of first PLS
plotVar(pls1, cutoff = 0.5, title = "(a) RNA vs Methylation",
        legend = c("RNA", "Methylation"),
        var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(2,2),
        col = c('darkorchid', 'lightgreen'))

# plot features of second PLS
plotVar(pls2, cutoff = 0.5, title = "(b) RNA vs Accessibility",
        legend = c("RNA", "Accessibility"),
        var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(2,2),
        col = c('darkorchid', 'lightgreen'))

# plot features of third PLS
plotVar(pls3, cutoff = 0.5, title = "(c) Methylation vs Accessibility",
        legend = c("Methylation", "Accessibility"),
        var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(2,2),
        col = c('darkorchid', 'lightgreen'))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cor(pls1$variates$X, pls1$variates$Y) # calculate correlation of RNA and methylation
cor(pls2$variates$X, pls2$variates$Y) # calculate correlation of RNA and accessibility
cor(pls3$variates$X, pls3$variates$Y) # calculate correlation of methylation and accessibility


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
design = matrix(0.5, ncol = length(X), nrow = length(X), # for square matrix filled with 0.5s
                dimnames = list(names(X), names(X)))
diag(design) = 0 # set diagonal to 0s

basic.mbspls.model = block.spls(X, indY = 1, # generate basic model
                                ncomp = 5, 
                                design = design)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
choice.ncomp <- 3 
choice.keepX <- list(rna = rep(50, 3), # 50 features per each component per dataframe
                     methylation = rep(50, 3), 
                     accessibility = rep(50, 3))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
final.mbspls.model = block.spls(X, indY = 1,  # generate final model using "tuned" parameters
                                ncomp = choice.ncomp, 
                                keepX = choice.keepX,
                                design = design)


## ---- out.width = "90%", fig.cap = "FIGURE 2: Sample plot for sPLS2 performed on the gastrulation data. Samples are projected into the space spanned by the components yielded from the RNA dataset."----
plotIndiv(final.mbspls.model, ind.names = FALSE,
          group = as.factor(cell_metadata$lineage), 
          pch = as.factor(cell_metadata$stage),
          col.per.group = color.mixo(1:8), 
          legend = TRUE, legend.title = 'Lineage', legend.title.pch = 'Stage',
          blocks = 1)


## ---- echo = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(422)

samples <- sample(1:dim(data$rna)[1], 30)

X1 <- data$rna[samples, ] 
X2 <- data$met_genebody[samples, ] 
X3 <- data$acc_genebody[samples, ] 
X.arrow <- list(rna = X1, methylation = X2, accessibility = X3)

final.mbspls.model.arrow = block.spls(X.arrow, indY = 1, 
                                ncomp = choice.ncomp, 
                                keepX = choice.keepX,
                                design = design)


## ---- out.width = "90%", fig.cap = "FIGURE 3: Arrow plot from the sPLS2 performed on the gastrulation data. The star indicates the location of the centroid of a sample across all the three datsets. The tip of each arrow shows the location of that same sample in the space spanned by the components associated to a specific dataset."----
symbols <- list(rna = 1, methylation = 6, accessibility = 7)

plotArrow(final.mbspls.model.arrow, ind.names = FALSE,
          group = as.factor(cell_metadata$lineage[samples]),
          pch = symbols, pch.size = 3)


## ---- out.width = "90%", fig.cap = "FIGURE 4: Correlation circle plot from the sPLS2 performed on the gastrulation data"----------------------------------------------------
plotVar(final.mbspls.model, var.names = FALSE,
        legend = TRUE, cutoff = 0.5,
        pch = c(0,1,2))


## ---- echo = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------
X.2 = list(methylation = data$met_genebody, 
         accessibility = data$acc_genebody)

list.keepX = list(methylation = rep(15, 2), accessibility = rep(15,2))
list.keepY = c(rep(15, 2))

final.mbspls.model.circos = block.spls(X = X.2, Y = data$rna,
                                  ncomp = 2, keepX = list.keepX,
                                  keepY = list.keepY, design = 'full')


## ---- out.width = "90%", fig.cap = "FIGURE 5: Circos plot from multiblock sPLS performed on the gastrulation data The plot represents the correlations greater than 0.8 between variables of different types, represented on the side quadrants", results = FALSE----
circosPlot(final.mbspls.model.circos, group = cell_metadata$lineage, cutoff = 0.8,
           Y.name = 'rna')


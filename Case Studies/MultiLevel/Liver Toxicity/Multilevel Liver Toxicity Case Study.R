## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE,
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


## -------------------------------------------------------------------------------------------------------------------
# generate fake subject feature to mimic a repeated measures design
repeat.indiv <- c(1, 2, 1, 2, 1, 2, 1, 2, 3, 3, 4, 3, 4, 3, 4, 4, 5, 6, 5, 5,
                  6, 5, 6, 7, 7, 8, 6, 7, 8, 7, 8, 8, 9, 10, 9, 10, 11, 9, 9,
                  10, 11, 12, 12, 10, 11, 12, 11, 12, 13, 14, 13, 14, 13, 14,
                  13, 14, 15, 16, 15, 16, 15, 16, 15, 16)

design <- data.frame(sample = repeat.indiv) # load this into a dataframe

summary(as.factor(repeat.indiv)) # 16 rats, 4 measurements each


## ---- echo = FALSE, eval = FALSE------------------------------------------------------------------------------------
## list.keepX <- c(seq(20, 50, 5))
## list.keepY <- c(3:10)
## 
## # undergo the tuning process to determine the optimal number of variables
## tune.spls.liver <- tune.splslevel(X, Y, ncomp = 2, # calculate for first 10 components
##                                   multilevel = design,
##                                   test.keepX = list.keepX,
##                                   test.keepY = list.keepY,
##                                   mode = "canonical",
##                                   already.tested.X = 20,
##                                   already.tested.Y = 5) # use balanced error rate of dist measure
## 
## #dist = 'max.dist', # use max.dist measure
## # measure = "BER"
## 
## splsTuning <- list(ncomp = 2,
##                    keepX = c(20, 50),
##                    keepY = c(5, 10))
## 
## save(splsTuning, file = "splsTuning.RData")


## ---- echo = FALSE--------------------------------------------------------------------------------------------------
load("splsTuning.RData")
optimal.ncomp <- splsTuning$ncomp
optimal.keepX <- splsTuning$keepX
optimal.keepY <- splsTuning$keepY


## -------------------------------------------------------------------------------------------------------------------
optimal.ncomp # optimal number of components
optimal.keepX # optimal number of features to use for each component for the X dataframe
optimal.keepY # optimal number of features to use for each component for the Y dataframe

spls.liver.multilevel <- spls(X, Y, # generate a tuned sPLS model
                              multilevel = design,
                              ncomp = optimal.ncomp,
                              keepX = optimal.keepX, 
                              keepY = optimal.keepY,
                              mode = 'canonical')


## ---- fig.cap = "FIGURE 1: Sample plot for sPLS2 performed on the liver.toxicity data after controlling for repeated measures. Samples are projected into the space spanned by the averaged components of both datasets."----

plotIndiv(spls.liver.multilevel, rep.space = "XY", # project onto averaged components
          group = liver.toxicity$treatment$Time.Group, # use time as colour
          pch = as.factor(liver.toxicity$treatment$Dose.Group), # use dosage as symbol
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')


## ---- fig.cap = "FIGURE 2: Correlation circle plot from the sPLS2 performed on the liver.toxicity data after controlling for repeated measures."----
plotVar(spls.liver.multilevel, 
        var.names = TRUE, 
        cex = c(2,2))


## ---- eval = FALSE--------------------------------------------------------------------------------------------------
## cim(spls.liver.multilevel)


## ---- fig.cap = "FIGURE 3: Clustered Image Map from the sPLS2 performed on the liver.toxicity data after controlling for repeated measures. The plot displays the similarity values between the **X** and **Y** variables selected across two dimensions, and clustered with a complete Euclidean distance method.", echo = FALSE----
knitr::include_graphics("Figures/CIM.png")


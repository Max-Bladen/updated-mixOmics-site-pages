library(mixOmics)
coul <- color.mixo(1:3)
plot(0, type="n", xlim=c(0,100), ylim=c(-80,100), axes=FALSE,
     xlab="",ylab="", main="mixOmics methods overview", cex.main = 1.2)
box()

# PCA
rect(xleft = 20, ybottom = 75, xright = 40, ytop = 95, col=coul[1])
text(5, 85, "PCA", col = coul[1])

# PLS
rect(xleft = 20, ybottom = 50, xright = 40, ytop = 70, col=coul[1])
rect(xleft = 43, ybottom = 50, xright = 63, ytop = 70, col=coul[1])
text(5, 60, "PLS", col = coul[1])

# CCA
rect(xleft = 20, ybottom = 25, xright = 40, ytop = 45, col=coul[1])
rect(xleft = 43, ybottom = 25, xright = 63, ytop = 45, col=coul[1])
text(5, 35, "CCA", col = coul[1])

# PLS-DA
rect(xleft = 20, ybottom = 0, xright = 40, ytop = 20, col=coul[1])
rect(xleft = 43, ybottom = 0, xright = 45, ytop = 20, col=coul[2])
text(5, 10, "PLS-DA", col = coul[2])

# Multiblock PLS
rect(xleft = 20, ybottom = -25, xright = 40, ytop = -5, col=coul[1])
rect(xleft = 43, ybottom = -25, xright = 63, ytop = -5, col=coul[1])
points(x=64, y=-15, pch=16, col=coul[3], cex=0.5)
points(x=65.5, y=-15, pch=16, col=coul[3], cex=0.5)
points(x=67, y=-15, pch=16, col=coul[3], cex=0.5)
rect(xleft = 68, ybottom = -25, xright = 88, ytop = -5, col=coul[1])
text(5, -15, "Multiblock PLS", col = coul[1])

# DIABLO
rect(xleft = 20, ybottom = -50, xright = 40, ytop = -30, col=coul[1])
rect(xleft = 43, ybottom = -50, xright = 63, ytop = -30, col=coul[1])
points(x=64, y=-40, pch=16, col=coul[3], cex=0.5)
points(x=65.5, y=-40, pch=16, col=coul[3], cex=0.5)
points(x=67, y=-40, pch=16, col=coul[3], cex=0.5)
rect(xleft = 68, ybottom = -50, xright = 88, ytop = -30, col=coul[1])
rect(xleft = 91, ybottom = -50, xright = 93, ytop = -30, col=coul[2])
text(5, -40, "DIABLO", col = coul[2])
text(5, -48, "*same samples", col = coul[2], cex=0.75)

# MINT
rect(xleft = 20, ybottom = -75, xright = 40, ytop = -55, col=coul[1])
rect(xleft = 43, ybottom = -75, xright = 63, ytop = -55, col=coul[1])
points(x=64, y=-65, pch=16, col=coul[3], cex=0.5)
points(x=65.5, y=-65, pch=16, col=coul[3], cex=0.5)
points(x=67, y=-65, pch=16, col=coul[3], cex=0.5)
rect(xleft = 68, ybottom = -75, xright = 88, ytop = -55, col=coul[1])
rect(xleft = 91, ybottom = -75, xright = 93, ytop = -55, col=coul[2])
text(5, -65, "MINT", col = coul[2])
text(5, -73, "*same variables", col = coul[2], cex=0.75)

# legend
rect(xleft = 75, ybottom = 93, xright = 77, ytop = 95, col=coul[1])
text(90, 94, "Quantitative")
rect(xleft = 75, ybottom = 86, xright = 77, ytop = 88, col=coul[2])
text(90, 87, "Qualitative")
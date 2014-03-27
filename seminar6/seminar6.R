library(limma)
library(lattice)

prDat <- read.table("GSE4051_data.tsv")
str(prDat, max.level = 0)

prDes <- readRDS("GSE4051_design.rds")
str(prDes)

prepareData <- function(myGenes) {
  miniDat <- t(prDat[myGenes, ])
  miniDat <- suppressWarnings(data.frame(gExp = as.vector(miniDat),
                                         gene = rep(colnames(miniDat), each = nrow(miniDat))))
  miniDat <- suppressWarnings(data.frame(prDes, miniDat))
  miniDat
}

m <- 1000
n <- 3
x <- matrix(rnorm(m * n), nrow = m)

obsVars <- apply(x, 1, var)
summary(obsVars)

wtDes <- subset(prDes, gType == "wt")
str(wtDes)

wtDat <- subset(prDat, select = prDes$gType == "wt")
str(wtDat, max.level = 0)

wtDesMat <- model.matrix(~devStage, wtDes)
str(wtDesMat)

wtFit <- lmFit(wtDat, wtDesMat)
wtEbFit <- eBayes(wtFit)

topTable(wtEbFit)
(dsHits <- topTable(wtEbFit, coef = grep("devStage", colnames(coef(wtEbFit)))))

prepareDataWT <- function(myGenes) {
  miniDat <- t(wtDat[myGenes, ])
  miniDat <- suppressWarnings(data.frame(gExp = as.vector(miniDat),
                                         gene = rep(colnames(miniDat), each = nrow(miniDat))))
  miniDat <- suppressWarnings(data.frame(wtDes, miniDat))
  miniDat
}
hitsDat <- prepareDataWT(rownames(dsHits)[c(3,6,9)])

makeStripplot <- function(myData, ...) {
  stripplot(gExp ~ devStage | gene, myData,
            group = gType, jitter.data = TRUE,
            auto.key = TRUE, type = c('p', 'a'), grid = TRUE, ...)
}

makeStripplot(hitsDat)

hitsPval <- topTable(wtEbFit, number = Inf, coef = grep("devStage", colnames(coef(wtEbFit))), p.value = 1e-05)
dim(hitsPval)

hitsPval[63, c("F", "adj.P.Val", "devStageP6")]

P2Hits <- topTable(wtEbFit, coef = "devStageP2", number = Inf, sort = "none")
P10Hits <- topTable(wtEbFit, coef = "devStageP10", number = Inf, sort = "none")

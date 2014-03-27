library(lattice)
prDat <- read.table("GSE4051_data.tsv")
str(prDat, max.level = 0)

prDes <- readRDS("GSE4051_design.rds")
str(prDes)

(luckyGenes <- c("1419655_at","1438815_at"))

prepareData <- function(probesetIds) {
  dataSubset <- data.frame(
    sidChar = rep(colnames(prDat), length(probesetIds)),
    gExp = as.vector(t(prDat[probesetIds,])),
    gene = rep(probesetIds, each = nrow(prDes))
    )
}
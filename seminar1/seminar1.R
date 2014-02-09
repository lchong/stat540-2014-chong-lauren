# Load the dataset
prDat <- read.table("GSE4051_MINI.txt", header = TRUE, row.names = 1)

# Subset the data
weeDat <- subset(prDat, poisonFang > 7.5)
nrow(weeDat)
addmargins(with(weeDat, table(devStage, gType)))

# Print rows which are below 0.1 quantile for eggBomb expression
prDat[prDat$eggBomb < quantile(prDat$eggBomb, 0.1),]
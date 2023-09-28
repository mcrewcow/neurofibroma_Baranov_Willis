head(deeann)
deeann_control <- subset(deeann, subset = orig.ident == 'control')
deeann_treated <- subset(deeann, subset = orig.ident == 'treated')

cellchat <- createCellChat(object = deeann_control, group.by = "cell type")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.2, raw.use = FALSE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_control <- cellchat

cellchat <- createCellChat(object = deeann_treated, group.by = "cell type")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.2, raw.use = FALSE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_treated <- cellchat

ht1 <- netAnalysis_signalingRole_heatmap(cellchat_control, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_control, pattern = "incoming")
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cellchat_treated, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_treated, pattern = "incoming")
ht1 + ht2

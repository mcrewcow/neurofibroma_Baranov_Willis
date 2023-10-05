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


strwidth <- function(x) {0.5}

netVisual_chord_cell(cellchat_control, signaling = 'THBS')
netVisual_chord_cell(cellchat_treated, signaling = 'THBS')

netVisual_chord_cell(cellchat_control, signaling = 'CNTN')
netVisual_chord_cell(cellchat_treated, signaling = 'CNTN')

netVisual_chord_cell(cellchat_control, signaling = 'PERIOSTIN')
netVisual_chord_cell(cellchat_treated, signaling = 'PERIOSTIN')

netVisual_chord_cell(cellchat_control, signaling = 'NT')
netVisual_chord_cell(cellchat_treated, signaling = 'NT')

netVisual_chord_cell(cellchat_control, signaling = 'CD86')
netVisual_chord_cell(cellchat_treated, signaling = 'CD86')

netVisual_chord_cell(cellchat_control, signaling = 'NRG')
netVisual_chord_cell(cellchat_treated, signaling = 'NRG')

netVisual_chord_cell(cellchat_control, signaling = 'RELN')
netVisual_chord_cell(cellchat_treated, signaling = 'RELN')

netVisual_chord_cell(cellchat_control, signaling = 'NECTIN')
netVisual_chord_cell(cellchat_treated, signaling = 'NECTIN')

netVisual_chord_cell(cellchat_control_no_endo, signaling = 'CXCL')
netVisual_chord_cell(cellchat_treated_no_endo, signaling = 'CXCL')

netVisual_chord_cell(cellchat_control_no_endo, signaling = 'CD23')
netVisual_chord_cell(cellchat_treated_no_endo, signaling = 'CD23')

netVisual_chord_cell(cellchat_control_no_endo, signaling = 'CD6')
netVisual_chord_cell(cellchat_treated_no_endo, signaling = 'CD6')

netVisual_chord_cell(cellchat_control_no_endo, signaling = 'ALCAM')
netVisual_chord_cell(cellchat_treated_no_endo, signaling = 'ALCAM')

netVisual_chord_cell(cellchat_control_no_endo, signaling = 'KIT')
netVisual_chord_cell(cellchat_treated_no_endo, signaling = 'KIT')

netVisual_chord_cell(cellchat_control_no_endo, signaling = 'LIFR')
netVisual_chord_cell(cellchat_treated_no_endo, signaling = 'LIFR')

netVisual_chord_cell(cellchat_control_no_endo, signaling = 'LIGHT')
netVisual_chord_cell(cellchat_treated_no_endo, signaling = 'LIGHT')

netVisual_chord_cell(cellchat_control_no_endo, signaling = 'CDH')
netVisual_chord_cell(cellchat_treated_no_endo, signaling = 'CDH')

netVisual_chord_cell(cellchat_control_no_endo, signaling = 'IL6')
netVisual_chord_cell(cellchat_treated_no_endo, signaling = 'IL6')

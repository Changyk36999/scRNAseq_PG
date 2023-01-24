# reticulate::py_install(packages = 'umap-learn')
rm(list=ls())

library(Seurat)
library(CellChat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(Hmisc)
library(circlize)
library(Matrix)
library(viridis)
library(patchwork)


CellChatDB <- CellChatDB.mouse # set CellChatDB <- CellChatDB.human if working on the human dataset
interaction_input <- CellChatDB$interaction
complex_input <- CellChatDB$complex
cofactor_input <- CellChatDB$cofactor
geneInfo <- CellChatDB$geneInfo

cc_ctrl = readRDS(file = "cellchat/PG_CTRL_cellChat.rds")
cc_ir = readRDS(file = "cellchat/PG_IR_cellChat.rds")


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
#
#        Preview plots
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

paths_ctrl = cc_ctrl@netP$pathways
paths_ir = cc_ir@netP$pathways

de_paths_ctrl = setdiff(paths_ctrl, paths_ir)
de_paths_ir = setdiff(paths_ir, paths_ctrl)

venn::venn(x = list(paths_ctrl, paths_ir), snames = c("CTRL", "IR"))


intpaths = c("NT", "NRG")
netAnalysis_signalingRole_scatter(cc_ctrl)+
  netAnalysis_signalingRole_scatter(cc_ir)
netAnalysis_signalingRole_scatter(cc_ctrl, signaling = c("NRG") )+
netAnalysis_signalingRole_scatter(cc_ir, signaling = c("NRG") )
netAnalysis_signalingRole_scatter(cc_ctrl, signaling = c("NT") )+
  netAnalysis_signalingRole_scatter(cc_ir, signaling = c("NT") )


netAnalysis_contribution(cc_ctrl, signaling = c("NT"))+
  netAnalysis_contribution(cc_ir, signaling = c("NT"))
  
netAnalysis_contribution(cc_ctrl, signaling = c("NRG"))+
  netAnalysis_contribution(cc_ir, signaling = c("NRG"))
  
netAnalysis_contribution(cc_ctrl, signaling = c("EGF"))+
  netAnalysis_contribution(cc_ir, signaling = c("EGF"))

netAnalysis_contribution(cc_ctrl, signaling = c("EGF"))+
  netAnalysis_contribution(cc_ir, signaling = c("EGF"))


netVisual_heatmap(cc_ctrl, color.heatmap = "Reds")



pdf("CTRL_unique_pathways.pdf", width = 5, height = 5)
netAnalysis_signalingRole_heatmap(cc_ctrl, pattern = "outgoing", height = 4, signaling =de_paths_ctrl, width = 5, font.size = 8)
dev.off()
  
pdf("IR_unique_pathways.pdf", width = 5, height = 5)
netAnalysis_signalingRole_heatmap(cc_ir, pattern = "outgoing", height = 4, signaling =de_paths_ir, width = 5, font.size = 8)
dev.off()


netAnalysis_contribution(cc_ctrl, signaling = de_paths_ctrl[1])+
netAnalysis_contribution(cc_ctrl, signaling = de_paths_ctrl[2])+
netAnalysis_contribution(cc_ctrl, signaling = de_paths_ctrl[3])+
netAnalysis_contribution(cc_ctrl, signaling = de_paths_ctrl[4])+
netAnalysis_contribution(cc_ctrl, signaling = de_paths_ctrl[5])+
netAnalysis_contribution(cc_ctrl, signaling = de_paths_ctrl[6])+
netAnalysis_contribution(cc_ctrl, signaling = de_paths_ctrl[7])+
netAnalysis_contribution(cc_ctrl, signaling = de_paths_ctrl[8])+
netAnalysis_contribution(cc_ctrl, signaling = de_paths_ctrl[9])+
  netAnalysis_contribution(cc_ctrl, signaling = de_paths_ctrl[10])+
  netAnalysis_contribution(cc_ctrl, signaling = de_paths_ctrl[11])


netAnalysis_contribution(cc_ir, signaling = de_paths_ir[1])+
  netAnalysis_contribution(cc_ir, signaling = de_paths_ir[2])+
  netAnalysis_contribution(cc_ir, signaling = de_paths_ir[3])+
  netAnalysis_contribution(cc_ir, signaling = de_paths_ir[4])+
  netAnalysis_contribution(cc_ir, signaling = de_paths_ir[5])+
  netAnalysis_contribution(cc_ir, signaling = de_paths_ir[6])+
  netAnalysis_contribution(cc_ir, signaling = de_paths_ir[7])+
  netAnalysis_contribution(cc_ir, signaling = de_paths_ir[8])+
  netAnalysis_contribution(cc_ir, signaling = de_paths_ir[9])

plotGeneExpression(cc_ir, signaling = "ncWNT", enriched.only = F)
plotGeneExpression(cc_ctrl, signaling = "ncWNT", enriched.only = F)


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
#
#        MERGE CELLCHAT OBJECTS
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
object.list = list("CTRL" = cc_ctrl,  "IR" = cc_ir)

cellchat = mergeCellChat(object.list = list(cc_ctrl, cc_ir), add.names = c("CTRL", "IR"))


# Part I: Predict general principles of cell-cell communication
# Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

# Compare the number of interactions and interaction strength among different cell populations
# Differential number of interactions or interaction strength among different cell populations


# The differential number of interactions or interaction strength in the
# cell-cell communication network between two datasets can be visualized using
# circle plot, where red (or blue ) colored edges represent increased (or
# decreased ) signaling in the second dataset compared to the first one.

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# We can also show differential number of interactions or interaction strength
# in a greater details using a heatmap. The top colored bar plot represents the
# sum of column of values displayed in the heatmap (incoming signaling). The
# right colored bar plot represents the sum of row of values (outgoing
# signaling). In the colorbar, red (or blue ) represents increased (or decreased
# ) signaling in the second dataset compared to the first one.


gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2


# The differential network analysis only works for pairwise datasets. If there
# are more datasets for comparison, we can directly show the number of
# interactions or interaction strength between any two cell populations in each
# dataset.
#
# To better control the node size and edge weights of the inferred networks
# across different datasets, we compute the maximum number of cells per cell
# group and the maximum number of interactions (or interaction weights) across
# all datasets.

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


levels(cellchat@idents$CTRL)
group.cellType <- c(rep("Epithelial", 5), "Stromal", "Endothelial", "B-cells", rep("T-cells", 5), rep("Macro-DCs", 2), "NK-cells")
group.cellType <- factor(group.cellType, levels = c("Epithelial","Stromal", "Endothelial", "B-cells", "T-cells","Macro-DCs", "NK-cells"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)


# Compare the major sources and targets in 2D space

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)



unique(cellchat@idents$CTRL)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Acinar")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T-cells CD4+CD8+")

pdf("cellchat/combined/netAnalysis_perCell_vstreatment.pdf", width = 5.5, height = 4)
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Acinar")
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T-cells CD4+CD8+")
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T-cells CD4+")
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Etv1+")
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Myoepithelial")
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "B-cells")
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "DCs")
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macrophages")
netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Intercalated duct")
#patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()




#Part II: Identify the conserved and context-specific signaling pathways

# CellChat performs joint manifold learning and classification of the inferred
# communication networks based on their functional and topological similarity.
# NB: Such analysis is applicable to more than two datasets.
#
# Functional similarity: High degree of functional similarity indicates major
# senders and receivers are similar, and it can be interpreted as the two
# signaling pathways or two ligand-receptor pairs exhibit similar and/or
# redundant roles. NB: Functional similarity analysis is not applicable to
# multiple datsets with different cell type composition.
#
# Structural similarity: A structural similarity was used to compare their
# signaling network structure, without considering the similarity of senders and
# receivers. NB: Structural similarity analysis is applicable to multiple
# datsets with the same cell type composition or the vastly different cell type
# composition.
#
# Here we can run the manifold and classification learning analysis based on the
# functional similarity because the two datasets have the the same cell type
# composition.


library(reticulate)

#cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)



library(ComplexHeatmap)
i=1
pathway.union <- intersect(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 18)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 18)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union[1:10], title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union[1:10], title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union[1:5], title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union[1:5], title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))



# Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
levels(cellchat@idents$joint)

netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1,2),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:2),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = c(8:16), targets.use = c(2),  comparison = c(1, 2), angle.x = 45, remove.isolate = F)
netVisual_bubble(cellchat, sources.use = c(8:16), targets.use = c(1),  comparison = c(1, 2), angle.x = 45, remove.isolate = F)
netVisual_bubble(cellchat, sources.use = c(7), targets.use = c(1,2),  comparison = c(1, 2), angle.x = 45, remove.isolate = F)
netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(8:16),  comparison = c(1, 2), angle.x = 45)

#Moreover, we can identify the upgulated (increased) and down-regulated
#(decreased) signaling ligand-receptor pairs in one dataset compared to the
#other dataset. This can be done by specifying max.dataset and min.dataset in
#the function netVisual_bubble. The increased signaling means these signaling
#have higher communication probability (strength) in one dataset compared to the
#other dataset.

pdf("acinar_to_tcells_cvsi.pdf", width = 8.5, height = 3.5)
gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(8:16),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in IR", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(8:16),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in IR", angle.x = 45, remove.isolate = T)
gg1 + gg2
dev.off()


pdf("etv1_to_tcells_cvsi.pdf", width = 8.5, height = 3.5)
gg1 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(8:16),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in IR", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(8:16),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in IR", angle.x = 45, remove.isolate = T)
gg1 + gg2
dev.off()


pdf("tcells_to_acinar_cvsi.pdf", width = 8.5, height = 2)
gg1 <- netVisual_bubble(cellchat, sources.use = c(8:16), targets.use = c(1),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in IR", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(8:16), targets.use = c(1),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in IR", angle.x = 45, remove.isolate = T)
gg1 + gg2
dev.off()

pdf("mecs_to_acinar_cvsi.pdf", width = 8.5, height = 6)
gg1 <- netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(1:2),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in IR", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(1:2),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in IR", angle.x = 45, remove.isolate = T)
gg1 + gg2
dev.off()

pdf("stromal_to_acinar_cvsi.pdf", width = 8.5, height = 6)
gg1 <- netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:2),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in IR", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:2),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in IR", angle.x = 45, remove.isolate = T)
gg1 + gg2
dev.off()


#The above method for identifying the upgulated and down-regulated signaling is
#perfomed by comparing the communication probability between two datasets for
#each L-R pair and each pair of cell groups. Alternative, we can identify the
#upgulated and down-regulated signaling ligand-receptor pairs based on the
#differential gene expression analysis. Specifically, we perform differential
#expression analysis between two biological conditions (i.e., NL and LS) for
#each cell group, and then obtain the upgulated and down-regulated signaling
#based on the fold change of ligands in the sender cells and receptors in the
#receiver cells. Such analysis can be done as follows.


# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "IR"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.2, thresh.fc = 0.2, thresh.p = 0.05)
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "IR", ligand.logFC = 0.2, receptor.logFC = 0.1)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "CTRL",ligand.logFC = -0.1, receptor.logFC = -0.1)

# Since the signaling genes in the net.up and net.down might be complex with
# multi-subunits, we can do further deconvolution to obtain the individual
# signaling genes.
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)



pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1), targets.use = c(5:16), comparison = c(1, 2),  angle.x = 90, remove.isolate = T, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1), targets.use = c(5:16), comparison = c(1, 2),  angle.x = 90, remove.isolate = T, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

dev.off()
computeEnrichmentScore(net.down, species = 'mouse')
computeEnrichmentScore(net.up, species = 'mouse')



netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "count", targets.use = 1, weight.scale = T, vertex.label.cex = 0.75)
netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "weight", targets.use = 1, weight.scale = T, vertex.label.cex = 0.75)

netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "count", sources.use = 1, weight.scale = T, vertex.label.cex = 0.75)
netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "weight", sources.use = 1, weight.scale = T, vertex.label.cex = 0.75)


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "count", targets.use = 2, weight.scale = T, vertex.label.cex = 0.75)
netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "weight", targets.use = 2, weight.scale = T, vertex.label.cex = 0.75)

netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "count", sources.use = 2, weight.scale = T, vertex.label.cex = 0.75)
netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "weight", sources.use = 2, weight.scale = T, vertex.label.cex = 0.75)



par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "count", targets.use = 13, weight.scale = T, vertex.label.cex = 0.75)
netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "weight", targets.use = 13, weight.scale = T, vertex.label.cex = 0.75)

netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "count", sources.use = 13, weight.scale = T, vertex.label.cex = 0.75)
netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "weight", sources.use = 13, weight.scale = T, vertex.label.cex = 0.75)


pdf("ID to all.pdf", width = 8.5, height = 6)
gg1 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:16),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in IR", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:16),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in IR", angle.x = 45, remove.isolate = T)
gg1 + gg2
dev.off()





pdf(paste0("DiffLRpairs_celltypes.pdf"), width = 10, height = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(unique(cellchat@idents$joint))) {
  netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "count", sources.use = i, weight.scale = T, vertex.label.cex = 0.75)
  netVisual_diffInteraction(object = cellchat, comparison = c(1,2), measure = "count", targets.use = i, weight.scale = T, vertex.label.cex = 0.75)
}
dev.off()



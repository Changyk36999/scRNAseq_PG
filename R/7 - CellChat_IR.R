rm(list=ls())
# devtools::install_github("sqjin/CellChat")

library(Seurat)
library(CellChat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(Hmisc)
library(circlize)
library(Matrix)
library(viridis)

options(stringsAsFactors = FALSE)

dir.create("./cellchat/plots_IR/", showWarnings = F)

CellChatDB <- CellChatDB.mouse # set CellChatDB <- CellChatDB.human if working on the human dataset
interaction_input <- CellChatDB$interaction
complex_input <- CellChatDB$complex
cofactor_input <- CellChatDB$cofactor
geneInfo <- CellChatDB$geneInfo

# write.csv(interaction_input, file = "cellchat/interaction_input_CellChatDB.csv")
# write.csv(complex_input, file = "cellchat/complex_input_CellChatDB.csv")
# write.csv(cofactor_input, file = "cellchat/cofactor_input_CellChatDB.csv")
# write.csv(geneInfo, file = "cellchat/geneInfo_input_CellChatDB.csv")


seu <- readRDS("/gstore/project/marta235_projects/scRNAseq_PG/results/Parotid - Ctrl vs IR - annotated (filtered cells).rds")
DimPlot(seu)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
#
#        Prepare the cell chat object for the treatment of interest
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

treatment = "IR"

# input files
expMat = GetAssayData(seu, assay = "RNA", slot = "data")
#rownames(expMat)[grep(pattern = "H2-BI|H2-Ea-ps", x = rownames(expMat))]
metadata = seu@meta.data
cell.use = rownames(metadata)[metadata$Treatment == treatment]

# filter group of interest (control or IR)
data.input = expMat[, cell.use]
meta = metadata[cell.use, ]
unique(meta$cellID) # check the cell labels

# Prepare a cellchat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellID")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "cellID") # set "labels" as default cell identity
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# use CellChatDB.mouse 
CellChatDB <- CellChatDB.mouse 
#showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB # use the default CellChatDB

# remove H2-BI and H2-Ea-ps from interaction database because they cause errors - developers will fix in subsequent update
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-which(CellChatDB.use[["interaction"]]$ligand == "H2-BI") ,]
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-which(CellChatDB.use[["interaction"]]$ligand == "H2-Ea-ps") ,]


cellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database



# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
#
#        Determine significant interactions
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

#if (future::supportsMulticore()) {
#  future::plan(future::multicore)
#} else {
#  future::plan(future::multisession)
#}

# nCores = parallel::detectCores()
# future::plan("multisession", workers = nCores)

cellchat <- identifyOverExpressedGenes(cellchat, thresh.p = 0.1)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, population.size = TRUE, seed.use = 3215)
# cellchat <- computeCommunProb(cellchat,  population.size = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 2)

saveRDS(cellchat, file = paste0("cellchat/PG_", treatment,"_cellChat.rds"))

# cellchat = readRDS("cellchat/PG_IR_cellChat.rds")

# df.net <- subsetCommunication(cellchat)
levels(cellchat@idents)
df_fromAcinar <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(1:16)) # from Acinar and ETV1+
df_toAcinar <- subsetCommunication(cellchat, sources.use = c(1:16), targets.use = c(1,2)) # from Acinar and ETV1+
df_selPaths <- subsetCommunication(cellchat, signaling = c("EGF", "NRG", "NT", "EPHA"))
df_all <- subsetCommunication(cellchat)

write.csv(df_fromAcinar, file = paste0("cellchat/PG_",treatment,"_cellChat_LR_fromAcinar.csv"), row.names = F)
write.csv(df_toAcinar, file = paste0("cellchat/PG_",treatment,"_cellChat_LR_toAcinar.csv"), row.names = F)
write.csv(df_selPaths, file = paste0("cellchat/PG_",treatment,"_cellChat_LR_selPaths.csv"), row.names = F)
write.csv(df_all, file = paste0("cellchat/PG_",treatment,"_cellChat_LR_allInteractions.csv"), row.names = F)


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
#
#        Pathways represented by significant interactions
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

cellchat <- computeCommunProbPathway(cellchat, thresh = 0.05)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
#par(mfrow = c(1,2), xpd=TRUE)
#netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
#netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#dev.off()
dev.off()

mat <- cellchat@net$weight
pdf(paste0("cellchat/celllchat_PG_",treatment,"_celltypes.pdf"), width = 20, height = 20)
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


## Part III: Visualization of cell-cell communication network
# Hierarchy plot: USER should define vertex.receiver, which is a numeric vector
# giving the index of the cell groups as targets in the left part of hierarchy
# plot. This hierarchical plot consist of two components: the left portion shows
# autocrine and paracrine signaling to certain cell groups of interest (i.e, the
# defined vertex.receiver), and the right portion shows autocrine and paracrine
# signaling to the remaining cell groups in the dataset. Thus, hierarchy plot
# provides an informative and intuitive way to visualize autocrine and paracrine
# signaling communications between cell groups of interest. For example, when
# studying the cell-cell communication between fibroblasts and immune cells,
# USER can define vertex.receiver as all fibroblast cell groups.
#
# Chord diagram: CellChat provides two functions netVisual_chord_cell and
# netVisual_chord_gene for visualizing cell-cell communication with different
# purposes and different levels. netVisual_chord_cell is used for visualizing
# the cell-cell communication between different cell groups (where each sector
# in the chord diagram is a cell group), and netVisual_chord_gene is used for
# visualizing the cell-cell communication mediated by mutiple ligand-receptors
# or signaling pathways (where each sector in the chord diagram is a ligand,
# receptor or signaling pathway.)
#
# Explnations of edge color/weight, node color/size/shape: In all visualization
# plots, edge colors are consistent with the sources as sender, and edge weights
# are proportional to the interaction strength. Thicker edge line indicates a
# stronger signal. In the Hierarchy plot and Circle plot, circle sizes are
# proportional to the number of cells in each cell group. In the hierarchy plot,
# solid and open circles represent source and target, respectively. In the Chord
# diagram, the inner thinner bar colors represent the targets that receive
# signal from the corresponding outer bar. The inner bar size is proportional to
# the signal strength received by the targets. Such inner bar is helpful for
# interpreting the complex chord diagram. Note that there exist some inner bars
# without any chord for some cell groups, please just igore it because this is
# an issue that has not been addressed by circlize package.

dev.off()
pathways.show <- c("NRG") 
levels(cellchat@idents)
vertex.receiver = seq(1:2) # a numeric vector corresponding to acinar and ETV1 Cells
pathways.show.all <- cellchat@netP$pathways



pdf(paste0(treatment, "_PG_AcinarEtv1_circleplot.pdf"), width = 6, height = 6)
netVisual_aggregate(cellchat, signaling = pathways.show.all,  
                    vertex.receiver = vertex.receiver, layout = "circle", sources.use = c(1,2))


netVisual_aggregate(cellchat, signaling = pathways.show.all,  
                    vertex.receiver = vertex.receiver, layout = "circle", targets.use = c(1,2))

dev.off()



# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling 
# to fibroblast and the right portion shows signaling to immune cells 
# netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
# Circle plot
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# chord plot
# netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "NRG", color.heatmap = "Reds")+
  netVisual_heatmap(cellchat, signaling = "NT", color.heatmap = "Reds")+
  netVisual_heatmap(cellchat, signaling = "EPHA", color.heatmap = "Reds")+
  netVisual_heatmap(cellchat, signaling = "CEACAM", color.heatmap = "Reds")+
  netVisual_heatmap(cellchat, signaling = "THBS", color.heatmap = "Reds")

#> Do heatmap based on a single object
netAnalysis_contribution(cellchat, signaling = pathways.show)



# ::::::::::: AUTOMATICALLY SAVE ALL SIGNIFICANT PLOTS ::::::::::::::::::::::::: #

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,2)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i],# vertex.receiver = vertex.receiver, 
            layout = "circle", out.format = "pdf", show.legend = F, vertex.label.cex = 0.5)
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i], font.size = 5, font.size.title = 7)
  ggsave(filename=paste(pathways.show.all[i], "_LR_contribution.pdf"), plot=gg, 
         width = 3, height = 2, units = 'in', dpi = 300, path = "cellchat/plots/")
}
dev.off()



netAnalysis_contribution(cellchat, signaling = c("NT"))+
  netAnalysis_contribution(cellchat, signaling = c("NRG")) #+
netAnalysis_contribution(cellchat, signaling = c("CEACAM")) +
  netAnalysis_contribution(cellchat, signaling = c("THBS")) 


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
#
#        Additional visualizations by cell type & pathway
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
levels(cellchat@idents)
#Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
# Let's generate some plots for EGF, NT, NGF, & NRG pathways

pdf(paste0(treatment, "PG_Cell_type_interactions_bubble.pdf"), width = 6,height = 5)
netVisual_bubble(cellchat, sources.use = c(5:8), targets.use = c(1,2),  remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1,2), targets.use = c(1:16), remove.isolate = FALSE)
dev.off()


pdf(paste0(treatment, "PG_Cell_type_interactions_bubble.pdf"), width = 6,height = 6)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:16), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:16), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:16), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:16), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:16), remove.isolate = FALSE)
dev.off()


pdf(paste0(treatment, "PG_Pathway_expression_violin.pdf"), width = 6,height = 5)
plotGeneExpression(cellchat, signaling = "EGF")
plotGeneExpression(cellchat, signaling = "NRG")
plotGeneExpression(cellchat, signaling = "NT")
plotGeneExpression(cellchat, signaling = "EPHA")
plotGeneExpression(cellchat, signaling = "VEGF")
plotGeneExpression(cellchat, signaling = "NOTCH")
dev.off()



# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
#
#       Identify signaling roles (e.g., dominant senders, receivers) of 
#           cell groups as well as the major contributing signaling
#
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pdf(file = paste0(treatment, "PG_SenderReceiver_pathways.pdf"), width = 8, height = 2.5)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all, width = 8, height = 2.5, font.size = 10)
dev.off()

pdf(paste0(treatment, "_PG_aggregatedHeatmap.pdf"), width = 20, height = 30)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 35, signaling = pathways.show.all, width = 12, font.size = 8)+
  netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 35, signaling = pathways.show.all, width = 12, font.size = 8)
dev.off()

pdf(paste0(treatment, "_PG_aggregatedHeatmap_selPaths.pdf"), width = 5, height = 3)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 3, signaling = c("NRG", "NT", "THBS", "EPHA", "CEACAM"), width = 5, font.size = 8)+
  netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 3, signaling = c("NRG", "NT", "THBS", "EPHA", "CEACAM"), width = 5, font.size = 8)
dev.off()

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 15, signaling = pathways.show.all[1:48])
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 15, signaling = pathways.show.all[1:48])
ht3 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 15, signaling = pathways.show.all[50:96])
ht4 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 15, signaling = pathways.show.all[50:96])

pdf(file = paste0(treatment, "PG_SenderReceiver_pathways_aggregated.pdf"), width = 12, height = 15)
ht1 + ht2
ht3 + ht4
dev.off()


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_scatter(cellchat)
netAnalysis_signalingRole_scatter(cellchat, signaling = c("NRG", "NT") )


# transfer to subdirectory to organize outputs
file.copy(from = list.files(pattern = ".pdf", path = "./", full.names = T), to = "cellchat/plots_IR//")
unlink(list.files(pattern = ".pdf", path = "./", full.names = T))


saveRDS(cellchat, file = paste0("cellchat/PG_", treatment,"_cellChat.rds")) # save final rds

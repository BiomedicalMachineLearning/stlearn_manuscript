# Predicting CCIs using CellChat for comparison with stlearn; point is to show
# that CellChat will predict interactions for cell types which never spatially
#                                                                      co-occur.
#
#       INPUT: * /Volumes/GML001-Q1851/Brad/breast_ClusterCCIResults.h5ad
#       OUTPUT: * data/breast/cluster/cellchat*

################################################################################
                        # Environment setup #
################################################################################
# TODO: NOTE must be run in folder with README entitled StLearn Reproduce
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/stLearn/"
setwd(work_dir)

library(CellChat)
library(reticulate)

library(stringr)
library(ggplot2)
# Using reticulate with anndata package to read data #
library(reticulate)
ad <- import('anndata')

data_dir <- '/Volumes/GML001-Q1851/Brad/'
data_dir2 <- 'data/'
out_dir <- 'data/breast/cluster/'

################################################################################
                          # Loading the dataset #
################################################################################
ad_ <- ad$read_h5ad(paste0(data_dir, 'breast_ClusterCCIResults.h5ad'))

# Getting the normalised expression #
expr <- t(ad_$to_df())

# Getting the cell meta information #
meta <- as.data.frame(ad_$obs)
meta[,'cell_type'] <- as.character( as.integer(meta[,'cell_type']) )

# Creating the cellchat object #
cellchat <- createCellChat(object = expr, meta = meta, group.by = "cell_type")

################################################################################
                # Loading the CellChat database #
################################################################################
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

################################################################################
                  # Preprocessing for CellChat CCI #
################################################################################
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

# First calling DE genes between groups #
cellchat <- identifyOverExpressedGenes(cellchat)

# Now identifying which LRs are 'over-expressed' #
cellchat <- identifyOverExpressedInteractions(cellchat)

################################################################################
                  # Inference of CellChat CCI #
################################################################################
# Inferring CCIs #
cellchat <- computeCommunProb(cellchat) #nboot defines permutations. 

# Filter out the cell-cell communication if there are only few number of cells 
# in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 5)

# Calculating overall interactions across LR pairs #
cellchat <- aggregateNet(cellchat)

################################################################################
              # Visualising overall CCI results #
################################################################################
groupSize <- as.numeric(table(cellchat@idents))

# Saving the overall CCI network #
# plot_name <- paste0(out_plots, data_set, '_overall_network.pdf')
# pdf(plot_name)
# par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
#                  label.edge= F, title.name = "Interaction weights/strength")
# dev.off()

# Overall interactions per cell type #
# plot_name <- paste0(out_plots, data_set, '_per-cell_overall_network.pdf')
# pdf(plot_name)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# dev.off()

##### Saving the int_df ######
int_df <- cellchat@net$count
write.table(int_df, paste0(out_dir, 'cell_chat_ints.txt'), 
            sep='\t', quote=F)

### also saving the cell chat object ###
saveRDS(cellchat, file=paste0(out_dir,'cellchat.rds'))








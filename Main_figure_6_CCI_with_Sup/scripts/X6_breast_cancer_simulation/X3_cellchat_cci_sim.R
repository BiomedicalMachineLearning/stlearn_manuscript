#  Running CellChat on the simulated data.
#
#             INPUT:  * data/sim_data/spatialsim_v2.h5ad
#             OUTPUT: * data/sim_data/methods_out/cell_chat_ints.txt

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
# Using reticulate with anndata package to read data, set to conda env #
Sys.setenv(RETICULATE_PYTHON = "~/opt/miniconda3/envs/STI/bin/python/")
library(reticulate)
ad <- import('anndata')
sc <- import("scanpy")

data_dir <- 'data/sim_data/'
out_dir <- 'data/sim_data/methods_out/'

################################################################################
                          # Loading the dataset #
################################################################################
ad_ <- ad$read_h5ad(paste0(data_dir, 'spatialsim_v2.h5ad'))

# Getting the normalised expression #
expr <- t(ad_$to_df())

# Getting the cell meta information #
decon <- ad_$obsm['deconvolution']
cts<- as.character(apply(decon, 1, 
                  function(vals) {sample(names(vals)[vals==max(vals)], 1)}))
ad_$obs['cell_type'] <- cts
# Saving with the cell type information! #
ad_$write_h5ad(paste0(data_dir, 'simdata.h5ad'), compression='gzip')

###### Creating the cellchat object #########
meta <- as.data.frame(ad_$obs)
meta[,'cell_type'] <- as.character( as.integer(meta[,'cell_type']) )

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
ct_labels <- as.data.frame(ad_$obs)[,'cell_type']
ints <- colnames(int_df)
int_labels <- meta[,'cell_type']
ct_set <- unlist(lapply(ints, function(val) {ct_labels[int_labels==val][1]}))
colnames(int_df) <- ct_set
rownames(int_df) <- ct_set
write.table(int_df, paste0(out_dir, 'cell_chat_ints.txt'), 
            sep='\t', quote=F)






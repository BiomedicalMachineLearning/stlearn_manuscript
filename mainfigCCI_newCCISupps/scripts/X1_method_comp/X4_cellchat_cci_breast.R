#(21/09/21) -> Predicting CCIs using CellChat for comparison with stlearn; point
#                is to show that CellChat will predict interactions for cell
#                types which never spatially co-occur.
#
#                INPUT: * /Volumes/GML001-Q1851/Brad/
#                               breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad
#                        * data/label_transfer_bc.csv
#                OUTPUT: * data/method_comp/breast/cellchat*
#                        * plots/method_comp/cellchat*

################################################################################
                        # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/stLearn/"
setwd(work_dir)

Sys.setenv(RETICULATE_PYTHON = "~/opt/miniconda3/envs/STI/bin/python/")
library(CellChat)
library(reticulate)

library(stringr)
library(ggplot2)
library(reticulate)
ad <- import('anndata')

data_dir <- '/Volumes/GML001-Q1851/Brad/'
data_dir2 <- 'data/'
st_db_dir <- '/Users/uqbbalde/Desktop/Uni_Studies/myPython/stlearn_latest/stLearn/stlearn/tools/microenv/cci/databases/'
out_dir <- 'data/method_comp/breast/cellchat_'
out_plots <- 'plots/method_comp/cellchat_'

data_set <- 'breast'

################################################################################
                          # Loading the dataset #
################################################################################
ad_ <- ad$read_h5ad(paste0(data_dir, 'breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad'))

# Getting the normalised expression #
expr <- t(ad_$to_df())
rownames(expr) <- str_replace_all(rownames(expr), prefix, '')

# Getting the cell meta information #
meta <- as.data.frame(ad_$obs)
label_transfer = read.table(paste0(data_dir2, 'label_transfer_bc.csv'),
                     sep='\t', row.names=1)
bcs <- rownames(label_transfer)

meta[bcs, 'cell_type'] <- label_transfer[,'predicted.id']

# Creating the cellchat object #
cellchat <- createCellChat(object = expr, meta = meta, group.by = "cell_type")

################################################################################
                # Loading the CellChat database #
################################################################################
CellChatDB <- CellChatDB.human

# Checking overlap with current stlearn NATMI database #
st_db <- read.table(paste0(st_db_dir, 'connectomeDB2020_lit.txt'), 
                    sep='\t', header=T)
st_lrs <- str_c(st_db[,1], st_db[,2], sep='_')

cc_lrs <- as.character(CellChatDB$interaction$interaction_name)

overlap <- intersect(st_lrs, cc_lrs) # 545 overlap
cc_only <- setdiff(cc_lrs, st_lrs) # 1476
st_only <- setdiff(st_lrs, cc_lrs) # 1748

######## This overall suggests that will need to add the connectomeDB to 
######## CellChat for a fair comparison of methods.
#### Instructions on how to do that are here:
# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Update-CellChatDB.html
interaction_input <- CellChatDB$interaction
gene_info <- CellChatDB$geneInfo

allowed_genes <- gene_info[,'Symbol']
allowed_genes <- allowed_genes[order(allowed_genes)]

### Need to update the above with the stlearn interactions ####
lrs <- st_only
stlearn_input <- list("interaction_name"=c(), "pathway_name"=c(), "ligand"=c(),
                      "receptor"=c(), "agonist"=c(), "antagonist"=c(),
                      "co_A_receptor"=c(), "co_I_receptor"=c(), "evidence"=c(),
                      "annotation"=c(), "interaction_name_2"=c())
wrong <- c()
for (i in 1:length(lrs)) {
  lr_ <- lrs[i]
  l_ <- str_split(lr_, '_')[[1]][1]
  r_ <- str_split(lr_, '_')[[1]][2]
  
  if (!(l_%in%allowed_genes)) {wrong <- c(wrong, l_)}
  if (!(r_%in%allowed_genes)) {wrong <- c(wrong, r_)}
  if (!(l_%in%allowed_genes) || !(r_%in%allowed_genes)) {next}
  
  stlearn_input[["interaction_name"]] <- c(stlearn_input[["interaction_name"]], lr_)
  stlearn_input[["pathway_name"]] <- c(stlearn_input[["pathway_name"]], "")
  stlearn_input[["ligand"]] <- c(stlearn_input[["ligand"]], l_)
  stlearn_input[["receptor"]] <- c(stlearn_input[["receptor"]], r_)
  for (j in 5:length(stlearn_input)) {
    stlearn_input[[j]] <- c(stlearn_input[[j]], "")
  }
}
wrong <- unique(wrong)
stlearn_input <- as.data.frame(stlearn_input)
rownames(stlearn_input) <- stlearn_input[,1]

# Add this to the database...
new_interactions <- rbind(interaction_input, stlearn_input)

# Setting this new database to the cellchat object #
CellChatDB$interaction <- new_interactions
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
plot_name <- paste0(out_plots, data_set, '_overall_network.pdf')
pdf(plot_name)
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
#                  label.edge= F, title.name = "Interaction weights/strength")
dev.off()

# Overall interactions per cell type #
plot_name <- paste0(out_plots, data_set, '_per-cell_overall_network.pdf')
pdf(plot_name)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                   edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

################################################################################
# Visualising per LR CCI results #
################################################################################
df.net <- subsetCommunication(cellchat) # Retrieves the results as a dataframe.

############### Checking if LRs of interest are present ########################
int_lrs <- c('MCAM_MCAM', 'NECTIN4_NECTIN4', 
             'GPC3_IGF1R', 'SPP1_ITGB1', 'TGM2_ADGRG1')

sig_lrs <- as.character(df.net$interaction_name[df.net$pval<.05])

sig_int_lrs <- intersect(int_lrs, sig_lrs)
print(sig_int_lrs)

print("IL34_CSF1R" %in% cc_lrs)
print("IL34_CSF1R" %in% sig_int_lrs)
# NOTE: was not significant & was present in the CC LR database #

############### Visualising the networks for LRs of interest ########################
for (i in 1:length(sig_int_lrs)){
  lr_ <- sig_int_lrs[i]
  
  ##### Circle plot #######
  plot_name <- paste0(out_plots, data_set, '_', lr_, '_network.pdf')
  pdf(plot_name)
  netVisual_individual(cellchat, signaling = "", pairLR.use = lr_, layout = "circle")
  dev.off()
  
  ##### Chord diagrams #######
  # Couldn't get working:
  # https://github.com/sqjin/CellChat/issues/235
  
  plot_name <- paste0(out_plots, data_set, '_', lr_, '_chord.pdf')
  par(xpd=T)
  pdf(plot_name)
  netVisual_individual(cellchat, signaling = "", pairLR.use = lr_, layout = "chord",
                       )
  dev.off()
}

#### Overall view equivalent to what is generated for stLearn ####
pairLR.use <- as.data.frame(int_lrs)
colnames(pairLR.use) <- c('interaction_name')

# This function dosn't seem to work.... #
netVisual_bubble(cellchat, pairLR.use = pairLR.use, remove.isolate = T, 
                 group=c(1,2,2)
)

# Saving the results #
write.table(df.net, file=paste0(out_dir, data_set, '_interaction_data.txt'),
            quote=F, sep='\t')
saveRDS(cellchat, file=paste0(out_dir, data_set, '.rds'))
#saveRDS(cellchat, file=paste0(out_dir, data_set, '_10000.rds'))








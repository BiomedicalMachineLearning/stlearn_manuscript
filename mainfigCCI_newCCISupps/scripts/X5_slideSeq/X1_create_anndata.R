# Loads in the slide-seq data from Seurat & converts to anndata, subsequently
# saving as .h5ad for later loading in python. Also performs the Seurat label
# transfer onto this data as per the vignette, using downloaded reference scRNA-seq data.
# 
#       INPUT: * data/slideSeq/mouse_hippocampus_reference.rds
#                 -> Downloaded from here:
#                 https://www.dropbox.com/s/cs6pii5my4p3ke3/mouse_hippocampus_reference.rds?dl=0
#              * Uses Seurat to download data of hippocampus.
#       OUTPUT:
#              * /Volumes/GML001-Q1851/Brad/slideSeq/hipp.h5ad

###############################################################################
                        # Environment Setup #
###############################################################################
# TODO change this to your directory
work_dir <- '/Users/uqbbalde/Desktop/Uni_Studies/projects/stlearn_manuscript/mainfigCCI_newCCISupps/'
setwd(work_dir)
data_dir <- '/Volumes/GML001-Q1851/Brad/slideSeq/'
out_dir <- data_dir

library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
Sys.setenv(RETICULATE_PYTHON = "~/opt/miniconda3/envs/STI/bin/python/")
library(reticulate)
np <- import('numpy')
pd <- import('pandas')
sc <- import('scanpy')
pil <- import('PIL')

InstallData("ssHippo")

###############################################################################
                      # Loading the data #
###############################################################################
slide.seq <- LoadData("ssHippo")

# Just visualising to make sure data loaded correctly #
plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
### Looks good!!

###############################################################################
                    # Converting to anndata #
###############################################################################
counts <- as.data.frame(slide.seq@assays$Spatial@counts)
counts <- t(counts)

# Retrieving the spatial locations #
coords <- slide.seq@images$image@coordinates[,c('x','y')]

####### Creating the anndata ########
counts_py <- r_to_py( counts )
counts_df_py <- pd$DataFrame(counts_py, 
                             index=rownames(counts), columns=colnames(counts))
data = sc$AnnData(counts_df_py)

data$obs["imagecol"] = coords['y'] 
data$obs["imagerow"] = coords['x'] 

###############################################################################
     # Loading in reference scRNA-seq data & performing label transfer #
###############################################################################
### Preprocessing ###
slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = T)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)

### Loading reference ###
ref <- readRDS(paste0(data_dir, 'mouse_hippocampus_reference.rds'))

### Running label transfer ###
anchors <- FindTransferAnchors(reference = ref, query = slide.seq, normalization.method = "SCT",
                               npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = ref$celltype, prediction.assay = TRUE,
                                  weight.reduction = slide.seq[["pca"]], dims = 1:50)
slide.seq[["predictions"]] <- predictions.assay

slide.seq$predicted.id <- GetTransferPredictions(slide.seq, score.filter=.75)
Idents(slide.seq) <- "predicted.id"
SpatialDimPlot(slide.seq, cells.highlight = CellsByIdentities(object = slide.seq, 
                                                              idents = c("CA3 Principal cells",
                                                                         "Dentate Principal cells", 
                                                                         "Endothelial tip")), 
                facet.highlight = TRUE)

SpatialDimPlot(slide.seq)

###############################################################################
    # Adding label transfered cell type labels to anndata & saving #
###############################################################################
print( all( data$obs_names$values==names(slide.seq$predicted.id) ) )
data$obs['cell_type'] <- slide.seq$predicted.id

## Saving the anndata, will require further formatting in python ##
#data$write_h5ad(paste0(out_dir, 'hipp.h5ad'), compression='gzip')
data$write_h5ad(paste0(out_dir, 'hipp_rep.h5ad'), compression='gzip')














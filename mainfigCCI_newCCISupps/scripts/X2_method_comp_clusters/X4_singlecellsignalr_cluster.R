# Runs SingleCellSignalR on the breast data using the cluster, then
#                saves the overall CCI information inferred and by what LR pairs.
#
#                INPUT: * /Volumes/GML001-Q1851/Brad/
#                                                  breast_ClusterCCIResults.h5ad
#                OUTPUT: * data/breast/cluster/singlecellsignalr*

################################################################################
                    # Environment setup #
################################################################################
# TODO: NOTE must be run in folder with README entitled StLearn Reproduce
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/stLearn/"
setwd(work_dir)
library(SingleCellSignalR)
library(stringr)
library(ggplot2)
# Using reticulate to load the anndata #
library(reticulate)
sc <- import('scanpy')

data_dir <- '/Volumes/GML001-Q1851/Brad/'
out_dir <- 'data/breast/cluster/'

################################################################################
                        # Loading the data #
################################################################################
# Actually, need to load the equivalent LRs from the breast data #
breast_ad <- sc$read_h5ad(paste0(data_dir, 'breast_ClusterCCIResults.h5ad'))
expr <- as.data.frame( t(breast_ad$to_df()) )
clusters <- as.numeric(breast_ad$obs$cell_type)
cell_types <- as.character(breast_ad$obs$cell_type)
c.names <- c()
for (i in 1:length(unique(clusters))) {
  c.names <- c(c.names, cell_types[clusters==i][1])
}

################################################################################      
                      # Inferring the CCIs #
################################################################################
# The latter produced more interactions, including between same cell types #
genes <- rownames(expr)
signal1 <- cell_signaling(expr, genes, clusters, c.names=c.names,
                         species='homo sapiens', write=F, int.type='paracrine')
signal2 <- cell_signaling(expr, genes, clusters, c.names=c.names, 
                          species='homo sapiens', write=F, int.type='autocrine')
signal <- list() # Binding all results
all_ccis <- unique(c(names(signal1), names(signal2)))
for (i in 1:length(all_ccis)) {
  cci <- all_ccis[i]
  if ((cci %in% names(signal1)) & (cci %in% names(signal2))) {
    signal[[cci]] <- rbind(signal1[[cci]], signal2[[cci]])
  } else if (cci %in% names(signal1)) {
    signal[[cci]] <- signal1[[cci]]
  } else if (cci %in% names(signal2)) {
    signal[[cci]] <- signal2[[cci]]
  }
}

################################################################################      
                          # Summarising results #
################################################################################
# Great! Saving the results. #
saveRDS(signal, file=paste0(out_dir, 'singlecellsignalr.rds'))
# signal <- readRDS( paste0(out_dir, 'singlecellsignalr.rds') )

####### Now also doing the necessary counting if LRs per CCI! #######
# CCIs & CCI counts per LR #
n_ccis <- length(all_ccis)
int_df <- as.data.frame(matrix(0, length(c.names), length(c.names)),
                        row.names=c.names,)
colnames(int_df) <- c.names
for (i in 1:n_ccis) {
  cci <- all_ccis[i]
  rowi <- as.integer( str_split(cci, '-')[[1]][1] ) + 1
  coli <- as.integer( str_split(cci, '-')[[1]][2] ) + 1
  lr_df <- signal[[i]]
  n_lrs <- dim(lr_df)[1]
  int_df[rowi, coli] <- n_lrs
}

# Creating & saving DF #
write.table(int_df, file=paste0(out_dir, 'singlecellsignalr_ints.txt'), 
            sep='\t', quote=F)











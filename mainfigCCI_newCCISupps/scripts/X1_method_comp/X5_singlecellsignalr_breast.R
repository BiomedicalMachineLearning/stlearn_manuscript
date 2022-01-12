# (21/10/21) -> Runs SingleCellSignalR on the breast
#               data using the cell type labels, then saves the overall
#               CCI information inferred and by what LR pairs. 
# 
# INPUT: * /Volumes/GML001-Q1851/Brad/breast_LR&CCIResults.h5ad
# OUTPUT: * data/method_comp/breast/singlecellsignalr/*

################################################################################
# Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/stLearn/"
setwd(work_dir)

Sys.setenv(RETICULATE_PYTHON = "~/opt/miniconda3/envs/STI/bin/python/")
library(SingleCellSignalR)
library(stringr)
library(ggplot2)
library(reticulate)
sc <- import('scanpy')

data_dir <- 'data/method_comp/breast/'
data_dir2 <- '/Volumes/GML001-Q1851/Brad/'
st_db_dir <- '/Users/uqbbalde/Desktop/Uni_Studies/myPython/stlearn_latest/stLearn/stlearn/tools/microenv/cci/databases/'
out_dir <- 'data/method_comp/breast/singlecellsignalr/'

data_set <- 'breast'

################################################################################
                        # Loading the data #
################################################################################
# Actually, need to load the equivalent LRs from the breast data #
breast_ad <- sc$read_h5ad(paste0(data_dir2, 'breast_LR&CCIResults.h5ad'))
expr <- as.data.frame( t(breast_ad$to_df()) )
clusters <- as.numeric(breast_ad$obs$cell_type)
cell_types <- as.character(breast_ad$obs$cell_type)
c.names <- c()
for (i in 1:length(unique(clusters))) {
  c.names <- c(c.names, cell_types[clusters==i][1])
}
st_lrs <- names(breast_ad$uns[['per_lr_cci_cell_type']])

################################################################################      
                    # Updating the database #
################################################################################
# No explicit instructions on how to update the database given here:
# https://bioconductor.org/packages/release/bioc/vignettes/SingleCellSignalR/inst/doc/UsersGuide.html
# But based on the code, it appears it uses the dataframe LRdb within the functions 
# to infer CCI:
# https://github.com/SCA-IRCM/SingleCellSignalR/blob/master/R/cell_signaling.R
# HENCE; will overwrite this database & continue vignette to see if you can override. 

st_genes <- as.data.frame( str_split(st_lrs, '_') )
st_db <- as.data.frame(list('ligand'=as.character(st_genes[1,]), 
                            'receptor'=as.character(st_genes[2,]),
                            'source'=rep('NATMI', length(st_lrs)), 
                            'PMIDs'=rep('', length(st_lrs)))
)

# Overriding the DB #
LRdb <- st_db

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

# Making sure all interactions in the st_db #
sig_lrs <- c()
for (i in 1:length(signal)){
  sig_lrs <- c(sig_lrs, 
               str_c(signal[[i]][,1], signal[[i]][,2], sep='_'))
}
sig_lrs <- unique(sig_lrs)
print(all(!is.na(match(sig_lrs, st_lrs))))
# Great!

################################################################################      
                          # Summarising results #
################################################################################
# Great! Saving the results. #
saveRDS(signal, file=paste0(out_dir, 'singlecellsignalr.rds'))
# signal <- readRDS( paste0(out_dir, 'singlecellsignalr.rds') )

####### Now also doing the necessary counting! #######
# CCIs & CCI counts per LR #
ccis <- character(length(st_lrs))
cci_counts <- integer(length(st_lrs))
cci_scores <- double(length(st_lrs))
for (i in 1:length(st_lrs)) {
  ccis_i <- c()
  cci_scores_i <- c()
  lr <- st_lrs[i]
  for (j in 1:length(signal)) {
    lrs <- str_c(signal[[j]][,1], signal[[j]][,2], sep='_')
    if (lr %in% lrs) {
      cci <- str_replace_all(names(signal)[j], '-', '--')
      ccis_i <- c(ccis_i, cci)
      cci_scores_i <- c(cci_scores_i, signal[[j]][which(lrs==lr),'LRscore'])
    }
  }
  ccis[i] <- paste(ccis_i, collapse=',')
  cci_counts[i] <- length(ccis_i)
  cci_scores[i] <- paste(cci_scores_i, collapse=',')
}
hist(cci_counts)

# Creating & saving DF #
cci_df <- as.data.frame(list('n_cci_sig_scsignalr'=cci_counts,
                             'ccis_scsignalr'=ccis, 
                             'ps_scsignalr'=cci_scores) 
                        )
rownames(cci_df) <- st_lrs
write.table(cci_df, file=paste0(out_dir, 'scsignalr_summary.txt'), 
            sep='\t', quote=F)











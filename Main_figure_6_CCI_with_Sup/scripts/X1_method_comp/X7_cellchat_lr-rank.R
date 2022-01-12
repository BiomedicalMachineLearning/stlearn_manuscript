# (15/10/21) -> Examining the LR results from CellChat, in
# particular the ranking of LR pairs in-terms of the CCIs
# the LRs are mean to facilitate.
# 
# INPUT: * data/method_comp/breast/cellchat*
# OUTPUT: ... 

################################################################################
                        # Environment setup #
################################################################################
work_dir <- "/Users/uqbbalde/Desktop/Uni_Studies/projects/stLearn/"
setwd(work_dir)

Sys.setenv(RETICULATE_PYTHON = "~/opt/miniconda3/envs/STI/bin/python/")
library(CellChat)
library(stringr)
library(ggplot2)
library(reticulate)
sc <- import('scanpy')

data_dir <- 'data/method_comp/breast/'
data_dir2 <- '/Volumes/GML001-Q1851/Brad/'
st_db_dir <- '/Users/uqbbalde/Desktop/Uni_Studies/myPython/stlearn_latest/stLearn/stlearn/tools/microenv/cci/databases/'
out_dir <- 'data/method_comp/breast/'

data_set <- 'breast'

################################################################################
                         # Loading the data #
################################################################################
cellchat <- readRDS(paste0(data_dir, 'cellchat_', data_set, '.rds'))

# The stlearn LR pairs used, so can make sure comparing same thing #
# st_db <- read.table(paste0(st_db_dir, 'connectomeDB2020_lit.txt'), 
#                     sep='\t', header=T)
# st_lrs <- str_c(st_db[,1], st_db[,2], sep='_')

# Actually, need to load the equivalent LRs from the breast data #
breast_ad <- sc$read_h5ad(paste0(data_dir2, 'breast_LR&CCIResults.h5ad'))
lr_summary <- as.data.frame( breast_ad$uns[['lr_summary']] )
st_lrs <- rownames( lr_summary )

################################################################################
       # Extracting the no. of interaction per LR pair from CellChat #
################################################################################
##### Checking what stlearn lrs were tested in cellchat #####
cc_lrs <- rownames(cellchat@DB$interaction)
cc_intersect <- intersect(st_lrs, cc_lrs)
print(length(st_lrs)) # 2293, 549
print(length(cc_intersect)) # 2293, 549
# So did test all LR pairs

###### Seeing the celltype interaction counts of stlearn LRs in cell chat #####
df.net <- subsetCommunication(cellchat, thresh=1 #pval thresh, need all data.
                              ) # dataframe with one row per CCI.
ls <- df.net[,'ligand']
rs <- df.net[,'receptor']
cc_int_lrs <- str_c(ls, rs, sep='_')
df.net['lrs'] <- cc_int_lrs

cc_lrs_match <- match(cc_int_lrs, st_lrs)
cc_lrs_bool <- !is.na(cc_lrs_match)

# Subsetting to stlearn lrs for fair comparison #
st.df.net <- df.net[cc_lrs_bool,]
print(max(st.df.net['pval'])) # Making sure only contains sig interactions.
cc_st_lrs <- cc_int_lrs[cc_lrs_bool]

# Count the no. of interactions per stlearn LR, get all LR p-values & CCIs #
st_lr_ints <- integer(length(st_lrs))
cc_lr_ps <- character(length(st_lrs))
lr_ccis <- character(length(st_lrs)) # Used to check spatial enrichment afterward.
for (i in 1:length(st_lrs)){
  st_lr <- st_lrs[i]
  st_lr_bool <- cc_st_lrs==st_lr
  sig_bool <- st.df.net[,'pval'] < .05
  sig_st_lr_bool <- st_lr_bool & sig_bool
  st_lr_ints[i] <- sum(sig_st_lr_bool)
  pvals <- st.df.net[sig_st_lr_bool,'pval']
  pvals[pvals==0] <- .01 #due to n_perms = 100
  cc_lr_ps[[i]] <- paste(-log10(pvals), collapse=',')
  ccis <- str_c(st.df.net[sig_st_lr_bool,'source'], 
                st.df.net[sig_st_lr_bool,'target'], sep='--')
  print(length(ccis)==st_lr_ints[i])
  lr_ccis[i] <- paste(ccis, collapse=',')
}
n_top <- 100
#print(st_lr_ints[order(-st_lr_ints)][1:50])
print(st_lrs[order(-st_lr_ints)][1:n_top])

############# Looking at the intersection for the top ranked LRs ###############
top_intersect <- intersect(st_lrs[1:n_top], st_lrs[order(-st_lr_ints)][1:n_top])
print(length(top_intersect))
print(top_intersect) # 27 overlap, that seems quite reasonable actually !!!!!!!

################################################################################
      # Extracting the no. of CCI interactions for stLearn CCI #
################################################################################
# TODO, should add this into the lr_summary when running CCI by default!!!!!!!!!
lr_pval_dfs <- breast_ad$uns[['per_lr_cci_pvals_cell_type']]
lr_int_counts <- integer(length(st_lrs))
st_lr_ps <- character(length(st_lrs))
st_lr_ccis <- character(length(st_lrs))
for (i in 1:length(st_lrs)) {
  lr <- st_lrs[i]
  if (lr %in% names(lr_pval_dfs)) {
    pval_df <- py_to_r( lr_pval_dfs[[lr]] )
    sig_bool <- pval_df<.05
    n_ints <- sum(sig_bool)
    lr_int_counts[i] <- n_ints
    #st_lr_ps[i] <- paste(unlist( pval_df ), collapse=',')
    st_ccis <- c()
    pval_str <- ''
    for (k in 1:dim(sig_bool)[1]) {
      for (l in 1:dim(sig_bool)[2]) {
        if (sig_bool[k,l]) {
          st_ccis <- c(st_ccis, 
                      paste(rownames(pval_df)[k], colnames(pval_df)[l], sep='--'))
          pval <- pval_df[k,l]
          if (pval == 0) {pval <- .01}
          pval <- -log10(pval)
          if (pval_str==''){ pval_str <- paste0(pval_str, pval) }
          else { pval_str <- paste(pval_str, pval, sep=',') }
        }
      }
    }
    st_lr_ccis[i] <- paste(st_ccis, collapse=',')
    st_lr_ps[i] <- pval_str
    print(length(st_ccis)==lr_int_counts[i])
  } else {
    lr_int_counts[i] <- 0
  }
}

# TODO figure out why there are less interaction dataframes than lrs...
print(sum(breast_ad$uns[['lr_summary']][,'n_spots_sig'] > 3))
# Only tested for LRs with a minimum no. of significant spots is why!

######## Adding the CCI interactions to the lr_summary!!! #########
lr_summary['cci_tested'] <- lr_summary[,'n_spots_sig'] > 3
lr_summary['n_cci_sig'] <- lr_int_counts
lr_summary['n_cci_sig_cellchat'] <- st_lr_ints
lr_summary['ps'] <- st_lr_ps
lr_summary['ps_cellchat'] <- cc_lr_ps
lr_summary['ccis'] <- st_lr_ccis
lr_summary['ccis_cellchat'] <- lr_ccis

####### Saving this for visualisation in python ########
write.table(lr_summary, sep='\t', quote=F,
            file=paste0(out_dir, 'cci_summary.txt'))










# Index
Detailed documentation of each script in terms of input/outputs.

scripts/X1_method_comp/
        -> Contains the scripts to compare stLearn CCI with other methods.
            In this comparison, a consistent LR database from NATMI was used
            for each tool (with the exception of cellphonedb, which has a bug 
            for updating the database). 

    X1_lr-analysis_breast.py -> Performs the LR testing on other 
                            visium data which will be in-use for the paper;
                            i.e. the breast cancer & the skin dataset.
                            NOTE: this is the local version, edited extensively
                                                                     on tinaroo.
                            
                            INPUT: * GML001-Q1851/Jon/
                                         Human_Breast_Cancer_Block_A_Section_1/*
                            OUTPUT: * GML001-Q1851/Brad/
                                          bgPerLR_lrsubset_spotsubset_noBgs.h5ad
    
    X2_stlearn_cci_breast.py -> Perform stlearn CCI permutation & creating the 
                CCI networks and spatial CCI visualisations from stlearn analysis.
                (Plots the GPC3_IGF1R results for main figure).
                
                INPUT:  * /Volumes/GML001-Q1851/Brad/
                                   breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad
                        * data/breast/label_transfer_bc.csv     
                OUTPUT: * plots/X1_method_comp/stlearn_breast_*
                        * /Volumes/GML001-Q1851/Brad/breast_LR&CCIResults.h5ad
    
    X3_stlearn_diagnostics.py -> Visualisation of the stlearn diagnostic plots 
                                for the breast data, namely LR ranked lists, 
                                correlation between LR expression/abundance and 
                                significance, and correlation of cell type 
                                abundance & interactions.

                INPUT: * /Volumes/GML001-Q1851/Brad/breast_LR&CCIResults.h5ad
                OUTPUT: * plots/X1_X3_diagnostics/*
                        
    X4_cellchat_cci_breast.R -> Predicting CCIs using CellChat for comparison with stlearn; point
                is to show that CellChat will predict interactions for cell
                types which never spatially co-occur.

                INPUT: * /Volumes/GML001-Q1851/Brad/
                               breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad
                        * data/breast/label_transfer_bc.csv
                OUTPUT: * data/breast/cellchat*
                        * plots/X1_method_comp/cellchat*
    
    X5_natmi_breast.py -> Runs the NATMI pipeline.
    
                        INPUT: * /Volumes/GML001-Q1851/Brad/
                                                       breast_LR&CCIResults.h5ad
                        OUTPUT: * data/breast/natmi/
    
    X6_cellchat_lr-rank.R -> Examining the LR results from CellChat, 
                        in particular the ranking of LR pairs in-terms of the 
                        CCIs the LRs are mean to facilitate.
                        
                  INPUT: * data/breast/cellchat*
                         * /Volumes/GML001-Q1851/Brad/breast_LR&CCIResults.h5ad
                  OUTPUT: * data/breast/cci_summary.txt    
               
    X7_singlecellsignalr_breast.R -> Runs SingleCellSignalR on the breast
                        data using the cell type labels, then saves the overall
                        CCI information inferred and by what LR pairs. 
                        
                        INPUT: * /Volumes/GML001-Q1851/Brad/
                                                       breast_LR&CCIResults.h5ad
                        OUTPUT: * data/breast/singlecellsignalr/*
    
    X8_cellphonedb_cci_breast.py -> Prepares the breast data to be run with 
                                                         cellphonedb, then runs.
                  Here for tutorial:
                  https://github.com/Teichlab/cellphonedb
                  Apparently inputting pre-normalised data:
                  https://github.com/Teichlab/cellphonedb/issues/12

                INPUT: * /Volumes/GML001-Q1851/Brad/
                                 breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad
                       * data/label_transfer_bc.csv
                OUTPUT: * data/breast/cellphonedb/*
                        * data/breast/cellphonedb/db/*
                                            -> Running native database
                        * data/breast/cellphonedb/st_db/*
                                            -> Running stlearn database               

    X9_squidpy_cpdb_breast.py -> Trying the squidpy implimentation of cpdb, 
            since from the API
            (https://squidpy.readthedocs.io/en/latest/api/squidpy.gr.ligrec.html#squidpy.gr.ligrec)
            can clearly update the database. Following tutorial here:
            https://squidpy.readthedocs.io/en/latest/auto_examples/graph/compute_ligrec.html

            INPUT: * /Volumes/GML001-Q1851/Brad/breast_LR&CCIResults.h5ad
            OUTPUT: * data/breast/squidpy/*
    
    X10_method_comp_vis.py -> Multi-method comparison of CCI ranking based 
                            on no. of significant LR pairs. Will focus on pair
                            which is highly ranked in the other methods but low-
                            ranked by stLearn as example, and show the squidpy 
                            colocalisation analysis of the two cell types.
                            
              INPUT: * data/method_comp/breast/cci_summary.txt
                     * data/method_comp/breast/singlecellsignalr/
                                                       scsignalr_summary.txt
                     * data/method_comp/breast/natmi/natmi_cci_summary.txt
                     * data/method_comp/breast/squidpy/squidpy_cci_summary.txt
              OUTPUT: * plots/multi_comp_v2/*

scripts/X2_method_comp_clusters/
        -> Contains scripts for comparing the methods when using breast cluster
            information; in this comparison, the native LR database of each tool
            is used. 

    X1_stlearn_cci_cluster.py -> Adds in the stLearn clustering information to 
                                    the anndata generate from SME clustering.

         INPUT:  * /Volumes/GML001-Q1851/Brad/
                                   breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad
                 * data/breast/bc_clusters.csv
         OUTPUT: * /Volumes/GML001-Q1851/Brad/breast_ClusterCCIResults.h5ad
                 * plots/X2_method_comp_clusters/
    
    X2_cellphonedb_cci_cluster.py -> Prepares the breast data to be run with 
                           cellphonedb, then runs for the breast cluster annots.

        Refer here for Laura reference on cellphonedb:
            https://github.com/ellefeg/RNASeq/blob/master/scRNASeq_101/
                                                           Seurat2CellphoneDB.md
        Also here for tutorial:
            https://github.com/Teichlab/cellphonedb
        Apparently inputting pre-normalised data:
            https://github.com/Teichlab/cellphonedb/issues/12

        INPUT: * /Volumes/GML001-Q1851/Brad/breast_ClusterCCIResults.h5ad
        OUTPUT: * data/breast/cluster/cellphonedb/*
                * data/breast/cluster/cellphonedb/db/*

    X3_cellchat_cci_cluster.R -> Predicting CCIs using CellChat for 
                                 comparison with stlearn; point is to show 
                                 that CellChat will predict interactions for 
                                 cell types which never spatially co-occur.
    
       INPUT: * /Volumes/GML001-Q1851/Brad/breast_ClusterCCIResults.h5ad
       OUTPUT: * data/breast/cluster/cellchat*

    X4_singlecellsignalr_cluster.R -> Runs SingleCellSignalR on the breast data 
                                    using the cluster, then saves the overall 
                                    CCI information inferred and by what LR 
                                    pairs.
               INPUT: * /Volumes/GML001-Q1851/Brad/breast_ClusterCCIResults.h5ad
               OUTPUT: * data/breast/cluster/singlecellsignalr*

    X5_natmi_cluster.py -> Runs the NATMI pipeline on the breast clusters.

        INPUT: * /Volumes/GML001-Q1851/Brad/breast_ClusterCCIResults.h5ad
        OUTPUT: * data/breast/cluster/natmi/

    X6_squidpy_cpdb_cluster.py -> Trying the squidpy implimentation of cpdb, 
                                since from the API. Following tutorial here:
        https://squidpy.readthedocs.io/en/latest/auto_examples/graph/compute_ligrec.html

            INPUT: * /Volumes/GML001-Q1851/Brad/breast_ClusterCCIResults.h5ad
            OUTPUT: * data/breast/cluster/
            
    X7_method_cci_vis.py -> Visualises the results from running the 
                            cluster-based CCI with them different methods.

             INPUT: * /Volumes/GML001-Q1851/Brad/breast_ClusterCCIResults.h5ad
             OUTPUT: * plots/X2_method_comp_clusters/*_chordplot.pdf

scripts/X3_skin_analysis/
    -> Analysis of the skin ST-seq data.
    
    X1_preprocess_skin.py -> Performs the preprocessing on the skin data.

                            INPUT: * /Volumes/GML001-Q1851/Jon/C1_skin_raw/
                            OUTPUT: * /Volumes/GML001-Q1851/Brad/skin_analysis/

    X2_run_skin_LRs.py -> Loads in the preprocessed skin data & performs the LR 
                            analysis on this using literature reviewed LRs.
                     NOTE: run on HPC, should be reasonably fast locally though.

            INPUT:  * /Volumes/GML001-Q1851/Brad/skin_analysis/skin_ppd.h5ad
            OUTPUT: * /Volumes/GML001-Q1851/Brad/skin_analysis/skin_lrs.h5ad

    X3_gen_naiive_bg.py -> Generates a uniform background of randomly selected 
                                                                          pairs.

            INPUT: * /Volumes/GML001-Q1851/Brad/skin_analysis/skin_ppd.h5ad
            OUTPUT: * /Volumes/GML001-Q1851/Brad/skin_analysis/naiive_bg.pkl

    X4_vis_IL34-CSF1R_results.py -> Visualising the IL34-CSF1R results and the 
                                            background statistics as an example.

            INPUT: * /Volumes/GML001-Q1851/Brad/skin_analysis/skin_lrs.h5ad
                   * /Volumes/GML001-Q1851/Brad/skin_analysis/naiive_bg.pkl
            OUTPUT: * plots/X3_skin_analysis/IL34_CSF1R*

    X5_vis_high-pair_results.py -> Visualising the highest expressing LR 
                                   pair (COL1A1-DDR1) results and the
                                   background statistics as an example.

            INPUT: * /Volumes/GML001-Q1851/Brad/skin_analysis/skin_lrs.h5ad
                   * /Volumes/GML001-Q1851/Brad/skin_analysis/naiive_bg.pkl
            OUTPUT: * plots/X3_skin_analysis/COL1A1_DDR1*

scripts/X4_seqFish/
    -> Analysis of the mouse SVZ seqFISH+ data.

    X1_create_anndata.py -> Creating anndata object from the seqFISH data.

                    INPUT: * /Volumes/GML001-Q1851/Brad/seqFISH/cortex_svz_*.csv
                    OUTPUT: * /Volumes/GML001-Q1851/Brad/seqFISH/svz.h5ad

    X2_run_st-lr-cci.py -> Running the stlearn LR-CCI analyses on the slideSeq.

         INPUT: * /Volumes/GML001-Q1851/Brad/slideSeq/slideSeq/hipp.h5ad
         OUTPUT: * /Volumes/GML001-Q1851/Brad/slideSeq/slideSeq/hipp_LR-CCI.h5ad

    X3_vis_st-lr-cci.ipynb -> Visualises results from stlearn LR-CCI pipeline on 
                                the SVZ SeqFish+ data.

        INPUT: * /Volumes/GML001-Q1851/Brad/seqFISH/svz_LR-CCI.h5ad
        OUTPUT: * plots/X4_seqFISH/*

    X4_grid_data.py -> Runs gridding function on the seqFISH+ data, essentially 
                        pseudobulking by neighbourhoods to improve speed with 
                        running CCI. Then runs stlearn LR-CCI analysis.

        INPUT: * /Volumes/GML001-Q1851/Brad/seqFISH/svz.h5ad
        OUTPUT: * /Volumes/GML001-Q1851/Brad/seqFISH/svz_gridded.h5ad

    X5_vis_grid.ipynb -> Visualises the results from the gridding analysis on 
                                                           the SVZ seqFISH+ data

        INPUT: * /Volumes/GML001-Q1851/Brad/seqFISH/svz_gridded.h5ad
        
        OUTPUT: * plots/X4_seqFISH/grid_*

scripts/X5_slideSeq/
    -> Analysis of the mouse hippocampus slide-seq data.

    X1_create_anndata.R -> Loads in the slide-seq data from Seurat & converts to 
                            anndata, subsequently saving as .h5ad for later 
                            loading in python. Also performs the Seurat label 
                            transfer onto this data as per the vignette, using 
                            downloaded reference scRNA-seq data.
 
           INPUT: * data/slideSeq/mouse_hippocampus_reference.rds
                     -> Downloaded from here:
                     https://www.dropbox.com/s/cs6pii5my4p3ke3/mouse_hippocampus_reference.rds?dl=0
                  * Uses Seurat to download data of hippocampus.
           OUTPUT:
                  * /Volumes/GML001-Q1851/Brad/slideSeq/hipp.h5ad

    X2_format_anndata.py ->  Loads in the slide-seq data & subsequently formats 
                the spatial information and normalises data using scanpy 
                workflow prior to downstream lr-cci analysis.

             INPUT: * /Volumes/GML001-Q1851/Brad/slideSeq/hipp.h5ad
             OUTPUT: * /Volumes/GML001-Q1851/Brad/slideSeq/hipp.h5ad

    X3_rn_st-lr-cci.py -> Running the stlearn LR-CCI analyses on the slideSeq.

         INPUT: * /Volumes/GML001-Q1851/Brad/slideSeq/slideSeq/hipp.h5ad
         OUTPUT: * /Volumes/GML001-Q1851/Brad/slideSeq/slideSeq/hipp_LR-CCI.h5ad

    X4_vis_st-lr-cci.ipynb -> Visualising results from stLearn LR-CCI analysis 
                                                           on the slideSeq data.

        INPUT: * /Volumes/GML001-Q1851/Brad/slideSeq/hipp_LR-CCI.h5ad
        OUTPUT: * plots/X5_slideSeq/*
                * /Volumes/GML001-Q1851/Brad/slideSeq/hipp_LR-CCI_clustered.h5ad
                
    X5_grid_st-lr-cci.py -> Runs stLearn CCI analysis on the slide-seq data 
                                                                 after gridding.

        INPUT: * /QRISdata/Q1851/Brad/slideSeq/hipp_LR-CCI_clustered.h5ad
        
        OUTPUT: * /QRISdata/Q1851/Brad/slideSeq/hipp_gridded_45-45_LR-CCI.h5ad
                * plots/X5_slideSeq/
                
    X6_vis_gridded-st-lr.ipynb -> Visualising the results from gridding the data 
                                                    & performing the LR analyis.

        INPUT: * /Volumes/GML001-Q1851/Brad/slideSeq/
                                                  hipp_gridded_45-45_LR-CCI.h5ad
        OUTPUT: * plots/X5_slideSeq/*     

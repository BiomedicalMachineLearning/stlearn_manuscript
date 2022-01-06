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









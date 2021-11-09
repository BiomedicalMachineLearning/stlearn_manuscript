# Index
Detailed documentation of each script in terms of input/outputs.

scripts/X1_method_comp/

    X1_lr-analysis_breast.py -> Performs the LR testing on other 
                            visium data which will be in-use for the paper;
                            i.e. the breast cancer & the skin dataset.
                            NOTE: this is the local version, edited extensively
                                                                     on tinaroo.
                            
                            INPUT: GML001-Q1851\Jon\
                                         Human_Breast_Cancer_Block_A_Section_1/*
                                   GML001-Q1851\Jon\C1_skin_raw/*  
                            OUTPUT: /scratch/user/uqbbalde/stLearn/
                                                            data/bg_eval/
                                                              breast/* or skin/*
    
    X2_stlearn_cci_breast.py -> Creating the equivalent plots as CellChat
                for stlearn, in order to do a side-by-side comparison, which 
                if considered in conjunction with the colocalisation results &
                the spatial plots of location of the interactions, should show 
                how stlearn can reduce the false-positive-rate by taking into
                account spatial information.
                (Plots the GPC3_IGF1R results for main figure).
                
                INPUT:  * /Volumes/GML001-Q1851/Brad/
                               breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad
                        * data/label_transfer_bc.csv     
                OUTPUT: * plots/method_comp/stlearn_breast_*
                        * /Volumes/GML001-Q1851/Brad/breast_LR&CCIResults.h5ad
                        
    X3_stlearn_lr-spatial_breast.py -> For selected LRs, visualising 
                the spatial LR scores along with the cell types as a way to show
                that the stlearn result is correct, where the CellChat 
                prediction is incorrect because it lacks the spatial context.
                
                INPUT:  * /Volumes/GML001-Q1851/Brad/breast_LR&CCIResults.h5ad
                OUTPUT: * plots/method_comp/stlearn_breast_spatial_*
                
    X4_cellchat_cci_breast.R -> Predicting CCIs using CellChat for 
                comparison with stlearn; point is to show that CellChat will
                predict interactions for cell types which never spatially 
                co-occur.
                
                INPUT: * /Volumes/GML001-Q1851/Brad/
                               breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad
                        * data/label_transfer_bc.csv
                OUTPUT: * data/method_comp/breast/cellchat*
                        * plots/method_comp/cellchat*
                        
    X5_singlecellsignalr_breast.R -> Runs SingleCellSignalR on the breast
                        data using the cell type labels, then saves the overall
                        CCI information inferred and by what LR pairs. 
                        
                        INPUT: * /Volumes/GML001-Q1851/Brad/
                                                       breast_LR&CCIResults.h5ad
                        OUTPUT: * data/method_comp/breast/singlecellsignalr/*
                        
    X6_natmi_breast.py -> Runs the NATMI pipeline.
    
                        INPUT: * /Volumes/GML001-Q1851/Brad/
                                                       breast_LR&CCIResults.h5ad
                        OUTPUT: * data/method_comp/breast/natmi/
    
    X7_cellchat_lr-rank.R -> Examining the LR results from CellChat, 
                        in particular the ranking of LR pairs in-terms of the 
                        CCIs the LRs are mean to facilitate.
                        
                  INPUT: * data/method_comp/breast/cellchat*
                         * /Volumes/GML001-Q1851/Brad/breast_LR&CCIResults.h5ad
                  OUTPUT: * data/method_comp/breast/cci_summary.txt           
    
    X8_method_comp_vis.py -> Multi-method comparison of CCI ranking based 
                            on no. of significant LR pairs. Will focus on pair
                            which is highly ranked in the other methods but low-
                            ranked by stLearn as example, and show the squidpy 
                            colocalisation analysis of the two cell types. Then
                            repeat this for every CCI and predicted LR pair in 
                            order to plot distributions of squidpy spatial 
                            AUC to show stLearn enriches for cell types spatially
                            when compared to other methods. 
                            
                  INPUT: * data/method_comp/breast/cci_summary.txt
                         * data/method_comp/breast/singlecellsignalr/
                                                           scsignalr_summary.txt
                     * data/method_comp/breast/natmi/natmi_cci_summary.txt
                  OUTPUT: * plots/multi_comp_v2/*
    
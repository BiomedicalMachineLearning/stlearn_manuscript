# Code for reproducing stLearn paper

## Abstract: 
Spatial Transcriptomics (ST) is an emerging technology that adds spatial dimensionality to the genome-wide transcriptional profiling of cells in undissociated tissues. From a sample, a typical ST technology can generate three spatial data types, namely morphological imaging, physical distance between spatial data points, and gene expression values. We developed three new algorithms, collectively referred to as stLearn, utilising all the three data types to study tissue biology. First, stLearn corrects for technical noise in spatial sequencing data by using tissue image features (optional), physical distance, and gene experession profiles (SME method). This way, stLearn can also impute missing data that can increase tissue coverage of spatial gene expression data. SME significantly improves clustering analysis to find specific spatial patterns not detectable by existing methods. Second, we present a new method, pseudo-time-space (PSTS), to model the spatiotemporal relationship of cellular transcriptional states across a tissue. Using immunofluorescence imaging and engineered mouse model of central nervous system injury, we validated PSTS results in reconstructing the spatial trajectory of microglia activation following brain insult. We also assessed the diagnostic potential of PSTS in studying breast cancer progression. Third, we developed a spatial interaction test, which integrates ligand-receptor expression and spatial neighbourhood information to find highly interactive regions within a tissue across thousands of ligand-receptor pairs. We thoroughly benchmark, assess false discovery, and validate the interactions in skin and breast cancer tissues. Together, the three algorithms that we developed, as implemented in the comprehensive and fast stLearn software, allow for the elucidation of biological processes within healthy and diseased tissues.
 
## Main parts:
 
####   Mapping cell types/

- [Main Figure 2, 3 - clustering](https://github.com/BiomedicalMachineLearning/stlearn_manuscript/blob/main/Main_figure_2_SME/stSME_clustering.ipynb) - [ipynb]

- [Main Figure 2, 3 - comparison](https://github.com/BiomedicalMachineLearning/stlearn_manuscript/blob/main/Main_figure_2_SME/stSME_comparison.ipynb) - [ipynb]

####  Spatio-temporal trajectories/
                    
- [Figure 4 - Traumatic brain injury](https://github.com/BiomedicalMachineLearning/stlearn_manuscript/blob/main/Main_figure_4_PSTS_TBI/TBI_related_fig.ipynb) - [ipynb]
- [Figure 5 - Breast cancer progression](https://github.com/BiomedicalMachineLearning/stlearn_manuscript/blob/main/Main_figure_5_PSTS_Cancer_progression/PSTS_tutorial_BCBA.ipynb) - [ipynb]
- [Figure 5 - Neuronal development](https://github.com/BiomedicalMachineLearning/stlearn_manuscript/blob/main/Main_figure_5_PSTS_Neuronal_development/PSTS_Embryo.ipynb) - [ipynb]

####  Cell-cell interactions/

- [Figure 6 - CCI](https://github.com/BiomedicalMachineLearning/stlearn_manuscript/tree/main/Main_figure_6_CCI_with_Sup) - [tutorial]

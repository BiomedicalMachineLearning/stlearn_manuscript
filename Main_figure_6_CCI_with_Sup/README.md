# StLearn Reproduce
Scripts to reproduce key results from the stLearn CCI benchmarking analysis. 

NOTE that assumes the working directory for all scripts is within the cloned
                                                               github directory:

    stlearn_reproduce_results/

To install python package dependencies, ideally in a conda environment:

    pip install -r requirements.txt
    
Will also need to install development branch of stlearn:

    git clone https://github.com/BiomedicalMachineLearning/stLearn.git
    git fetch
    git checkout brad_dev
    
Then add the path to the cloned github to your system python path.

For information regarding R dependencies, please see:

    docs/rsession_info.md

## Global index:

    docs/
       
       index.md -> Contains granular documentation of each scripts in scripts/
       
    scripts/
    
        X1_method_comp/ -> Contains the scripts to compare stLearn CCI with
                            other methods. 
                            
        X2_method_comp_clusters/ -> Contains scripts for comparing the methods 
                                    when using breast cluster information.
        
        X3_skin_analysis/ -> Scripts for analysis of the skin ST-seq data.
                
        utils/ -> Utility scripts used across many scripts.

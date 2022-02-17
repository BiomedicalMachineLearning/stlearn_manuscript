"""
Runs stLearn LR-CCI on the simulated spatial data.

      INPUT: * data/sim_data/spatialsim_v2.h5ad
      OUTPUT: * data/sim_data/spatialsim_v2.h5ad
                                             -> Adds the stlearn-LR-CCI results!
"""

################################################################################
                        # Environment Setup #
################################################################################
import numpy as np
import scanpy as sc
import stlearn as st

out_dir = 'data/sim_data/'
data_dir = out_dir

################################################################################
                        # Loading the data #
################################################################################
spatial_data = sc.read_h5ad(data_dir+'spatialsim_v2.h5ad')

################################################################################
                    # Running stLearn LR-CCI #
################################################################################
## First getting appropriate normalisation ##
spatial_counts = spatial_data.raw.to_adata().X.astype(int)
spatial_data.X = spatial_counts
sc.pp.normalize_total(spatial_data)

## Loading LRs ##
lrs = st.tl.cci.load_lrs()

## Running the LR analysis ##
spatial_data = spatial_data[:,
            [gene for gene in spatial_data.var_names if '_' not in gene]].copy()
st.tl.cci.run(spatial_data, lrs,
                  min_spots = 50,#50, #Filter out any LR pairs with no scores for less than min_spots
                  distance=0, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=100,#10000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=2, # Number of CPUs for parallel. If None, detects & use all available.
                  )

# spatial_data.write_h5ad(out_dir+'spatialsim_v2.h5ad', compression='gzip')

decon = spatial_data.obsm['deconvolution']
max_indices = np.apply_along_axis(np.argmax, 1, decon.values)
cell_type = [decon.columns.values[index] for index in max_indices]
spatial_data.obs['cell_type'] = cell_type
spatial_data.obs['cell_type'] = spatial_data.obs['cell_type'].astype('category')
spatial_data.uns['cell_type'] = decon

# spatial_data.write_h5ad(out_dir+'spatialsim_v2.h5ad', compression='gzip')

## Running CCI ##
st.tl.cci.run_cci(spatial_data, 'cell_type', # Spot cell information either in data.obs or data.uns
                  min_spots=0, # Minimum number of spots for LR to be tested.
                  spot_mixtures=True, # If True will use the label transfer scores,
                                      # so spots can have multiple cell types if score>cell_prop_cutoff
                  cell_prop_cutoff=0.1, # Spot considered to have cell type if score>0.2
                  sig_spots=False, # Only consider neighbourhoods of spots which had significant LR scores.
                  n_perms=10#100 # Permutations of cell information to get background, recommend ~1000
                 )

spatial_data.write_h5ad(out_dir+'spatialsim_v2.h5ad', compression='gzip')



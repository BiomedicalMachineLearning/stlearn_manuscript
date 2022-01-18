import stlearn as st
import numpy
import os

def setUp():
    os.chdir('/Users/uqbbalde/Desktop/Uni_Studies/projects/stLearn/')

def load():
    # read in visium dataset downloaded from: support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Breast_Cancer_Block_A_Section_2
    data_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/data/Visium8_A1_Hybrid_treated/'
    data = st.Read10X(data_dir)
    data.var_names_make_unique()
    data.var_names = numpy.array([var_name.replace('_', '-')
                                  for var_name in data.var_names])
    # Adding in the label transfer information #
    # TODO

    #st.pp.filter_genes(data, min_cells=int(0.00 * data.n_obs))
    st.pp.normalize_total(data)
    # st.pp.log1p(data)
    # st.pp.scale(data)
    #st.em.run_pca(data, n_comps=50)

    return data

def load_breast(rdm_dir='/Volumes/GML001-Q1851/Jon/', add_raw=False):
    """ Loads the breast-cancer data; need to have RDM mounted beforehand.
    """
    data_dir = rdm_dir+'Human_Breast_Cancer_Block_A_Section_1/'
    data = st.Read10X(data_dir)
    st.add.image(adata=data,
                 imgpath=data_dir+"spatial/tissue_hires_nobg.png",
                 library_id="V1_Breast_Cancer_Block_A_Section_1", visium=True)
    # preprocessing
    st.pp.filter_genes(data, min_cells=3)

    # Adding raw counts #
    if add_raw:
        data.raw = data

    # Normalising to total counts #
    st.pp.normalize_total(data)

    return data


def load_skin(rdm_dir='/Volumes/GML001-Q1851/Jon/'):
    """ Loads skin cancer data, see \
    scripts/jon_code/stlearn-cci_rank-skin-inhouse-between.ipynb for details, this \
                                                          was used in the paper.
    """
    data_dir = rdm_dir+'C1_skin_raw/'
    data = st.Read10X(data_dir)
    st.add.image(adata=data,
                 imgpath=data_dir+"spatial/tissue_hires_image.png",
                 library_id="Visium_LP2_map_C1_trimTSOpolyA", visium=True)
    st.pp.filter_genes(data, min_cells=3)
    st.pp.normalize_total(data)

    return data




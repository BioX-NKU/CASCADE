import h5py
import torch
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from .utils import *
    
def makesnap(train_path,w_path,f_path):
    """
    Generate snap files of different cell types from data in the form of h5ad files.
    
    Parameters
    --------------
    train_path:
        path to the h5ad file of the train dataset you want to read. 
    w_path:
        path to the output snap files.
    f_path:
        path to the h5py file.
    """
    
    setup_seed(2022)
    traindata = sc.read(train_path)
    traindata_concat = traindata.copy()
    celltypes = traindata_concat.obs.cell_type.value_counts().index
    typelist = pd.DataFrame(traindata_concat.obs.cell_type.value_counts())
    typelist.to_csv(w_path+'typelist.csv')

    f = h5py.File(f_path, 'r')
    f.keys()
    
    for celltype in celltypes:
        
        f1 = h5py.File(w_path+'%s.snap' % (celltype), 'w')
        sub_set = traindata_concat[np.where(traindata_concat.obs.cell_type.values == celltype)[0],].copy()

        data_use = sub_set.X
        data_use_sparse = sparse.csr_matrix(data_use)
        idx = data_use_sparse.tocoo().row + 1
        idy = data_use_sparse.tocoo().col + 1
        count = np.array(data_use_sparse.data, dtype='uint32')
        f1['AM/5000/count'] = count
        f1['AM/5000/idx'] = idx
        f1['AM/5000/idy'] = idy
        f1['AM/5000/binChrom'] = f['AM/5000/binChrom'][:max(idy)]
        f1['AM/5000/binStart'] = f['AM/5000/binStart'][:max(idy)]

        for index in list(f["FM"].keys()):
            f1["FM/" + index] = f['FM/' + index][()]

        del f1["FM/barcodeLen"]
        del f1["FM/barcodePos"]
        f1["FM/barcodeLen"] = np.array(range(data_use.shape[0])) + 10000
        f1["FM/barcodePos"] = (np.array(range(data_use.shape[0])) + 1) * 10
        f1.close()


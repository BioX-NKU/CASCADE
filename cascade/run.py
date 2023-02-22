from __future__ import print_function

import copy
import random
import scipy

import pandas as pd
import scanpy as sc
import numpy as np
import episcanpy as epi
import anndata as ad
import torch

import torch.nn.functional as F
from sklearn import preprocessing


from .model import *
from .data_preprocessing import *
from .utils import *

def run(train_path,test_path,work_path,device):

    setup_seed(2022)
    test_data = sc.read(test_path)
    test_data_copy = test_data.copy()
    
    # Read into the training set
    train_ = sc.read(train_path)
    train_data_truth = train_.copy()
    celltypes = train_data_truth.obs.cell_type.value_counts().index
    
    # Read into the simulation training set
    adata_ = sc.read(work_path+'%s.h5ad' % (celltypes[0]))
    adata_.obs['cell_type'] = celltypes[0]
    adata_concat = adata_.copy()
    for celltype in celltypes[1:celltypes.size + 1]:
        adata_ = sc.read(work_path+'%s.h5ad' % (celltype))
        adata_.obs['cell_type'] = celltype
        adata_concat = sc.concat([adata_concat, adata_])
    adata_concat.var = pd.DataFrame({'peak': test_data.var.peak[0:adata_concat.X.shape[1]]}, columns={'peak'})
    
    # Combine training set and simulation training set
    adatas = [adata_concat,train_data_truth]
    train_data = ad.concat(adatas)
    train_data.var = pd.DataFrame({'peak': test_data.var.peak[0:train_data.X.shape[1]]}, columns={'peak'})
    
    # Data preprocessing
    train_data_loader, test_data_loader,inp_size,number_of_classes,le = data_preprocessing(train_data, test_data_copy)

    # Make annotation
    
    cf_model = Classifier(output_dim=number_of_classes, input_size=inp_size)
    
    cf_model.train(train_data_loader,device)  
    pred_labels = cf_model.predict(test_data_loader,le,device)
    return pred_labels
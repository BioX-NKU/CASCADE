import numpy as np
import copy
import random
import scanpy as sc
import torch
from torch.autograd import Variable
from torch.utils.data import DataLoader
from sklearn import preprocessing
from scipy.sparse import isspmatrix
import torch.nn.functional as F
import episcanpy as epi
import anndata as ad
import scipy
import pandas as pd
from .utils import *

def data_preprocessing(train_data,test_data):
    '''
    The inputs are training datasets and test datasets with Anndata type,
    and it will return the preprocessing train_data_loader and tes_data_loader
    '''
    maskrate = 0.05
    setup_seed(2022)
    traindata_copy = train_data.copy()
    testdata_copy = test_data.copy()
    adatas_all = [traindata_copy,testdata_copy]
    adata_all = ad.concat(adatas_all,label='dataset')
    #filter_features
    fpeak = 0.01
    epi.pp.binarize(adata_all)
    epi.pp.filter_features(adata_all, min_cells=np.ceil(fpeak*adata_all.shape[0]))
    
    #TF-IDF transformation
    tfidf_res = tfidf3(adata_all.X.T).T
    adata_all.X = tfidf_res.copy()
    
    
    train_set = adata_all[adata_all.obs['dataset'] == '0']
    test_set = adata_all[adata_all.obs['dataset'] == '1']
    
    le = preprocessing.LabelEncoder()
    train_target = le.fit_transform(train_set.obs['cell_type'])
    train_set.obs['cluster'] = train_target
    
    #mask the training set
    train_set1 = mask(train_set,maskrate)
    train_set2 = mask(train_set,maskrate)
    train_set3 = mask(train_set,maskrate)
    train_set4 = mask(train_set,maskrate)
    train_set5 = mask(train_set,maskrate)
    train_sets = [train_set1,train_set2,train_set3,train_set4,train_set5]
    train_set = ad.concat(train_sets, label="dataset")

    
    inp_size = train_set.shape[1]
    number_of_classes = train_set.obs.cell_type.value_counts().shape[0]
    
    y_train = [int(x) for x in train_set.obs['cluster'].to_list()]

    try:
        train_dense = np.asarray(train_set.X.todense());
        test_dense = np.asarray(test_set.X.todense());
    except:
        train_dense = np.asarray(train_set.X);
        test_dense = np.asarray(test_set.X); 

    data_and_labels = []
    test_data = []
    for k in range(len(y_train)):
        data_and_labels.append([train_dense[k], y_train[k]])
        try:
            test_data.append([test_dense[k]])
        except:
            pass;

    test_data_loader = DataLoader(test_data, batch_size=32,shuffle=None, sampler=None,batch_sampler=None, num_workers=4, collate_fn=None, pin_memory=True)
        
    train_data_loader = DataLoader(data_and_labels, batch_size=32, shuffle=True, sampler=None,batch_sampler=None, num_workers=4, collate_fn=None, pin_memory=True)
    
    return train_data_loader, test_data_loader,inp_size,number_of_classes,le
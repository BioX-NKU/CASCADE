import numpy as np
import copy
import random
import scanpy as sc
import torch
from sklearn import preprocessing
from scipy.sparse import isspmatrix
import episcanpy as epi
import anndata as ad
import scipy
import pandas as pd
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.metrics import f1_score,accuracy_score,precision_score,cohen_kappa_score

def setup_seed(seed):
    """
    Set random seed.

    Parameters
    ----------
    seed
        Number to be set as random seed for reproducibility.

    """
    np.random.seed(seed)
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.cuda.manual_seed(seed)



def tfidf3(count_mat): 
    model = TfidfTransformer(smooth_idf=False, norm="l2")
    model = model.fit(np.transpose(count_mat))
    model.idf_ -= 1
    tf_idf = np.transpose(model.transform(np.transpose(count_mat)))
    return scipy.sparse.csr_matrix(tf_idf)

def mask(traindata,rate):
    traindata_copy = traindata.copy()
    if isspmatrix(traindata.X):
        non_zero = np.array(traindata_copy.X.nonzero())
        index_list = random.sample(list(range(non_zero.shape[1])), int(len(list(range(non_zero.shape[1])))*rate))
        mask_index_0 = non_zero[:,index_list]
        mask_index = tuple(mask_index_0.tolist())
        traindata_copy.X[mask_index] =0
    else:
        non_zero = np.array(np.where(traindata_copy.X != 0))
        index_list = random.sample(list(range(non_zero.shape[1])), int(len(list(range(non_zero.shape[1])))*rate))
        mask_index_0 = non_zero[:,index_list]
        mask_index = tuple(mask_index_0.tolist())
        traindata_copy.X[mask_index] =0
    return traindata_copy


def evaluate_metrics(gt_labels,pred_labels):
    '''
    input the ground truth labels and the predicted labels,
    and it will return the evaluate metrics
    '''
    acc = np.sum(gt_labels==pred_labels)/len(pred_labels)
    kappa = cohen_kappa_score(gt_labels, pred_labels)
    f1_macro = f1_score(gt_labels, pred_labels, average='macro')
    f1_weighted = f1_score(gt_labels, pred_labels, average='weighted')

    return acc,kappa,f1_macro,f1_weighted
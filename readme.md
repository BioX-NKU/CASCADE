[![PyPI](https://img.shields.io/pypi/v/epicascade.svg)](https://pypi.org/project/epicascade)
[![Downloads](https://pepy.tech/badge/epicascade)](https://pepy.tech/project/epicascade)


# Accurate Annotation for Differentiating and Imbalanced Cell Types in Single-cell Chromatin Accessibility Data

## Installation

Install CASCADE from PYPI

```
pip install epicascade
```

You can also install CASCADE from GitHub via

```
git clone git://github.com/BioX-NKU/CASCADE.git
cd CASCADE
python setup.py install
```

The dependencies will be automatically installed along with CASCADE.

## Quick Start

### Input:

**h5ad file** Files from the training set scCAS data and files from the scCAS data that need to be annotated.

### Output:

**pred_labels**: Array object which contains cell type annotation results.

### Using tutorial:

First, in order to simulate the dataset, we need to split the training set into different cell types and generate the corresponding snap files.
Note that in this step, we need a snap file as the underlying file to generate snap files, but the content of the file will not affect the content of the final generated snap file. This underlying file is available at [here](https://www.dropbox.com/s/muypr5w5ab7580p/GSE99172.snap?dl=0). Users can download it and provide the file path in ``cascade.makesnap``.

```python
import epicascade as cascade
cascade.makesnap(train_path,work_path,underlying_file_path)
```

Then, we use the simATAC package in R to simulate the data and obtain h5ad file of the simulated data.

```R
library(simATAC)
library(Matrix)
library(SingleCellExperiment)
library(zellkonverter)

simulation<-function(work_path){
    rate <- 6
    typelist = read.csv(paste(work_path,sprintf('typelist.csv',order),sep = ""))
    for (i in c(1:dim(typelist)[1])){
        goal_num = ceiling(colSums(typelist['cell_type'])/dim(typelist)[1]*rate)
        if(goal_num-typelist[i,2]<=5){
            n_cell <- 5
        }
        else{
            n_cell <- goal_num-typelist[i,2]
        }
        celltype <- typelist[i,1]
        print(celltype)
        print(n_cell)
        count <- getCountFromh5(paste(work_path,sprintf('%s.snap',celltype),sep = ""))
        object <- simATACEstimate(t(count))
        sim <- simATACSimulate(object, nCells=n_cell)
        print(sim)
        writeH5AD(sim, file=paste(work_path,sprintf('%s.h5ad',celltype),sep = ""))
    	}
    }
 simulation(work_path)
```

Finally, we can obtain the predicted labels of the test set via ``cascade.run``.

```python
pred_labels = cascade.run(train_path,test_path,work_path,device)
```
What's more, we can get the scores of prediction if we have the ground truth labels using ``cascade.evaluate_metrics``.
```python
acc,kappa,f1_macro,f1_weighted = cascade.evaluate_metrics(ground_truth_labels,pred_labels)
```

The source datasets are available at [here](https://www.dropbox.com/sh/hu1h340i70ktfc4/AACoq5-PlCTSJY5UvI-i4mAYa?dl=0). 


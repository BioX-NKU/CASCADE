# Accurate Annotation for Differentiating and Imbalanced Cell Types in Single-cell Chromatin Accessibility Data

## Installation

install CASCADE from PYPI

```
pip install CASCADE
```

Install CASCADE from GitHub

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

**array** Array object with containing cell type annotation results

### Use processes:

First, in order to simulate the dataset, we need to split the training set into different cell types and generate the corresponding snap files.

```python
import CASCADE
CASCADE.split.makesnap(train_path,work_path)
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
        count <- getCountFromh5(paste(work_path,sprintf('%s.snap',order,celltype),sep = ""))
        object <- simATACEstimate(t(count))
        sim <- simATACSimulate(object, nCells=n_cell)
        print(sim)
        writeH5AD(sim, file=paste(work_path,sprintf('%s.h5ad',order,celltype),sep = ""))
    	}
    }
```

Finally, we combine the training set with the simulation set, via the run function in CASCADE after data preprocessing together with the test set, and feed training set into model to train, and obtain the annotation results of the given test set.

```python
annotation = runmodel(train_path,test_path,work_path)
```

The source datasets are available at [here](https://www.dropbox.com/sh/hu1h340i70ktfc4/AACoq5-PlCTSJY5UvI-i4mAYa?dl=0). 


library(simATAC)
library(Matrix)
library(SingleCellExperiment)
library(zellkonverter)


simh5ad<-function(order,path){
    rate <- 6
    typelist = read.csv(paste(path ,sprintf('/fold%s/typelist.csv',order),sep = ""))
    typelist_all = read.csv(paste(path,'/typelist_all.csv',sep = ""))
    for (i in c(1:dim(typelist)[1])){
        goal_num = ceiling(colSums(typelist_all['cell_type'])/dim(typelist)[1]*rate)
        if(goal_num-typelist[i,2]<=5){
            n_cell <- 5
        }
        else{
            n_cell <- goal_num-typelist[i,2]
        }
        celltype <- typelist[i,1]
        print(celltype)
        print(n_cell)
        count <- getCountFromh5(paste(path,sprintf('/fold%s/%s.snap',order,celltype),sep = ""))
        object <- simATACEstimate(t(count))
        sim <- simATACSimulate(object, nCells=n_cell)
        print(sim)
        writeH5AD(sim, file=paste(path,sprintf('/fold%s/%s.h5ad',order,celltype),sep = ""))

}
    
    }
Args <- commandArgs(T)
simh5ad(Args[1],Args[2])

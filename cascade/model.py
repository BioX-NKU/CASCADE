import torch.nn as nn
import torch
import numpy as np
from torch.utils.data import DataLoader
import torch.utils.data
import torch.optim as optim
import torch.nn.functional as F

from sklearn.metrics import f1_score,accuracy_score,precision_score,cohen_kappa_score
torch.set_default_tensor_type(torch.DoubleTensor)


class Classifier(nn.Module):
    
    def __init__(self, output_dim=None, input_size=None):
       """
       A MLP module to predict cell type 
       consisting of an input layer,two hidden layers,an output layer.
       
       Parameters
       -------------
       input_size: the features in the training set
       output_dim: the number of cell types in the training set
       """
        super(Classifier, self).__init__();
        self.inp_dim = input_size;
        self.out_dim = output_dim;

        # feed forward layers
        self.classifier_sequential = nn.Sequential(
                                        nn.Linear(self.inp_dim, 128),
                                        nn.ReLU(),
                                        nn.Dropout(0.2),

                                        nn.Linear(128, 64),
                                        nn.ReLU(),
                                        nn.Dropout(0.2),

                                        nn.Linear(64, output_dim)
                                       
                                        )

    def forward(self, x):

        out = self.classifier_sequential(x);

        return out
    
    def train(train_data_loader,device='cpu',epochs=40):
        
        self.to(device)
        # define loss_function,optimizer,lr_scheduler
        cf_criterion = nn.CrossEntropyLoss().to(device)
        cf_optimizer = torch.optim.Adam(params=self.parameters(),
                                    lr=0.01,
                                    betas=(0.9, 0.999),
                                    eps=1e-08,
                                    weight_decay=0.005,
                                    amsgrad=False
                                    )
        cf_lr_scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer=cf_optimizer, gamma=0.7)
        
        
        for epoch in range(1, epochs + 1): 
            
            train_loss = 0                 
            for data, label in train_data_loader:
                cf_optimizer.zero_grad()
                data = data.to(device)
                label = label.to(device)
                output = self(data.double())
                loss = cf_criterion(output, label)
                loss.backward()
                cf_optimizer.step()
                train_loss += loss.item()*data.size(0)
            train_loss = train_loss/len(train_data_loader.dataset)

            cf_lr_scheduler.step()


    def predict(test_data_loader,le,device='cpu'):
        
        pred_labels = []

        #predict
        with torch.no_grad():
            for data in test_data_loader:
                data = data.to(device)
                output = self(data.double())
                preds = torch.argmax(output, 1)
                pred_labels.append(preds.data.cpu().numpy())
        pred_labels = le.inverse_transform(np.concatenate(pred_labels))
        
        return pred_labels

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
                
                
                
                
                
                
                
                
                

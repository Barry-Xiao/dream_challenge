#!/usr/bin/env python
# coding: utf-8

# ## Load packages

# In[ ]:


import pandas as pd
import numpy as np
import os

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader

from sklearn import metrics
from joblib import dump, load


# ## Set seed

# In[ ]:


random_seed = 8022022 # or any of your favorite number 
torch.manual_seed(random_seed)
torch.cuda.manual_seed(random_seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(random_seed)


# ## Setting outcome

# In[ ]:


outcome = "was_preterm"


# ## Loading data

# In[ ]:


metadata = pd.read_csv('/input/metadata/metadata.csv', delimiter=',',index_col = 'specimen')
metadata_feature = ['participant_id','project','collect_wk']
metadata_selected = metadata[metadata_feature].sort_index()


# In[ ]:


CST = pd.read_csv('/input/community_state_types/cst_valencia.csv', delimiter=',',index_col = 'specimen')['CST'].sort_index()


# In[ ]:


alpha = pd.read_csv('/input/alpha_diversity/alpha_diversity.csv', delimiter=',',index_col = 'specimen').sort_index()


# In[ ]:


taxonomy = pd.read_csv('/input/taxonomy/taxonomy_relabd.family.csv', delimiter=',',index_col = 'specimen')
tax_feature = list(pd.read_csv('/selected_feature/tax_family_preterm.csv', delimiter=',').feature_selected)
taxonomy_selected = taxonomy[tax_feature].sort_index()


# In[ ]:


phylotype = pd.read_csv('/input/phylotypes/phylotype_relabd.1e_1.csv', delimiter=',',index_col = 'specimen')
phylo_feature = list(pd.read_csv('./selected_feature/phylo_.1_preterm.csv', delimiter=',').feature_selected)
phylotype_selected = phylotype[phylo_feature].sort_index()


# In[ ]:


mydata = pd.concat([metadata_selected,CST,alpha,taxonomy_selected,phylotype_selected], axis = 1).copy()


# In[ ]:


mydata["project"] = mydata["project"].astype('category')
mydata["CST"] = mydata["CST"].astype('category')


mydata.dtypes


# ## Define functions

# In[ ]:


def tensor_generator(data):
    X_data = data[data.columns[2:]].copy()
    X_data['CST'] = X_data['CST'].cat.codes
    X_feature = X_data.to_numpy().astype('float32')
    X_group = data['participant_id'].astype('category').cat.codes.values.reshape(-1,1)
    input_X = torch.from_numpy(np.hstack((X_feature,X_group)))

    return input_X


# In[ ]:


class MLP(nn.Module):
    def __init__(self, input_dim, hidden_dim1,hidden_dim2,drop_out):
        

        
        #inherit from super class
        super(MLP, self).__init__()
        
        #define layers
        
        self.fc1 = nn.Linear(input_dim, hidden_dim1)
        self.fc2 = nn.Linear(hidden_dim1,hidden_dim2)
        self.fc3 = nn.Linear(hidden_dim2,2)
        self.dropout = nn.Dropout(drop_out)

        
    def forward(self, x):
        
        X_feature = x[:,:-1]
        X_group = x[:,-1].long()
        
        X_feature = torch.tanh(self.fc1(X_feature))
        X_feature = self.dropout(X_feature)
        X_feature = torch.tanh(self.fc2(X_feature))
        X_feature = self.dropout(X_feature)
        X_feature = torch.tanh(self.fc3(X_feature))

        X_feature = F.softmax(X_feature, dim = 1)

        
        M = torch.zeros(X_group.max()+1, len(X_feature))
        M[X_group, torch.arange(len(X_feature))] = 1
        M = F.normalize(M, p=1, dim=1)
        X_feature = torch.mm(M, X_feature)

        
        return X_feature


# In[ ]:


# torch.cuda.is_available() checks and returns a Boolean True if a GPU is available, else it'll return False
is_cuda = torch.cuda.is_available()

# If we have a GPU available, we'll set our device to GPU. We'll use this device variable later in our code.
if is_cuda:
    device = torch.device("cuda")
    print("GPU is available")
else:
    device = torch.device("cpu")
    print("GPU not available, CPU used")


# In[ ]:


def test_metrics(model, data):
    
    val_X = tensor_generator(data)
    
    model.eval()
    

    out = model(val_X)
    predicted_props = out[:,0].detach().numpy()
    predicted_labels = 1*(predicted_props >0.5)
    
    
    result_tab = pd.DataFrame(data = {'participant': list(data.groupby('participant_id').first().index),
                                     'was_preterm':list(predicted_labels),
                                     'probability':list(predicted_props)})
    

    
    return result_tab


# ## Loading models

# In[ ]:


best_model = load('/usr/local/bin/test_submission.save')


# ## Results

# In[ ]:


test_metrics(best_model, mydata).to_csv('/output/predictions.csv',index = False)


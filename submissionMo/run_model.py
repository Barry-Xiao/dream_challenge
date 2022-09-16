import pandas as pd
import numpy as np
import random

import torch
import torch.nn as nn
import sklearn as sk

from sklearn import metrics
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import GridSearchCV

from collections import Counter,defaultdict, OrderedDict
from itertools import islice
from joblib import dump, load


class Model_Mtd(nn.Module):
    def __init__(self, input_size, output_size, hidden_dim, seq_len, n_layers, fc_size, dropoutrate):
        super(Model_Mtd, self).__init__()

        # Defining some parameters
        self.input_size  = input_size      # number of input node
        self.output_size = output_size     # number of output node
        self.seq_len     = seq_len         # seq_len: number of timepoints (collection period)
        self.fc_size     = fc_size         # size of the fully connected net
        self.n_layers    = n_layers        # number of LSTM/RNN layers
        self.hidden_dim  = hidden_dim      # hidden size of LSTM/RNN, also the size of fully connected NN 1
        
        self.gru = nn.GRU(input_size=input_size, hidden_size=hidden_dim, num_layers=n_layers, batch_first=True)
        self.fc_1 = nn.Linear(in_features=hidden_dim*seq_len, out_features=fc_size[0], bias=False)
        self.fc_2 = nn.Linear(in_features=fc_size[0], out_features=output_size, bias=False)

        # define dropout proportion to prevent overfitting
        self.dropout = nn.Dropout(dropoutrate)
        self.tanh = nn.Tanh()
        
    def forward(self, x, device):
        
        # Initializing hidden state for first input using method defined below
        batch_size = x.size(0)
        h0 = self.init_hidden(batch_size, device)
        #------------ RNN  ------------#
        # outp, hidden = self.rnn(x, h0)
        #------------ LSTM ------------#
        # c0 = self.init_hidden(batch_size, device)
        # outp, hidden = self.lstm(x, (h0, c0))
        #------------ GRU  ------------#
        outp, hidden = self.gru(x, h0)
            
        outp = outp.reshape(outp.shape[0], -1)  # reshaping the data for Dense layer next

        outp = self.fc_1(outp)
        outp = self.tanh(outp)   # relu
        outp = self.dropout(outp)# dropout
        outp = self.fc_2(outp)
        
        return outp, hidden
    
    def init_hidden(self, batch_size, device):
        # This method generates the first hidden state of zeros which we'll use in the forward pass
        hidden = torch.zeros(self.n_layers, batch_size, self.hidden_dim).to(device)
        # We'll send the tensor holding the hidden state to the device we specified earlier as well
        return hidden
    
class Model_pty(nn.Module):
    def __init__(self, input_size, output_size, hidden_dim, seq_len, n_layers, fc_size, dropoutrate):
        super(Model_pty, self).__init__()

        # Defining some parameters
        self.input_size  = input_size      # number of input node
        self.output_size = output_size     # number of output node
        self.seq_len     = seq_len         # seq_len: number of timepoints (collection period)
        self.fc_size     = fc_size         # size of the fully connected net
        self.n_layers    = n_layers        # number of LSTM/RNN layers
        self.hidden_dim  = hidden_dim      # hidden size of LSTM/RNN, also the size of fully connected NN 1
        
        self.gru = nn.GRU(input_size=input_size, hidden_size=hidden_dim, num_layers=n_layers, batch_first=True)
        self.fc_1 = nn.Linear(in_features=hidden_dim*seq_len, out_features=fc_size[0], bias=False)
        self.fc_2 = nn.Linear(in_features=fc_size[0], out_features=fc_size[1], bias=False)
        self.fc_3 = nn.Linear(in_features=fc_size[1], out_features=fc_size[2], bias=False)
        self.fc_4 = nn.Linear(in_features=fc_size[2], out_features=output_size, bias=False)
        # self.relu = nn.ReLU()
        self.tanh = nn.Tanh()
        # define dropout proportion to prevent overfitting
        self.dropout = nn.Dropout(dropoutrate)

    
    def forward(self, x, device):
        
        # Initializing hidden state for first input using method defined below
        batch_size = x.size(0)
        h0 = self.init_hidden(batch_size, device)
        
        #------------ RNN  ------------#
        # outp, hidden = self.rnn(x, h0)
        #------------ LSTM ------------#
        # c0 = self.init_hidden(batch_size, device)
        # outp, hidden = self.lstm(x, (h0, c0))
        #------------ GRU  ------------#
        outp, hidden = self.gru(x, h0)
        
        outp = outp.reshape(outp.shape[0], -1)  # reshaping the data for Dense layer next
        
        outp = self.tanh(outp)   # relu
        outp = self.dropout(outp)# dropout
        outp = self.fc_1(outp)   # first Dense
        outp = self.tanh(outp)   # relu
        outp = self.dropout(outp)# dropout
        outp = self.fc_2(outp)   # 2nd Dense
        outp = self.tanh(outp)   # relu
        outp = self.dropout(outp)# dropout
        outp = self.fc_3(outp)   # 3rd Output
        outp = self.tanh(outp)   # relu
        outp = self.dropout(outp)# dropout
        outp = self.fc_4(outp)   # 4th Ouuput
        
        return outp, hidden
    
    def init_hidden(self, batch_size, device):
        # This method generates the first hidden state of zeros which we'll use in the forward pass
        hidden = torch.zeros(self.n_layers, batch_size, self.hidden_dim).to(device)
        # We'll send the tensor holding the hidden state to the device we specified earlier as well
        return hidden

    
class Model_txy(nn.Module):
    def __init__(self, input_size, output_size, hidden_dim, seq_len, n_layers, fc_size, dropoutrate):
        super(Model_txy, self).__init__()

        # Defining some parameters
        self.input_size  = input_size      # number of input node
        self.output_size = output_size     # number of output node
        self.seq_len     = seq_len         # seq_len: number of timepoints (collection period)
        self.fc_size     = fc_size         # size of the fully connected net
        self.n_layers    = n_layers        # number of LSTM/RNN layers
        self.hidden_dim  = hidden_dim      # hidden size of LSTM/RNN, also the size of fully connected NN 1
        
        self.gru = nn.GRU(input_size=input_size, hidden_size=hidden_dim, num_layers=n_layers, batch_first=True)
        self.fc_1 = nn.Linear(in_features=hidden_dim*seq_len, out_features=fc_size[0], bias=False)
        self.fc_2 = nn.Linear(in_features=fc_size[0], out_features=fc_size[1], bias=False)
        self.fc_3 = nn.Linear(in_features=fc_size[1], out_features=fc_size[2], bias=False)
        self.fc_4 = nn.Linear(in_features=fc_size[2], out_features=output_size, bias=False)
        # self.relu = nn.ReLU()
        self.tanh = nn.Tanh()
        # define dropout proportion to prevent overfitting
        self.dropout = nn.Dropout(dropoutrate)

    
    def forward(self, x, device):
        
        # Initializing hidden state for first input using method defined below
        batch_size = x.size(0)
        h0 = self.init_hidden(batch_size, device)
        #------------ RNN  ------------#
        # outp, hidden = self.rnn(x, h0)
        #------------ LSTM ------------#
        # c0 = self.init_hidden(batch_size, device)
        # outp, hidden = self.lstm(x, (h0, c0))
        #------------ GRU  ------------#
        outp, hidden = self.gru(x, h0)
        
        outp = outp.reshape(outp.shape[0], -1)  # reshaping the data for Dense layer next
        
        outp = self.tanh(outp)   # relu
        outp = self.dropout(outp)# dropout
        outp = self.fc_1(outp)   # first Dense
        outp = self.tanh(outp)   # relu
        outp = self.dropout(outp)# dropout
        outp = self.fc_2(outp)   # 2nd Dense
        outp = self.tanh(outp)   # relu
        outp = self.dropout(outp)# dropout
        outp = self.fc_3(outp)   # 3rd Output
        outp = self.tanh(outp)   # relu
        outp = self.dropout(outp)# dropout
        outp = self.fc_4(outp)   # 4th Ouuput
        outp = self.tanh(outp)   # relu
        
        return outp, hidden
    
    def init_hidden(self, batch_size, device):
        # This method generates the first hidden state of zeros which we'll use in the forward pass
        hidden = torch.zeros(self.n_layers, batch_size, self.hidden_dim).to(device)
        # We'll send the tensor holding the hidden state to the device we specified earlier as well
        return hidden
    
    
def Data_Reshaper_Input(data, seq_length):
    
    numsubjects = len(np.unique(data['participant_id']))
    myvary = list(data.columns.values)[2:data.shape[1]]
    num_covariates = len(myvary)
    
    myinput = np.zeros((numsubjects, seq_length, num_covariates), dtype=np.float32)
    for i in range(num_covariates):
        data_wide = data.pivot_table(index=['participant_id'], columns='collect_period', values=myvary[i])
        data_wide = data_wide.sort_index(axis=1)
        data_wide = data_wide.fillna(0)
        tmpindex = data_wide._get_numeric_data().columns.values - 1
        tmpindex = tmpindex.astype(int)
        # time varying variables need to impute all and no records are denoted as 0
        for j in range(numsubjects):
                myinput[j,tmpindex,i] = data_wide.iloc[[j]]
    return myinput

    
    
def InputLoader(data_dir, feature_dir, meta_data, finalperiod):
    
    participant_id = meta_data['participant_id']
    collect_period = meta_data['collect_period']
   
    Input_data = pd.DataFrame(pd.read_csv(data_dir, delimiter=','))
    selectedfeature = pd.DataFrame(pd.read_csv(feature_dir, delimiter=','))
    Input_data = Input_data.iloc[:,selectedfeature['id']]
    Input_data = pd.concat([participant_id, collect_period, Input_data], axis=1)
        
    # Average within each collection period
    Input_data = Input_data.groupby(['participant_id', 'collect_period'], as_index = False).mean()
    
    #---- Input features reshaper ----#
    mydata_input = Data_Reshaper_Input(data=Input_data, seq_length=finalperiod)
    
    return mydata_input

def evaluate(model, device, myinput, finalperiod, cutoff=0.5):
    
    model.eval()
    
    # predicted labels
    myinput  = torch.from_numpy(myinput).float().to(device)
    myoutput_nn, hidden = model(myinput, device)
    myoutput_nn = myoutput_nn.reshape(([myinput.shape[0], finalperiod, 2]))
    output_prob = nn.functional.softmax(myoutput_nn, dim=2)
    mypredprob = output_prob[:,finalperiod-1,:].cpu().detach().numpy()
    mypred = 1*(mypredprob[:,0] > cutoff)
    
    return mypred, mypredprob





# task = "was_preterm"
finalperiod = 5
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')


#-------------------------------------------#
#---- Metadata                          ----#
#-------------------------------------------#

meta_dir      = '/input/metadata/metadata.csv'
alpha_dir     = '/input/alpha_diversity/alpha_diversity.csv'
cst_dir       = '/input/community_state_types/cst_valencia.csv'

meta_data = pd.read_csv(meta_dir, delimiter=',')
meta_data.replace('Unknown', np.nan, inplace=True)
meta_data = meta_data[['participant_id', 'collect_wk', 'age', 'race']]
    

alpha_data = pd.DataFrame(pd.read_csv(alpha_dir, delimiter=','))
cst_data = pd.DataFrame(pd.read_csv(cst_dir, delimiter=','))

meta_data = pd.concat([meta_data, alpha_data['shannon'], alpha_data['inv_simpson'], alpha_data['rooted_pd'], cst_data['CST']], axis=1)

for i in range(1,meta_data.shape[1]):
        if meta_data.iloc[:,i].dtypes == object:
            meta_data.iloc[:,i] = meta_data.iloc[:,i].astype('category').cat.codes + 1
            meta_data.iloc[:,i] = meta_data.iloc[:,i].astype('float64')
            
# create new variable 'collect_period'
meta_data['collect_period'] = 1
meta_data.loc[(meta_data['collect_wk']>=8)  & (meta_data['collect_wk']<=14),'collect_period'] = 2
meta_data.loc[(meta_data['collect_wk']>=15) & (meta_data['collect_wk']<=21),'collect_period'] = 3
meta_data.loc[(meta_data['collect_wk']>=22) & (meta_data['collect_wk']<=28),'collect_period'] = 4
meta_data.loc[(meta_data['collect_wk']>=29) & (meta_data['collect_wk']<=32),'collect_period'] = 5
meta_data.loc[(meta_data['collect_wk']>=33), 'collect_period']                                = 6

# remove when submit
meta_data = meta_data[meta_data['collect_period']<=finalperiod]

meta_data = meta_data.groupby(['participant_id', 'collect_period'], as_index = False).mean()


# scale the input features in this data set
columns = ['collect_wk', 'age', 'race', 'shannon', 'inv_simpson', 'rooted_pd', 'CST']
for col in columns:
    meta_data[col] = MinMaxScaler().fit_transform(np.array(meta_data[col]).reshape(-1,1))
    
participant_id = meta_data['participant_id']
collect_period = meta_data['collect_period']

meta_data_Input = Data_Reshaper_Input(data=meta_data, seq_length=finalperiod)

# load trained model
model_Mtd = load('/usr/local/bin/Mtd_waspreterm.save')

Mtd_pred, Mtd_prob = evaluate(model_Mtd, device, meta_data_Input, finalperiod, cutoff=0.5)


#-------------------------------------------#
#---- ptydata                          ----#
#-------------------------------------------#


pty_dir = '/input/phylotypes/phylotype_relabd.5e_1.csv'
pty_feature = './selected_feature/ptyfeature_was_preterm_dot5.csv'

pty_data_Input = InputLoader(pty_dir, pty_feature, meta_data, finalperiod)

# load trained model
model_pty = load('/usr/local/bin/pty_waspreterm.save')

pty_pred, pty_prob = evaluate(model_pty, device, pty_data_Input, finalperiod, cutoff=0.5)


#-------------------------------------------#
#---- txydata                           ----#
#-------------------------------------------#


txy_dir = '/input/taxonomy/taxonomy_relabd.genus.csv'
txy_feature_dir = './selected_feature/txyfeature_was_preterm_gen.csv'

txy_data_Input = InputLoader(txy_dir, txy_feature_dir, meta_data, finalperiod)

# load trained model
model_txy = load('/usr/local/bin/txy_waspreterm.save')

txy_pred, txy_prob = evaluate(model_txy, device, txy_data_Input, finalperiod, cutoff=0.5)


#-------------------------------------------#
#---- 2nd stage                         ----#
#-------------------------------------------#


# load('/usr/local/bin/S2Logistic.save')
L2Logistic_model = load('/usr/local/bin/L2logistic_waspreterm.save')

x_test = np.array(np.column_stack([Mtd_prob, pty_prob, txy_prob])).reshape(-1, 3*2)

final_prob = L2Logistic_model.predict_proba(x_test)[:,1]
final_pred = L2Logistic_model.predict(x_test)


result_tab = pd.DataFrame(data = {'participant': list(meta_data.groupby('participant_id').first().index),
                                  'was_preterm': list(final_pred),
                                  'probability': list(final_prob)})

result_tab.to_csv('/output/predictions.csv', index = False)

import pandas as pd
import numpy as np

import torch
import torch.nn as nn

from sklearn import metrics
from sklearn import preprocessing

import matplotlib.pyplot as plt

import random

import sklearn as sk

from itertools import islice

from collections import Counter,defaultdict, OrderedDict




# data directory
alpha_dir     = 'D:/DREAM/alpha_diversity/alpha_diversity.csv'
cst_dir       = 'D:/DREAM/community_state_types/cst_valencia.csv'
meta_dir      = 'D:/DREAM/preprocessed/metadata_imputed.csv'
# krdlong_dir   = '/Users/mli171/Desktop/JHU/3Summer2022_JHU/DREAM/training_data_2022-05-27/pairwise_distance/krd_distance_long.csv'
# krdwide_dir   = '/Users/mli171/Desktop/JHU/3Summer2022_JHU/DREAM/training_data_2022-05-27/pairwise_distance/krd_distance_wide.csv'
phylotype_dir = 'D:/DREAM/phylotypes/phylotype_relabd.1e0.csv'
taxonomy_dir  = 'D:/DREAM/taxonomy/taxonomy_relabd.family.csv'




meta_data = pd.DataFrame(pd.read_csv(meta_dir, delimiter=','))
meta_data = meta_data[['participant_id', 'project', 'delivery_wk', 'collect_wk', 'age_imp', 'race_imp']]


alpha_data = pd.DataFrame(pd.read_csv(alpha_dir, delimiter=','))
cst_data = pd.DataFrame(pd.read_csv(cst_dir, delimiter=','))
meta_data = pd.concat([meta_data, alpha_data['shannon'], alpha_data['inv_simpson'], alpha_data['rooted_pd'], cst_data['CST']], axis=1)


print(meta_data.shape)

for i in range(1,meta_data.shape[1]):
    if meta_data.iloc[:,i].dtypes == object:
        meta_data.iloc[:,i] = meta_data.iloc[:,i].astype('category').cat.codes + 1
        meta_data.iloc[:,i] = meta_data.iloc[:,i].astype('float64')

# create new variable collection period
meta_data['collect_period'] = 1
meta_data.loc[(meta_data['collect_wk']>=8)  & (meta_data['collect_wk']<=14),'collect_period'] = 2
meta_data.loc[(meta_data['collect_wk']>=15) & (meta_data['collect_wk']<=21),'collect_period'] = 3
meta_data.loc[(meta_data['collect_wk']>=22) & (meta_data['collect_wk']<=28),'collect_period'] = 4
meta_data.loc[(meta_data['collect_wk']>=29) & (meta_data['collect_wk']<=32),'collect_period'] = 5
meta_data.loc[(meta_data['collect_wk']>=33), 'collect_period']                                = 6

collect_period = meta_data['collect_period']
participant_id = meta_data['participant_id']

# create class label
meta_data['classlabel'] = 1*(meta_data['delivery_wk'] < 32)
#meta_data['was_preterm'] = 1*(meta_data['delivery_wk'] < 37)
#meta_data['was_early_preterm'] = 1*(meta_data['delivery_wk'] < 32)


# Filtered out observations with "collect_wk<=32" == "collect_period<=5" 
meta_data = meta_data[meta_data['collect_period']<=4]
# Average within each collection period
#meta_data = meta_data.groupby(['participant_id', 'collect_period'], as_index = False).mean()
print(meta_data.shape)



taxonomy_data = pd.DataFrame(pd.read_csv(taxonomy_dir, delimiter=','))
taxonomy_data = pd.concat([participant_id, collect_period, taxonomy_data], axis=1)

# Filtered out observations with "collect_wk<=32" == "collect_period<=5" 
taxonomy_data = taxonomy_data[taxonomy_data['collect_period']<=4]
# Average within each collection period
#taxonomy_data = taxonomy_data.groupby(['participant_id', 'collect_period'], as_index = False).mean()
print(taxonomy_data.shape)


# delete the temporary filter-used and ID columns
taxonomy_data = taxonomy_data.drop(["participant_id", "collect_period"], axis = 1)


phylotype_data = pd.DataFrame(pd.read_csv(phylotype_dir, delimiter=','))
phylotype_data = pd.concat([participant_id, collect_period, phylotype_data], axis=1)

# Filtered out observations with "collect_wk<=32" == "collect_period<=5" 
phylotype_data = phylotype_data[phylotype_data['collect_period']<=4]
# Average within each collection period
#phylotype_data = phylotype_data.groupby(['participant_id', 'collect_period'], as_index = False).mean()
print(phylotype_data.shape)


# delete the temporary filter-used and ID columns
phylotype_data = phylotype_data.drop(["participant_id", "collect_period"], axis = 1)





mydata = pd.concat([meta_data, phylotype_data, taxonomy_data], axis=1)



mydata






mydata_input = mydata.drop('classlabel', axis=1)
mydata_output = mydata[['participant_id', 'collect_period', 'classlabel']]
mydata_project = mydata[['participant_id', 'collect_period', 'project']]




def Data_Reshaper_Input(data, seq_length):
    num_samples = len(np.unique(data['participant_id']))
    myvary = list(data.columns.values)
    myvary.remove('participant_id')
    myvary.remove('collect_period')
    myvary.remove('project')
    myvary.remove('delivery_wk')
    myvary.remove('collect_wk')
    num_covariates = len(myvary)
    
    myinput = np.zeros((num_samples, num_covariates, seq_length), dtype=np.float32)
    
    for i in range(num_covariates):
        data_wide = data.pivot_table(index=['participant_id'], columns='collect_period', values=myvary[i], aggfunc='mean')
        data_wide = data_wide.sort_index(axis=1)
        data_wide = data_wide.fillna(0)
        tmpindex = data_wide._get_numeric_data().columns.values - 1
        tmpindex = tmpindex.tolist()
        # time varying variables need to impute all and no records are denoted as 0
        for j in range(num_samples):
                myinput[j,i,tmpindex] = data_wide.iloc[[j]]
    return myinput




mydata_input = Data_Reshaper_Input(mydata_input, seq_length=4)
print(mydata_input.shape)




# =============================================================================
# def Data_Reshaper_Output_ManytoMany(data, seq_length, num_covariates):
#     
#     num_samples = len(np.unique(data['participant_id']))
# 
#     data_wide = data.pivot_table(index=['participant_id'], columns='collect_period', values="classlabel")
#     data_wide = data_wide.sort_index(axis=1)
# 
#     myoutput = np.zeros((num_samples, num_covariates, seq_length), dtype=np.float32)
# 
#     myoutput[:,0,data_wide.columns.values-1] = data_wide
#     myoutput[:,1,data_wide.columns.values-1] = 1 - data_wide
#     myoutput[np.isnan(myoutput)] = 0
#     
#     return myoutput
# =============================================================================


def Data_Reshaper_Output(data, vals):
    
    data_wide = data.pivot_table(index=['participant_id'], columns='collect_period', values=vals)
    data_wide = data_wide.sort_index(axis=1)
    
    data_wide = data_wide.apply(lambda row: row.fillna(row.mean()), axis=1)

    myoutput = data_wide.iloc[:,1]
    
    return myoutput






mydata_output = Data_Reshaper_Output(mydata_output, "classlabel")
print(mydata_output.shape)

mydata_project = Data_Reshaper_Output(mydata_project, "project")
print(mydata_project.shape)







uniquenames, counts = np.unique(meta_data["participant_id"], return_counts=True)
subjects = list(uniquenames)
seq_max_len = max(counts)

print("# of subjects = ", len(subjects))
print("# of samples  = ", meta_data.shape[0])
print("# of taxonnomy features = ", len(list(taxonomy_data)))
print("# of phylotype features = ", len(list(phylotype_data)))


    
# set myseed=None to have complete random state
myseed = 0
random.seed(myseed)
ID_shuffle = random.sample(range(0,len(subjects)),len(subjects))
prop = [0.6,0.3,0.1]
splitID_train = ID_shuffle[0:len(subjects)] 
splitID_valid = ID_shuffle[(int(len(subjects)*(prop[0]+prop[1]))+2):len(subjects)]
splitID_test = ID_shuffle[(int(len(subjects)*(prop[0]+prop[1]))+2):len(subjects)]




# apply to each data sets
mytrain_input = mydata_input[splitID_train,:,:]
myvalid_input = mydata_input[splitID_valid,:,:]
mytest_input  = mydata_input[splitID_test,:,:]
mytrain_output = mydata_output[splitID_train]
myvalid_output = mydata_output[splitID_valid]
mytest_output  = mydata_output[splitID_test]

print(mytrain_input.shape)
print(myvalid_input.shape)
print(mytest_input.shape)
print(mytrain_output.shape)
print(myvalid_output.shape)
print(mytest_output.shape)



from tsai.all import *
from tsai.inference import load_learner
from IPython.display import clear_output


    
    
X, y, splits = combine_split_data([mytrain_input, myvalid_input], [mytrain_output, myvalid_output])
X2, y2, splits2 = combine_split_data([mytrain_input, mytest_input], [mytrain_output, mytest_output])

    
    #X[1,1,1] = np.nan
    
# =============================================================================
#     tfms  = [None, [Categorize()]]
#     dsets = TSDatasets(X, y, tfms=tfms, splits=splits, inplace=True)
#     dls   = TSDataLoaders.from_dsets(dsets.train,
#                                      dsets.valid, 
#                                      bs=[64, 128],
#                                      batch_tfms=[TSStandardize()], 
#                                      num_workers=0)
#     model = TST(dls.vars, dls.c, dls.len, dropout=0.3, fc_dropout=0.9)
#     learn = Learner(dls, model, metrics=accuracy)
#     learn.fit_one_cycle(30, lr_max=1e-3)
#     # dls = learn.dls
#     valid_dl = dls.valid
#     valid_probas, valid_targets, valid_preds = learn.get_preds(dl=valid_dl, with_decoded=True)
#     valid_probas, valid_targets, valid_preds
#     (valid_targets == valid_preds).float().mean()
# =============================================================================
    

    
bs = 64
tfms  = [None, [Categorize()]]
dsets = TSDatasets(X, y, tfms=tfms, splits=splits)
dls   = TSDataLoaders.from_dsets(dsets.train, dsets.valid, bs=[bs, bs*1])
valid_dl = dls.valid

dsets2 = TSDatasets(X2, y2, tfms=tfms, splits=splits2)
dls2   = TSDataLoaders.from_dsets(dsets2.train, dsets2.valid, bs=[bs, bs*1])
test_dl = dls2.valid


results = pd.DataFrame(columns=['arch', 'hyperparams', 'train loss', 'valid loss', 'Precision', 'Recall', 'AUC', 'accuracy', 'time'])


archs = [(LSTMPlus, {'n_layers':1, 'bidirectional':True, 'rnn_dropout':0.5, 'fc_dropout':0.5, 'feature_extractor':MultiConv1d(2377, kss=[1,3,5,7])})]


for i, (arch, k) in enumerate(archs):
    model = create_model(arch, dls=dls, **k)
    print(model.__class__.__name__)
    learn = Learner(dls, model, metrics=[Precision(),Recall(),RocAucBinary(),accuracy], loss_func = CrossEntropyLossFlat())
    start = time.time()
    learn.fit_one_cycle(15, 1e-3)
    elapsed = time.time() - start
    learn.plot_metrics(suptitle = arch.__name__)
    vals = learn.recorder.values[-1]
    results.loc[i] = [arch.__name__, k, vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], int(elapsed)]
    #valid_probas, valid_targets, valid_preds = learn.get_preds(dl=valid_dl, with_decoded=True, save_preds=None, save_targs=None)
    test_probas, test_targets, test_preds = learn.get_preds(dl=test_dl, with_decoded=True, save_preds=None, save_targs=None)
    results.iloc[i,[4,5,6,7]] = learn.recorder.values[0][1,2,3,4]
    clear_output()
    #display(results)


learn.save_all(path='D:/DREAM/models_earlypreterm', dls_fname='dls_earlypreterm', model_fname='model_earlypreterm', learner_fname='learner_earlypreterm')
learn = load_learner('D:/DREAM/submissionJiuyao/earlypreterm/learner_earlypreterm.pkl')



print(results.iloc[:,[0,2,3,4,5,6,7,8]].to_string())

print(results.iloc[:,[0,6]].to_string())


import pickle
pickle.dump(results, open('earlypreterm_15.pkl', 'wb'))

results_new = pickle.load(open('earlypreterm_15.pkl', 'rb'))
print(results_new.iloc[:,[0,1,6]].to_string())

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 11:51:09 2018

@author: mayasrikanth
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 10:57:26 2018

@author: mayasrikanth
"""
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import RepeatedKFold

def Ecv(y_val_pred, y_train_val):
    Ecv = 0
    for counter, element in enumerate(y_val_pred):
        if(element != y_train_val[counter]):
            Ecv += 1
    Ecv = Ecv / len(y_val_pred)
    return Ecv


def soft_margin_SVM(all_C):
    train_set = np.loadtxt(fname="features.train.txt")  
    test_set = np.loadtxt(fname="features.test.txt")
    
    # Create X_train
    X_train = train_set[:,1:]
    # Create y_train
    y_train = train_set[:,0]
    y_train[y_train == 1] = 1
    y_train[y_train == 5] = -1
    
    y_train_temp = []
    X_train_temp = []
      # Discard all other digits from X_train and y_train
    for counter, e in  enumerate(y_train):
        if(e == -1 or e == 1):
            y_train_temp.append(e)
            X_train_temp.append(X_train[counter])
    X_train_temp = np.array(X_train_temp)
    X_train = X_train_temp
    y_train_temp = np.array(y_train_temp)
    y_train = y_train_temp
    
    # Ecv = 0
    win_frequency = [0, 0, 0, 0, 0]
    rkf = RepeatedKFold(n_splits=10, n_repeats=100)
    c = 0
    for train_index, test_index in rkf.split(X_train):
        X_train_red, X_train_val = X_train[train_index], X_train[test_index]
        y_train_red, y_train_val = y_train[train_index], y_train[test_index]
        all_Ecvs = [0, 0, 0, 0, 0]
        for counter, model in enumerate(all_C):
            # fit CLF to X_train_red
            clf = SVC(gamma=1,coef0=1, C=model, kernel='poly', degree=2)
            X_train_red = np.array(X_train_red)
            y_train_red = np.array(y_train_red)
            clf.fit(X_train_red, y_train_red)
            # Use clf to predict on validation set and evaluate Ecv
            # Ecv_run = 0
            y_val_pred = clf.predict(X_train_val)
            all_Ecvs[counter] = Ecv(y_val_pred, y_train_val)
     
        win_index = all_Ecvs.index(min(all_Ecvs))
        c += 1
        print("Choosing the model: " + str(c))
        win_frequency[win_index] += 1
        
            
        # Ecv = Ecv + Ecv_run
    
    for counter, element in enumerate(win_frequency):
        chosen = element / 1000
        print("Model: " + str(all_C[counter]) + " chosen: " + str(chosen))
    # Ecv = Ecv/100
   # print("Model: " + str(choose_C))
    #print("Ecv averaged over 100 runs: " + str(Ecv))

def crossValidation():
    models = [0.0001, 0.001, 0.01, 0.1, 1]
    soft_margin_SVM(models)
       
        

    

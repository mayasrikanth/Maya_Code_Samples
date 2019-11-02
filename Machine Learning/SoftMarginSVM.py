#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 25 10:57:26 2018

@author: mayasrikanth
"""
import numpy as np
from sklearn.svm import SVC


def soft_margin_SVM(choose_C, choose_degree):
    train_set = np.loadtxt(fname="features.train.txt")  
    test_set = np.loadtxt(fname="features.test.txt")
    
    # Create X_train
    X_train = train_set[:,1:]
    # Create y_train
    y_train = train_set[:,0]
    y_train[y_train == 1] = 1
    y_train[y_train == 5] = -1
    
    ones = y_train[y_train == 1]
    fives = y_train[y_train == -1]
s    print("Ones: " + str(len(ones)))
    print("Fives: " + str(len(fives)))
    
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
    print("Length of training set: " + str(len(X_train)))
    

    
    # Creating SVC, poly kernel w/ Q = 2, C = 0.01
    clf = SVC(gamma=1, coef0=1, C=choose_C, kernel='poly', degree=choose_degree, decision_function_shape='ovr')
    X_train = np.array(X_train)
    y_train = np.array(y_train)
    clf.fit(X_train, y_train)
    
    # Create X_test
    X_test = test_set[:, 1:]
    # Create y_test
    y_test = test_set[:, 0]
    y_test[y_test == 1] = 1
    y_test[y_test == 5] = -1
    
    y_test_temp = []
    X_test_temp = []
    # Discard all other digits from X_train and y_train
    for counter, e in  enumerate(y_test):
        if(e == -1 or e == 1):
            y_test_temp.append(e)
            X_test_temp.append(X_test[counter])
    
    X_test_temp = np.array(X_test_temp)
    X_test = X_test_temp
    y_test_temp = np.array(y_test_temp)
    y_test = y_test_temp
    
    
    # Use clf to predict on test set 
    y_test_pred = clf.predict(X_test)
    
    # Evaluate Eout
    Eout = 0
    for counter, element in enumerate(y_test_pred):
        if(element != y_test[counter]):
            Eout += 1
    Eout = Eout / len(y_test_pred)
            
    
    # Use clf to predict on training set 
    y_train_pred = clf.predict(X_train)
    
    # Evaluate Ein
    Ein = 0
    for counter, element in enumerate(y_train_pred):
        if(element != y_train[counter]):
            Ein += 1
    Ein = Ein / len(y_train_pred)
    
    # print("Classifer: " + str(classifier) + " vs. all ")
    print("C: " + str(choose_C))
    print("Eout: " + str(Eout))
    print("Ein: " + str(Ein))
    support_v = np.sum(clf.n_support_)
    print("Amount of support vectors: " + str(support_v))
    
    
    
    

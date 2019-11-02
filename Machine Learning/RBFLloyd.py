#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 00:19:22 2018

@author: mayasrikanth
"""

import numpy as np
from sklearn.svm import SVC
from sklearn.cluster import KMeans


def generateTargets(X_train):
    y_train = []
    for e in X_train:
        x1 = e[1]
        x2 = e[2]
        add = 0.25 * np.sin(np.pi * x1)
        y_train.append(np.sign(x2 - x1 + add))
    y_train = np.array(y_train)
    return y_train

def generateInput(N):
    X_train = []
    for i in range(0, 100):
        x1 = np.random.uniform(-1, 1)
        x2 = np.random.uniform(-1, 1)
        X_train.append([1, x1, x2])
    
    X_train = np.array(X_train)
    return X_train

def discardRunKMeans(labels):
    for i in range(0,9):
        if(i not in labels):
            return True
    return False

def discardRunSVM(y_train, y_pred):
    # Find Ein for RBF kernel w/ hard-margin SVM
    Ein = 0
    for counter, e in enumerate(y_pred):
        if(e != y_train[counter]):
            Ein += 1
    Ein = Ein / len(y_pred)
    if(Ein != 0):
        return True
    else:
        return False

def phi_term(x, cluster):
    diff = np.subtract(x, cluster)
    norm = np.linalg.norm(diff)
    norm_sq = norm**2
    product = -2 * norm_sq
    phi_term = np.e**(product)
    # print("phi term: " + str(phi_term))
    return phi_term
    

def PseudoInverse(X_train, clusters, y_train):
    N = len(X_train)  # Amount of points
    K = len(clusters) # Amount of clusters
    Phi = [] 
    # Initialize Phi
    for row in range(0, N):
        temp_array = []
        for col in range(0,K + 1):
            if(col == 0):
                temp_array.append(1)
            else:
                temp_array.append(0)
        Phi.append(temp_array)
    # Populate Phi
    for row in range(0,N):
        x = X_train[row]
        for col in range(1, K + 1):
            cluster = clusters[col -1]
            Phi[row][col] = phi_term(x, cluster)
    Phi = np.matrix(Phi)     
    # print("The matrix phi: " + str(Phi))
    inverse = np.linalg.pinv(Phi)
    #y_transpose = np.matrix(np.transpose(y_train))
    y_train = np.transpose(np.matrix(y_train))
    
    weights = inverse * y_train
    weights = weights.tolist()

    return weights

def KMeans_predict(x, weights, clusters):
    h_x = 0
    # Adding bias term
    h_x = h_x + weights[0][0]
    for counter, cluster in enumerate(clusters,1):
        diff = np.subtract(x, cluster)
        norm = np.linalg.norm(diff)
        norm_sq = norm**2
        product = -2 * norm_sq
        exponent_term = np.e**(product)

        add = weights[counter][0] * exponent_term
        h_x = h_x + add
   
    return h_x

                         
# Hard-margin SVM and Kmeans clustering
def SVM_Kmeans():
    X_train = generateInput(100)
    y_train = generateTargets(X_train)
    
    # Initialize SVC for RBF
    # RBF Kernel, gamma = 1.5 
    clf = SVC(gamma=1.5, C=np.inf, kernel='rbf')
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_train)

    # Initialize kmeans machine
    kmeans = KMeans(n_clusters=9).fit(X_train)
    
    # Indexes into the array of clusters 
    labels = kmeans.labels_
    # All clusters 
    clusters = kmeans.cluster_centers_
    
    # If cluster becomes empty, discard run and repeat
    while(discardRunKMeans(labels) == True or discardRunSVM(y_train, y_pred)):
        X_train = generateInput(100)
        y_train = generateTargets(X_train)
        # Initialize KMeans 
        kmeans = KMeans(n_clusters=9, n_init=1, n_jobs=1).fit(X_train)
        labels = kmeans.labels_
        clusters = kmeans.cluster_centers_
        
    weights_Kmeans = PseudoInverse(X_train, clusters, y_train)
    
    
    # Find Ein for Kmeans 
    Ein_Kmeans = 0
    for counter, x in enumerate(X_train): 
        if(np.sign(KMeans_predict(x, weights_Kmeans, clusters)) != y_train[counter]):
            Ein_Kmeans += 1
    Ein_Kmeans = Ein_Kmeans/len(X_train)
    # print("Kmeans_Ein: " + str(Ein_Kmeans))
     

      
    # Eout test set
    X_test = generateInput(100)
    y_test = generateTargets(X_test)
    
    # Eout for RBF kernel w/ hard-margin SVM
    y_test_pred = clf.predict(X_test)
    SVM_Eout = 0
    for counter, e in enumerate(y_test_pred):
        if(e != y_test[counter]):
            SVM_Eout += 1
    SVM_Eout = SVM_Eout / len(y_test_pred)
    
    # Eout for Kmeans 
    Kmeans_Eout = 0
    for counter, x in enumerate(X_test):
        if(np.sign(KMeans_predict(x, weights_Kmeans, clusters)) != y_test[counter]):
            Kmeans_Eout += 1
    Kmeans_Eout = Kmeans_Eout/len(X_test)
    
    #print("SVM Eout: " + str(SVM_Eout))
    #print("Kmeans Eout: " + str(Kmeans_Eout))
    return(Kmeans_Eout, SVM_Eout, Ein_Kmeans)
    
    
def compare():
    SVM_win = 0
    Kmeans_win = 0 
    # SVM_Eout_avg = 0
    Kmeans_Ein_avg = 0
    Kmeans_Eout_avg = 0
    count_Ein_Kmeans = 0
    for i in range(0, 1000):
        (Kmeans_Eout, SVM_Eout, Ein_Kmeans) = SVM_Kmeans()
        if(Ein_Kmeans == 0):
            count_Ein_Kmeans += 1
        Kmeans_Ein_avg += Ein_Kmeans
        Kmeans_Eout_avg += Kmeans_Eout
        if(Kmeans_Eout < SVM_Eout):
            Kmeans_win += 1
        elif(SVM_Eout < Kmeans_Eout):
            SVM_win += 1
            
    count_Ein_Kmeans = count_Ein_Kmeans/1000
    #SVM_Eout_avg = SVM_Eout_avg / 1000
    Kmeans_Eout_avg = Kmeans_Eout_avg / 1000
    Kmeans_Ein_avg = Kmeans_Ein_avg / 1000
    
    
    Kmeans_win_freq = Kmeans_win / 1000
    SVM_win_freq = SVM_win / 1000
    print("Probability that Ein=0 for regular RBF: " + str(count_Ein_Kmeans))
    print("Kmeans wins: " + str(Kmeans_win_freq))
    print("SVM wins: " + str(SVM_win_freq))
    print("Kmeans Ein avg: " + str(Kmeans_Ein_avg))
    print("Kmeans Eout avg: " + str(Kmeans_Eout_avg))


# for K = 9
# Kmeans wins: 0.112
# SVM wins: 0.788
    
# for K = 12
# Kmeans wins: 0.138
# SVM wins: 0.756

        
    


def run():
    inseparable = 0
    Ein = 0
    
    # If RBF_Ein != 0, redo run
    for i in range(0, 1000):
        Ein = SVM_RBF()
        if(Ein != 0):
            inseparable += 1
            
    inseparable = inseparable/1000
    print("Dataset not separable by RBF kernel: " + str(inseparable))
        
        
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import sys
import math
import re

def read_sample(file):
    samples = np.loadtxt(file, dtype=bytes).astype(str)
    samples_ddg = dict()
    for i in samples:
        if float(i[1]) < 1:
            samples_ddg[i[0]] = 0
        else:
            samples_ddg[i[0]] = 1
    return samples_ddg

def read_feature(file):
    feafile = np.loadtxt(file, dtype=bytes).astype(str)
    samples_fea = dict()
    for i in feafile:
        samples_fea[i[0]] = i[1:]
    return samples_fea

if __name__ == '__main__':
    fil = sys.argv[1]
    train_data = fil+'/datasets/MPR233.txt'
    train_bio_fea = fil+'/features/bio_feature_RNA.txt'
    GB2 = fil+"/features/EWC.txt"
   
    
    bio_fea = fil+'/output/bio_feature/bio_feature.txt'
    gb2 = fil+"/output/energy_feature/GB2/EWC.txt"
    
    train = read_sample(train_data)
    ###########################################################################
    samp_bio_fea = read_feature(train_bio_fea)
    it_bio_fea = np.loadtxt(bio_fea, dtype=bytes).astype(str)
    
    train_bio_labels = []
    train_bio_samples = []
    test_bio_samples = []
    
    for site, ddg in train.items():
            train_bio_labels.append(ddg)
            train_bio_samples.append([float(x) for x in samp_bio_fea[site]])
     
    test_bio_samples.append(it_bio_fea[1:])      

    rf = RandomForestClassifier(n_estimators=500)
    rf.fit(train_bio_samples, train_bio_labels)
    prediction_bio = rf.predict_proba(test_bio_samples)
    
    ###########################################################################
    samp_GB_fea = read_feature(GB2)
    it_GB_fea = np.loadtxt(gb2, dtype=bytes).astype(str)
    
    train_GB_labels = []
    train_GB_samples = []
    test_GB_samples = []
    
    for site, ddg in train.items():
            train_GB_labels.append(ddg)
            train_GB_samples.append([float(x) for x in samp_GB_fea[site]])
     
    test_GB_samples.append(it_GB_fea[:])      

    rf = RandomForestClassifier(n_estimators=500)
    rf.fit(train_GB_samples, train_GB_labels)
    prediction_GB = rf.predict_proba(test_GB_samples)
    
    ###########################################################################
    prediction = 0.5* prediction_GB[:,1] + 0.5*prediction_bio[:,1]
     
    with open(fil+"/output/result.txt","a") as fil:
        fil.write("\n"+"Predicted score:"+str('%.3f' % prediction[0]))
        if prediction[0]>=0.320:
            fil.write("\n"+"Significant decrease:YES")
        else:
            fil.write("\n"+"Significant decrease:NO")
     
    
    
    
    
    
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
    train_data = fil+'/datasets/MPD276.txt'
    train_bio_fea = fil+'/features/bio_feature_DNA.txt'
    train_inter_intra = fil+"/features/EPI.txt"
    train_target = fil+"/features/ETOR.txt"
    
    bio_fea = fil+'/output/bio_feature/bio_feature.txt'
    inter_intra = fil+"/output/energy_feature/GB1/EPI.txt";
    target = fil+"/output/energy_feature/GB1/ETOR.txt"
      
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
    samp_inter_fea = read_feature(train_inter_intra)
    it_inter_fea = np.loadtxt(inter_intra, dtype=bytes).astype(str)
    
    train_inter_labels = []
    train_inter_samples = []
    test_inter_samples = []
    
    for site, ddg in train.items():
            train_inter_labels.append(ddg)
            train_inter_samples.append([float(x) for x in samp_inter_fea[site]])
     
    test_inter_samples.append(it_inter_fea[:])      

    rf = RandomForestClassifier(n_estimators=500)
    rf.fit(train_inter_samples, train_inter_labels)
    prediction_inter = rf.predict_proba(test_inter_samples)

    ###########################################################################
    samp_tar_fea = read_feature(train_target)
    it_tar_fea = np.loadtxt(target, dtype=bytes).astype(str)
    
    train_tar_labels = []
    train_tar_samples = []
    test_tar_samples = []
    
    for site, ddg in train.items():
            train_tar_labels.append(ddg)
            train_tar_samples.append([float(x) for x in samp_tar_fea[site]])
     
    test_tar_samples.append(it_tar_fea[:])      

    rf = RandomForestClassifier(n_estimators=500)
    rf.fit(train_tar_samples, train_tar_labels)
    prediction_tar = rf.predict_proba(test_tar_samples)
    ###########################################################################
    en_predition = 0.4*prediction_tar[:,1] + 0.6*prediction_inter[:,1]
    prediction = 0.6*en_predition + 0.4*prediction_bio[:,1]
    
    with open(fil+"/output/result.txt","a") as fil:
        fil.write("\n"+"Predicted score:"+str('%.3f' % prediction[0]))
        if prediction[0]>=0.471:
            fil.write("\n"+"Significant decrease:YES")
        else:
            fil.write("\n"+"Significant decrease:NO")
        


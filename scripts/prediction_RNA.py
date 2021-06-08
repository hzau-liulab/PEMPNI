import numpy as np
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
import sys
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.metrics import mean_squared_error
import math
import re

def read_feature(file):
    feafile = np.loadtxt(file, dtype=bytes).astype(str)
    samples_fea = dict()
    for i in feafile:
        samples_fea[i[0]] = [float(x) for x in i[1:]]
    return samples_fea

if __name__ == '__main__':
    fil = sys.argv[1]
    train_data = fil +'/datasets/MPR233.txt'
    train_bio_fea = fil+'/features/bio_feature_RNA.txt'
    GB2 = fil+"/features/EWC.txt"
      
    bio_fea = fil+'/output/bio_feature/bio_feature.txt'
    gb2 = fil+"/output/energy_feature/GB2/EWC.txt"
   
    train_samp_ddg = np.loadtxt(train_data, dtype=bytes).astype(str)
    
####################################################################################   
    samp_bio_fea = read_feature(train_bio_fea)
    it_bio_fea = np.loadtxt(bio_fea, dtype=bytes).astype(str)

    train_bio_X = []
    train_bio_Y = []
    for k in train_samp_ddg:
        train_bio_Y.append(float(k[1]))
        train_bio_X.append(samp_bio_fea[k[0]])

    test_bio_X = []
    test_bio_X.append(it_bio_fea[1:])
    
    model_line = RandomForestRegressor(n_estimators=500)
    model_line.fit(train_bio_X, train_bio_Y)
    prediction_bio = model_line.predict(test_bio_X)
    
##############################################################################
    samp_GB_fea = read_feature(GB2)
    it_GB_fea = np.loadtxt(gb2, dtype=bytes).astype(str)

    train_GB_X = []
    train_GB_Y = []
    for k in train_samp_ddg:
        train_GB_Y.append(float(k[1]))
        train_GB_X.append(samp_GB_fea[k[0]])

    test_GB_X = []
    test_GB_X.append(it_GB_fea[:])
    
    model_line = RandomForestRegressor(n_estimators=500)
    model_line.fit(train_GB_X, train_GB_Y)
    prediction_GB = model_line.predict(test_GB_X)
    
###############################################################################
    prediction = 0.5*prediction_GB + 0.5*prediction_bio

    with open(fil+"/output/result.txt","a") as fil:
        fil.write("\n"+str("Predicted affinity change:")+str('%.3f' % prediction[0]))
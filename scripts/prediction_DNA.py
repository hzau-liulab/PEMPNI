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
    train_data = fil +'/datasets/MPD276.txt'
    train_bio_fea = fil+'/features/bio_feature_DNA.txt'
    train_inter_intra = fil+"/features/EPI.txt"
    train_target = fil+"/features/ETOR.txt"
    
    bio_fea = fil+'/output/bio_feature/bio_feature.txt'
    inter_intra = fil+"/output/energy_feature/GB1/EPI.txt";
    target = fil+"/output/energy_feature/GB1/ETOR.txt"

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
    samp_inter_fea = read_feature(train_inter_intra)
    it_inter_fea = np.loadtxt(inter_intra, dtype=bytes).astype(str)

    train_inter_X = []
    train_inter_Y = []
    for k in train_samp_ddg:
        train_inter_Y.append(float(k[1]))
        train_inter_X.append(samp_inter_fea[k[0]])

    test_inter_X = []
    test_inter_X.append(it_inter_fea[:])
    
    model_line = RandomForestRegressor(n_estimators=500)
    model_line.fit(train_inter_X, train_inter_Y)
    prediction_inter = model_line.predict(test_inter_X)

###############################################################################
    samp_tar_fea = read_feature(train_target)
    it_tar_fea = np.loadtxt(target, dtype=bytes).astype(str)

    train_tar_X = []
    train_tar_Y = []
    for k in train_samp_ddg:
        train_tar_Y.append(float(k[1]))
        train_tar_X.append(samp_tar_fea[k[0]])

    test_tar_X = []
    test_tar_X.append(it_tar_fea[:])
    
    model_line = RandomForestRegressor(n_estimators=500)
    model_line.fit(train_tar_X, train_tar_Y)
    prediction_tar = model_line.predict(test_tar_X)

    en_predition = 0.4*prediction_tar + 0.6*prediction_inter
    prediction = 0.6*en_predition + 0.4*prediction_bio

    with open(fil+"/output/result.txt","a") as fil:
        fil.write("\n"+str("Predicted affinity change:")+str('%.3f' % prediction[0]))
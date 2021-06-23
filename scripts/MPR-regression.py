import numpy as np
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
import sys
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.metrics import mean_squared_error
import math
import re
import joblib

if __name__ == '__main__':
    workdir = sys.argv[1]
    fil = sys.argv[2]
    EWC_regression = joblib.load(workdir+'/models/MPR-EWC-Regression.m')
    Nonenergy_regression = joblib.load(workdir+'/models/MPR-Nonenergy-Regression.m')
    
    nonenergy_fea = np.loadtxt(workdir+fil+'/bio_feature/bio_feature.txt',dtype=bytes).astype(str)
    EWC_fea = np.loadtxt(workdir+fil+"/energy_feature/GB2/EWC.txt",dtype=bytes).astype(str)

    test_nonenergy = []
    test_nonenergy.append(nonenergy_fea[1:])
    prediction_nonenergy = Nonenergy_regression.predict(test_nonenergy)

    test_EWC = []
    test_EWC.append(EWC_fea[:])
    prediction_energy = EWC_regression.predict(test_EWC)

    prediction = 0.5*prediction_energy + 0.5*prediction_nonenergy

    with open(workdir+fil+"/result.txt","a") as fil:
        fil.write("\n"+str("Predicted affinity change:")+str('%.3f' % prediction[0]))
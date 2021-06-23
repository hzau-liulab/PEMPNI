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
    EPI_regression = joblib.load(workdir+'/models/MPD-EPI-Regression.m')
    ETOR_regression = joblib.load(workdir+'/models/MPD-ETOR-Regression.m')
    Nonenergy_regression = joblib.load(workdir+'/models/MPD-Nonenergy-Regression.m')
    
    nonenergy_fea = np.loadtxt(workdir+fil+'/bio_feature/bio_feature.txt',dtype=bytes).astype(str)
    EPI_fea = np.loadtxt(workdir+fil+"/energy_feature/GB1/EPI.txt",dtype=bytes).astype(str)
    ETOR_fea = np.loadtxt(workdir+fil+"/energy_feature/GB1/ETOR.txt",dtype=bytes).astype(str)

    test_nonenergy = []
    test_nonenergy.append(nonenergy_fea[1:])
    prediction_nonenergy = Nonenergy_regression.predict(test_nonenergy)

    test_EPI = []
    test_EPI.append(EPI_fea[:])
    prediction_EPI = EPI_regression.predict(test_EPI)

    test_ETOR = []
    test_ETOR.append(ETOR_fea[:])
    prediction_ETOR = ETOR_regression.predict(test_ETOR)

    predition_energy = 0.4*prediction_ETOR + 0.6*prediction_EPI
    prediction = 0.6*predition_energy + 0.4*prediction_nonenergy

    with open(workdir+fil+"/result.txt","a") as fil:
        fil.write("\n"+str("Predicted affinity change:")+str('%.3f' % prediction[0]))
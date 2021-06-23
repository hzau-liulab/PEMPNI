import numpy as np
from sklearn.ensemble import RandomForestClassifier
import sys
import math
import re
import joblib

if __name__ == '__main__':
    workdir = sys.argv[1]
    fil = sys.argv[2]
    EPI_classifier = joblib.load(workdir+'/models/MPD-EPI-Classification.m')
    ETOR_classifier = joblib.load(workdir+'/models/MPD-ETOR-Classification.m')
    Nonenergy_classifier = joblib.load(workdir+'/models/MPD-Nonenergy-Classification.m')
    
    nonenergy_fea = np.loadtxt(workdir+fil+'/bio_feature/bio_feature.txt',dtype=bytes).astype(str)
    EPI_fea = np.loadtxt(workdir+fil+"/energy_feature/GB1/EPI.txt",dtype=bytes).astype(str)
    ETOR_fea = np.loadtxt(workdir+fil+"/energy_feature/GB1/ETOR.txt",dtype=bytes).astype(str)
    
    test_nonenergy = []
    test_nonenergy.append(nonenergy_fea[1:])
    prediction_nonenergy = Nonenergy_classifier.predict_proba(test_nonenergy)
    
    test_EPI = []
    test_EPI.append(EPI_fea[:])
    prediction_EPI = EPI_classifier.predict_proba(test_EPI)

    test_ETOR = []
    test_ETOR.append(ETOR_fea[:])
    prediction_ETOR = ETOR_classifier.predict_proba(test_ETOR)
    
    ###########################################################################
    predition_energy = 0.4*prediction_ETOR[:,1] + 0.6*prediction_EPI[:,1]
    prediction = 0.6*predition_energy  + 0.4*prediction_nonenergy[:,1]
    
    with open(workdir+fil+"/result.txt","a") as fil:
        fil.write("\n"+"Predicted score:"+str('%.3f' % prediction[0]))
        if prediction[0]>=0.471:
            fil.write("\n"+"Significant decrease:YES")
        else:
            fil.write("\n"+"Significant decrease:NO")
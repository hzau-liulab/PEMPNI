import numpy as np
from sklearn.ensemble import RandomForestClassifier
import sys
import math
import re
import joblib

if __name__ == '__main__':
    workdir = sys.argv[1]
    
    EWC_classifier = joblib.load(workdir+'/models/MPR-EWC-Classification.m')
    Nonenergy_classifier = joblib.load(workdir+'/models/MPR-Nonenergy-Classification.m')
   
    nonenergy_fea = np.loadtxt(workdir+'/output/bio_feature/bio_feature.txt',dtype=bytes).astype(str)
    EWC_fea = np.loadtxt(workdir+"/output/energy_feature/GB2/EWC.txt",dtype=bytes).astype(str)
    
    test_nonenergy = []
    test_nonenergy.append(nonenergy_fea[1:])
    prediction_nonenergy = Nonenergy_classifier.predict_proba(test_nonenergy)

    test_EWC = []
    test_EWC.append(EWC_fea[:])
    prediction_energy = EWC_classifier.predict_proba(test_EWC)

    prediction = 0.5*prediction_energy[:,1] + 0.5*prediction_nonenergy[:,1]
     
    with open(workdir+"/output/result.txt","a") as fil:
        fil.write("\n"+"Predicted score:"+str('%.3f' % prediction[0]))
        if prediction[0]>=0.320:
            fil.write("\n"+"Significant decrease:YES")
        else:
            fil.write("\n"+"Significant decrease:NO") 

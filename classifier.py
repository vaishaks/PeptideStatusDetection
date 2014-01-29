"""
Creates a classifier for predicting the peptide status.
"""

import create_datasets as cd
import os
import numpy as np
import pandas as pd
import pylab as pl
from itertools import count
from sklearn import svm
from sklearn import cross_validation as cv
from sklearn.externals import joblib
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import roc_curve, auc, confusion_matrix

datasets_index = ["amylnset", "pafig", "zipper", "amylpred"]
global_cm = []
global_auc = []
amylpred_cm = []
amylpred_auc = []
fig_count = count()

def generic_svm_classifier(X, y, dataset="amylnset", n=6):
    # Split the data into training and test.
    X_train, X_test, y_train, y_test = cv.train_test_split(X, y, test_size=0.2, 
                                                       random_state=0)
    # Using GridSearchCV to find the best values for C and gamma
    C_range = 10.0 ** np.arange(-2, 4)
    gamma_range = 10.0 ** np.arange(-5, 1)
    param_grid = dict(gamma=gamma_range, C=C_range)
    skf = cv.StratifiedKFold(y=y_train, n_folds=3)
    grid = GridSearchCV(svm.SVC(), param_grid=param_grid, cv=skf)
    grid.fit(X_train, y_train)
    
    crossclf = svm.SVC(probability=True, **grid.best_params_)
    print "Cross-Validation score", cv.cross_val_score(crossclf, X_train, 
                                                        y_train, cv=5).mean()
    print "Independent test score", crossclf.fit(X_train, y_train).score(X_test
                                                                    , y_test)
    
    # Compute roc and auc
    probas_ = crossclf.predict_proba(X_test)
    fpr, tpr, thresholds = roc_curve(y_test, probas_[:, 1])
    roc_auc = auc(fpr, tpr)
    print "Area under the curve", roc_auc
    
    # Plot roc curve
    pl.figure(next(fig_count))
    pl.clf()
    pl.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
    pl.plot([0, 1], [0, 1], 'k--')
    pl.xlim([0.0, 1.0])
    pl.ylim([0.0, 1.0])
    pl.xlabel('False Positive Rate')
    pl.ylabel('True Positive Rate')
    pl.title('Receiver operating characteristic')
    pl.legend(loc="lower right")
    pl.show()
    
    # Confusion Matrix
    y_pred = crossclf.predict(X_test)
    cm = confusion_matrix(y_test, y_pred)
    print "The confusion matrix:"
    print cm
    
    if dataset == "amylpred" and n == 6:
        global_cm.append(cm)
        global_auc.append(roc_auc)
        amylpred_cm.append(cm)
        amylpred_auc.append(roc_auc)
    elif dataset == "amylpred" and n > 6:
        amylpred_cm.append(cm)
        amylpred_auc.append(roc_auc)
    else:
        global_cm.append(cm)
        global_auc.append(roc_auc)                
    
    clf = svm.SVC(**grid.best_params_) # Unpack the best params found
    clf.fit(X, y)
    # Save the model for future use. Saves computing time.
    joblib.dump(clf, "data/temp/"+dataset+str(n)+".pkl")

def create_confusion_csv():
    tp = []
    tn = []
    fp = []
    fn = []
    for x in global_cm:
        tp.append(x[0][0])
        tn.append(x[1][1])
        fp.append(x[1][0])
        fn.append(x[0][1])

    d = {"TP": tp, "TN": tn, "FP": fp, "FN": fn, "AUC": global_auc}
    confusion_dataframe = pd.DataFrame(d, index=datasets_index)
    confusion_dataframe.to_csv("data/temp/confusion.csv")
    print "Confusion csv created."
    
    tp = []
    tn = []
    fp = []
    fn = []
    for x in amylpred_cm:
        tp.append(x[0][0])
        tn.append(x[1][1])
        fp.append(x[1][0])
        fn.append(x[0][1])

    d = {"TP": tp, "TN": tn, "FP": fp, "FN": fn, "AUC": amylpred_auc}
    confusion_dataframe = pd.DataFrame(d, index=["amylpred6set", 
                                                "amylpred7set", "amylpred8set"])
    confusion_dataframe.to_csv("data/temp/amylpred_confusion.csv")    
    print "Amylpred Confusion csv created."    



n = 6
if os.path.exists("data/temp/amylnset.pkl"):
    print "Using a pre-trained classifier.."
    clf = joblib.load("data/temp/amylnset.pkl")

    # Testing on a standard dataset
    
else:
    # Training with amylnset data
    print "Creating amylnset.."
    cd.create_amylnset(n)
    print "Training the classifier.."
    # Extracting features and labels from the dataset.
    X = []
    y = []
    data = open("data/temp/amyl"+str(n)+"set.txt")
    for line in data:
        temp = line.rstrip().split()
        y.append(int(temp[1]))
        X.append(map(float, temp[2:]))
    data.close()

    generic_svm_classifier(X, y, "amylnset")

    # Training with pafig data
    print "Creating pafig dataset.."
    cd.create_pafig_data(n)
    print "Training the classifier.."
    # Extracting features and labels from the dataset.
    X = []
    y = []
    data = open("data/temp/pafig_hexpepset.txt")
    for line in data:
        temp = line.rstrip().split()
        y.append(int(temp[1]))
        X.append(map(float, temp[2:]))
    data.close()

    generic_svm_classifier(X, y, "pafig")
    
    # Training with zipper data
    print "Creating zipper dataset.."
    cd.create_zipper_data(n)
    print "Training the classifier.."
    # Extracting features and labels from the dataset.
    X = []
    y = []
    data = open("data/temp/zipper_hexpepset.txt")
    for line in data:
        temp = line.rstrip().split()
        y.append(int(temp[1]))
        X.append(map(float, temp[2:]))
    data.close()

    generic_svm_classifier(X, y, "zipper")
    
    # Training with amylpred data
    print "Creating amylpred dataset.."
    cd.create_amylpred_data(n)
    cd.create_amylpred_data(n+1)
    cd.create_amylpred_data(n+2)
    print "Training the classifier.."
    # Extracting features and labels from the dataset.
    X = []
    y = []
    data = open("data/temp/amylpred"+str(n)+"set.txt")
    for line in data:
        temp = line.rstrip().split()
        y.append(int(temp[1]))
        X.append(map(float, temp[2:]))
    data.close()

    generic_svm_classifier(X, y, "amylpred")
    
    # Extracting features and labels from the dataset.
    X = []
    y = []
    data = open("data/temp/amylpred"+str(n+1)+"set.txt")
    for line in data:
        temp = line.rstrip().split()
        y.append(int(temp[1]))
        X.append(map(float, temp[2:]))
    data.close()

    generic_svm_classifier(X, y, "amylpred", n+1)

    # Extracting features and labels from the dataset.
    X = []
    y = []
    data = open("data/temp/amylpred"+str(n+2)+"set.txt")
    for line in data:
        temp = line.rstrip().split()
        y.append(int(temp[1]))
        X.append(map(float, temp[2:]))
    data.close()

    generic_svm_classifier(X, y, "amylpred", n+2)
            
    create_confusion_csv()
# TODO Predicting the amyloidogenic regions in a protein sequence in fasta format.

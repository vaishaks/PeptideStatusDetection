"""
Creates a classifier for predicting the peptide status.
"""

import create_datasets as cd
import os
import numpy as np
import pylab as pl
from sklearn import svm
from sklearn import cross_validation as cv
from sklearn.externals import joblib
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import roc_curve, auc, confusion_matrix

def generic_svm_classifier(X, y, dataset="amylnset"):
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
    
    # Plot the confusion matrix
    pl.matshow(cm)
    pl.title('Confusion Matrix')
    pl.colorbar()
    pl.ylabel('Amyloidogenic')
    pl.xlabel('Non-Amyloidogenic')
    pl.show()
    
    clf = svm.SVC(**grid.best_params_) # Unpack the best params found
    clf.fit(X, y)
    # Save the model for future use. Saves computing time.
    joblib.dump(clf, "data/temp/"+dataset+".pkl")

n = int(raw_input("Enter the size of the window: "))
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
    print "Training the classifier.."
    # Extracting features and labels from the dataset.
    X = []
    y = []
    data = open("data/temp/amylpred_hexpepset.txt")
    for line in data:
        temp = line.rstrip().split()
        y.append(int(temp[1]))
        X.append(map(float, temp[2:]))
    data.close()

    generic_svm_classifier(X, y, "amylpred")

# TODO Predicting the amyloidogenic regions in a protein sequence in fasta format.


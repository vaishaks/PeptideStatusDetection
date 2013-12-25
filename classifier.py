"""
Creates a classifier for predicting the peptide status.
"""

import create_amylnset as ca
import os
import numpy as np
from sklearn import svm
from sklearn import cross_validation as cv
from sklearn.externals import joblib
from sklearn.grid_search import GridSearchCV

n = int(raw_input("Enter the size of the window: "))
if os.path.exists("data/temp/amyl"+str(n)+"pred.pkl"):
    print "Using a pre-trained classifier.."
    clf = joblib.load("data/temp/amyl"+str(n)+"pred.pkl")
else:
    print "Creating data.."
    ca.create_amylnset(n)
    print "Training the classifier.."
    # Extracting features and labels from the dataset.
    X= []
    y = []
    data = open("data/temp/amyl"+str(n)+"set.txt")
    for line in data:
        temp = line.rstrip().split()
        y.append(int(temp[1]))
        X.append(map(float, temp[2:]))
    data.close()
    
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
    crossclf = grid.best_estimator_
    print "Cross-Validation score", cv.cross_val_score(crossclf, X_train, y_train, cv=5).mean()
    print "Independent test score", crossclf.score(X_test, y_test)
    clf = svm.SVC(**grid.best_params_) # Unpack the best params found
    clf.fit(X, y)
    # Save the model for future use. Saves computing time.
    joblib.dump(clf, "data/temp/amyl"+str(n)+"pred.pkl")

# Predicting the amyloidogenic regions in a protein sequence in fasta format.
ip = open("data/input.fasta")
op = open("data/temp/output.txt", "w")
ip.readline()
seq = ""
for line in ip:
    seq += line.rstrip()
for i in xrange(len(seq)-n):
    window = seq[i:i+n]
    if clf.predict(ca.compute_features(window))[0]:
        op.write(window + " " + str(i) + " " + str(i+n) + "\n")
op.close()

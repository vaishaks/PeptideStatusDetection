"""
Creates a classifier for predicting the peptide status.
"""

import create_amylnset as ca
import os
import numpy as np
from sklearn import svm
from sklearn import cross_validation as cv
from sklearn.externals import joblib

n = int(raw_input("Enter the size of the window: "))
if os.path.exists("data/temp/amyl"+str(n)+"pred.pkl"):
    print "Using a pre-trained classifier.."
    clf = joblib.load("data/temp/amyl"+str(n)+"pred.pkl")
else:
    print "Creating data.."
    ca.create_amylnset(n)
    print "Training the classifier.."
    X= []
    y = []
    data = open("data/temp/amyl"+str(n)+"set.txt")
    for line in data:
        temp = line.rstrip().split()
        y.append(int(temp[1]))
        X.append(map(float, temp[2:]))
    data.close()
    X_train, X_test, y_train, y_test = cv.train_test_split(X, y, test_size=0.2, 
                                                       random_state=0)
    max_score = 0
    c = 1
    for i in np.arange(1.0, 1.4, 0.08):
        crossclf = svm.SVC(C=i).fit(X_train, y_train)
        scores = cv.cross_val_score(crossclf, X_train, y_train, cv=5)
        if scores.mean() > max_score:
            max_score = scores.mean()
            c = i
    g = 0
    for i in np.arange(0.0, 0.004, 0.0008):
        crossclf = svm.SVC(C=c, gamma=i).fit(X_train, y_train)
        scores = cv.cross_val_score(crossclf, X_train, y_train, cv=5)
        if scores.mean() > max_score:
            max_score = scores.mean()
            g = i

    print "Cross-Validation score", max_score      
    crossclf = svm.SVC(C=c, gamma=g).fit(X_train, y_train)      
    print "Independant test score", crossclf.score(X_test, y_test)

    clf = svm.SVC(C=c)
    clf.fit(X, y)
    joblib.dump(clf, "data/temp/amyl"+str(n)+"pred.pkl")

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

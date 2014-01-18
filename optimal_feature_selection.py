import os
import create_datasets as cd
from sklearn import svm
from sklearn import cross_validation as cv

n = int(raw_input("Enter the size of the window: "))
print "Creating data.."
cd.create_amylnset(n)

# Extracting features and labels from the dataset.
X= []
y = []
data = open("data/temp/amyl"+str(n)+"set.txt")
for line in data:
    temp = line.rstrip().split()
    y.append(int(temp[1]))
    X.append(map(float, temp[2:]))
data.close()

print len(X[0])
feature_scores = []
"""
# Computing scores taking single features at a time other than the correlation
for i in xrange(len(X)-n):
    clf = svm.SVC();
    feature_scores.append(cv.cross_val_score(clf, X[i], y, cv=5).mean())
    
print feature_scores[:10], len(feature_scores)
"""
import os
import numpy as np
import pandas as pd
import create_datasets as cd
from sklearn.naive_bayes import GaussianNB
from sklearn import cross_validation as cv

n = int(raw_input("Enter the size of the window: "))
if not os.path.exists("data/temp/amyl"+str(n)+"set.txt"):
    print "Creating data.."
    cd.create_amylnset(n)

# Extracting features and labels from the dataset.
X = []
y = []
data = open("data/temp/amyl"+str(n)+"set.txt")
for line in data:
    temp = line.rstrip().split()
    y.append(int(temp[1]))
    X.append(map(float, temp[2:]))
data.close()

X = np.array(X)
y = np.array(y)
Xtrans = X.T

feature_scores = []
# Computing scores taking single features at a time other than the correlation
for i in xrange(len(Xtrans)-n):
    # Using naive bayes because it's a faster algorithm, check for disadvantages
    clf = GaussianNB()
    feature_scores.append(cv.cross_val_score(clf, [[x] for x in Xtrans[i]], 
                                                y, cv=5).mean())
# Get the names of the features   
feature_names = []
aaindex = open("data/temp/amino_acid_index.txt")
for line in aaindex:
    feature_names.append(line.split()[0])

# Create a pandas series of the feature scores and names, sort and save it to csv
feature_series = pd.Series(feature_scores, index=feature_names)
feature_series.sort(ascending=False)
feature_series.to_csv("data/temp/feature_series.csv")

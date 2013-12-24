"""
Creates a classifier for predicting the peptide status.
"""

import create_amylnset as ca
from sklearn import svm
from sklearn import cross_validation as cv

n = raw_input("Enter the size of the window: ")
print "Creating data.."
ca.create_amylnset(int(n))
print "Training the classifier.."
X= []
y = []
data = open("data/temp/amyl"+str(n)+"set.txt")
for line in data:
    temp = line.rstrip().split()
    y.append(int(temp[1]))
    X.append(map(float, temp[2:]))

X_train, X_test, y_train, y_test = cv.train_test_split(X, y, test_size=0.2, 
                                                       random_state=0)
max_score = 0
c = 1
for i in xrange(1, 5, 1):
    crossclf = svm.SVC(C=i).fit(X_train, y_train)
    scores = cv.cross_val_score(crossclf, X_train, y_train, cv=5)
    if scores.mean() > max_score:
        max_score = scores.mean()
        c = i
        
crossclf = svm.SVC(C=c).fit(X_train, y_train)      
print "Classifier trained with accuracy", crossclf.score(X_test, y_test)

clf = svm.SVC(C=c)
clf.fit(X, y)

print "Enter the n-peptide sequence: ",
seq = raw_input()
print ""
if clf.predict(ca.compute_features(seq))[0]:
    print "The n-peptide is amyloidogenic."
else:
    print "The n-peptide is non-amyloidogenic."

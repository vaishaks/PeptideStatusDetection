from sklearn.naive_bayes import GaussianNB
from sklearn import cross_validation as cv

print "Testing with amyl6set"

data = open('data/temp/amyl6set.txt','r')
X= []
y = []

#data = open("data/temp/amyl"+str(n)+"set.txt")

for line in data:
	temp = line.rstrip().split()
	y.append(int(temp[1]))
	X.append(map(float, temp[2:]))
data.close()

X_train, X_test, y_train, y_test = cv.train_test_split(X, y, test_size=0.2, 
                                                       random_state=0)

gnb = GaussianNB()
gnb.fit(X_train, y_train)
print "Independent test score", gnb.score(X_test, y_test)

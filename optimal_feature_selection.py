import numpy as np
import pandas as pd
from sklearn.naive_bayes import GaussianNB
from sklearn import cross_validation as cv

def select_optimal_features(n, dataset="amylnset"):
    # Extracting features and labels from the dataset.
    X = []
    y = []
    try:
    	if dataset == "amylnset":
        	data = open("data/temp/temp_amyl"+str(n)+"set.txt")
    	elif dataset == "pafig":
        	data = open("data/temp/temp_pafig_hexpepset.txt")
    	elif dataset == "zipper":
        	data = open("data/temp/temp_zipper_hexpepset.txt")    
        elif dataset == "amylpred":
        	data = open("data/temp/temp_amylpred"+str(n)+"set.txt")     
    	else:
        	data = open("data/temp/temp_amyl"+str(n)+"set.txt")
    except IOError:
    	print "Error: Temporary dataset not found."
    	return
        
    for line in data:
        temp = line.rstrip().split()
        y.append(int(temp[1]))
        X.append(map(float, temp[2:]))
    data.close()
    X = np.array(X)
    y = np.array(y)

    feature_scores = []
    feature_id = []
    #return
    # Computing scores taking single features at a time other than the correlation
    print "Training one feature at a time."
    for i in xrange(len(X[0])-n):
        # Using naive bayes because it's a faster algorithm, check for disadvantages
        clf = GaussianNB()
        try:
            feature_scores.append(cv.cross_val_score(clf, np.array([[x[i]] for x in X]), y, cv=5).mean())
        except:
            print [i for i in xrange(len(X)) if len(X[i]) < 566]
        feature_id.append(i)
    # Get the names of the features   
    feature_names = []
    aaindex = open("data/temp/amino_acid_index.txt")
    for line in aaindex:
        feature_names.append(line.split()[0])

    # Sort features based on score and save to csv
    d = {"scores": feature_scores, "id": feature_id}
    feature_dataframe = pd.DataFrame(d, index=feature_names)
    feature_dataframe = feature_dataframe.sort('scores', ascending=False)
    feature_dataframe.to_csv("data/temp/"+dataset+"_feature_dataframe.csv")
    print ".csv created."

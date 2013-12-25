Exploiting heterogeneous features to improve in silico prediction of peptide status - amyloidogenic or non-amyloidogenic
=========================================================================================================================

## Includes:
* 'create_amylnset.py': Python script used to extract positive data regions and prepare dataset for classification from the protein sequences in fasta format.
* 'classifier.py': Classifier written in python script using scikit-learn to predict the amyloidogenic regions of size 'n' in a protein sequence.
* Data: Protein sequences, both positive and negative regions.

## Dependancies:
* Scikit-Learn
* Numpy
* Pylab

## Instructions:
* Make sure that all dependancies are pre-installed
* Clone the repository as .zip
* Extract the contents
* The protein sequence for which the amyloidogenic regions are to be predicted in .fasta format is to be renamed as input.fasta and placed in the 'data' directory.
* Run the 'classifier.py' script and input the window size.
* Clean the 'data/temp' directory if required.

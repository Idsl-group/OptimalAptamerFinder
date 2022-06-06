from src.Sequence import *

# Read sequencelib from pickle file
with open('/home/javier/PycharmProjects/OptimalAptamerFinder/results/2022-06-06 11:35:29.596923/sequencelib.pickle', 'rb') as f:
    sequencelib = pickle.load(f)
f.close()
import pandas as pd
import numpy as np
import pickle
import sys
import csv
import os

#run by passing in your mpraflow dictionary tsv file
myfile=sys.argv[1]
basename=os.path.splitext(myfile)[0]
outputfile=basename + ".pickle"

#read in association tsv file, with only columns id and barcode
assocs = pd.read_csv(myfile, sep='\t')

#turn into dictionary
bcdict = assocs.groupby('id')['barcode'].apply(set).to_dict()

#write to pickle file
with open(outputfile, 'wb') as fp:
    pickle.dump(bcdict, fp)
    print('dictionary saved to file')

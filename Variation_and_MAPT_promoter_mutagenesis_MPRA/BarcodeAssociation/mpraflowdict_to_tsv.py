import pandas as pd
import sys
import csv
import pickle
import os

#run by passing in your mpraflow dictionary
#should be something like *_filtered_coords_to_barcodes.pickle

myfile=sys.argv[1]
basename=os.path.splitext(myfile)[0]
outputfile=basename + ".tsv"

#function for opening pickle files
def read_pickle_file(file):
    pickle_data = pd.read_pickle(file)
    return pickle_data

#open my dictionary from mpraflow, and save as a dictionary object
mpraflowdict = open(myfile, "rb")
mpraflowdict2 = pd.read_pickle(mpraflowdict)

#turn my dictionary into a dataframe
dict_df = pd.DataFrame.from_dict(mpraflowdict2, orient='index')
dict_df.reset_index(inplace=True)
dict_df.rename(columns={'index': 'id'}, inplace=True)

#each column with a barcode in it was given a number in the previous step
#get a list of every column name, saved to a variable, removing the id column (only barcode columns)
mycolumns = list(dict_df.columns)
mycolumns.remove('id')

#pivot dataframe longer so that there's only 2 columns, with every barcode in a row and their matched ids
dict_df2 = dict_df.melt(id_vars='id', value_vars=mycolumns, var_name='barcode_id', value_name='barcode')
dict_df3 = dict_df2[['id', 'barcode']]

#remove rows with NA values
dict_df4=dict_df3.dropna()

#save dataframe as a tab separated file
dict_df4.to_csv(outputfile, header=True, index=False, sep="\t", quoting=csv.QUOTE_NONE, escapechar=' ')



import csv
import sys

myfile=sys.argv[1]
outputfile=sys.argv[2]

#open tsv file to pad
with open(myfile, 'r', newline='', encoding='utf-8') as tsvfile:
    reader=csv.reader(tsvfile, delimiter='\t')

    #get the number of expected columns from the length of the first row
    first_row = next(reader, None)
    colnumber = len(first_row)

    #saving output to a csv files
    with open(outputfile, 'w', newline='', encoding='utf-8') as csvfile:
        writer=csv.writer(csvfile)
        #write header row
        writer.writerow(first_row)

        #iterate through each row in the files
        #pad the row with 0's if the number of columns is less than expected
        for row in reader:
            row += ['0'] * (colnumber - len(row))
            #row = list(map(int,row))
            writer.writerow(row)

import csv
import pandas as pd

'''
Function to find orthologous and co-orthologous protein sequences using Best Bidirectional Hits method
input = blastp results sorted by increasing e-value for each query protein
output = file containing Best Bidirectional Hits; two files with in-paralogs for each of two species
'''

def orthologs(blastp_1vs2, blastp_2vs1, blastp_1vs1, blastp_2vs2):
#create empty dataframes for each blastp query
    df1 = pd.DataFrame(columns=['MP','GO','eval'])
    df2 = pd.DataFrame(columns=['GO','MP','eval'])
    df3 = pd.DataFrame(columns=['MP1','MP2','eval_ip','in_paralog'])
    df4 = pd.DataFrame(columns=['GO1','GO2','eval_ip','in_paralog'])

#for all files, choose only the match with the lowest e-value for each protein pair
    with open(blastp_1vs2) as my_file1:
        tsv_reader = csv.reader(my_file1, delimiter="\t")
        id = ''
        for row in tsv_reader: 
            if row[0] == id:
                continue
            else:
                new_row={'MP':row[0],'GO':row[1],'eval':row[10]}
                df1 = df1.append(new_row, ignore_index=True)
                id = row[0]

    with open(blastp_2vs1) as my_file2:
        tsv_reader = csv.reader(my_file2, delimiter="\t")
        id = ''
        for row in tsv_reader: 
            if row[0] == id:
                continue
            else:
                new_row={'GO':row[0],'MP':row[1],'eval':row[10]}
                df2 = df2.append(new_row, ignore_index=True)
                id = row[0]
     
    with open(blastp_1vs1) as my_file3:
        tsv_reader = csv.reader(my_file3, delimiter="\t")
        id = ''
        for row in tsv_reader: 
            #ignore rows that match with themselves
            if row[0] == row[1]:
                continue
            elif row[0] == id:
                continue
            else:
                new_row={'MP1':row[0],'MP2':row[1],'eval_ip':row[10],'in_paralog':False}
                df3 = df3.append(new_row, ignore_index=True)
                id = row[0]

    with open(blastp_2vs2) as my_file4:
        tsv_reader = csv.reader(my_file4, delimiter="\t")
        id = ''
        for row in tsv_reader: 
            #ignore rows that match with themselves
            if row[0] == row[1]:
                continue
            elif row[0] == id:
                continue
            else:
                new_row={'GO1':row[0],'GO2':row[1],'eval_ip':row[10],'in_paralog':False}
                df4 = df4.append(new_row, ignore_index=True)
                id = row[0]
    
    #find pairs appearing in both best hit files
    column_titles = ['MP','GO']
    df2_new=df2.reindex(columns=column_titles)

    BBH = pd.concat([df1,df2_new]).groupby(['GO','MP']).size()
    BBH = BBH[BBH == 2]

    #find in-paralogs 
    #...in species 1
    for index1, row1 in df3.iterrows():
        for index2, row2 in df1.iterrows():
            if row1['MP1'] == row2['MP'] and float(row1['eval_ip']) < float(row2['eval']):
                row1['in_paralog'] = True
            else:
                continue
            
    #...in species 2
    for index1, row1 in df4.iterrows():
        for index2, row2 in df2.iterrows():
            if row1['GO1'] == row2['GO'] and float(row1['eval_ip']) < float(row2['eval']):
                row1['in_paralog'] = True
            else:
                continue

    return BBH, df3, df4









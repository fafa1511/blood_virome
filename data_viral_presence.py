"""
Get dataframe for all viral contig presence file
"""
import glob
import csv
import sys
import functions as fct
import pandas as pd
import numpy as np

# find all file adundance
filex = glob.glob("Plasma*")
filex.sort()


# dico_ab_viral = fct.create_ab_viral_dico(abundance_viral)
# dico_data = {}
# dico_data['ab'] = []
# for valeur in dico_ab_viral.values():
#     dico_data['ab'].append(valeur) 
# print(dico_data['ab'][5])

dico_data = {}
for j in range(len(filex)):
    sample = filex[j].split("_dico_univ")[0]
    with open(filex[j],"r") as fi_viral:
 	    viral_presence = fi_viral.readlines()
    dico_viral_presence = fct.create_presence_dico(viral_presence)
    dico_data[sample] = []
    for valeur in dico_viral_presence.values():
        dico_data[sample].append(valeur) 


ind = []
for contig in dico_viral_presence:
    ind.append(contig)

data = pd.DataFrame(data = dico_data, index = ind)

data.to_csv('viral_presence_plasma.csv', sep = '\t')

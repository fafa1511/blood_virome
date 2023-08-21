"""
Dataframe for new tables of abundance
    - input: all new tables of abundance
    - output: a dataframe for all new tables of abundance

example: python /work_projet/phages/hao/blood_virome/scripts/data_ab.py
"""
import glob
import csv
import sys
import functions as fct
import pandas as pd
import numpy as np

# find all file adundance
filex = glob.glob("*_dico_new_ab*")
filex.sort()


# dico_ab_viral = fct.create_ab_viral_dico(abundance_viral)
# dico_data = {}
# dico_data['ab'] = []
# for valeur in dico_ab_viral.values():
#     dico_data['ab'].append(valeur) 
# print(dico_data['ab'][5])

dico_data = {}
for j in range(len(filex)):
    sample = filex[j].split("_dico_new_ab")[0]
    with open(filex[j],"r") as fi_ab:
 	    abundance_presence = fi_ab.readlines()
    dico_ab_presence = fct.create_new_ab_dico(abundance_presence)
    dico_data[sample] = []
    for valeur in dico_ab_presence.values():
        dico_data[sample].append(valeur) 


ind = []
for contig in dico_ab_presence:
    ind.append(contig)

data = pd.DataFrame(data = dico_data, index = ind)

data.to_csv('presence_abundance.csv', sep = '\t')



# -*- coding: utf-8 -*-
"""
list of present contig in files

    - input: presence.csv files, output path
    - output: list of contig present

example: list.py files presence.csv /work_projet/phages/hao/blood_virome/07-Presence/list_contig

"""
import glob
import csv
import sys


def lists(li):
    li_1 = []

    #for each contig
    for i in range(len(li)):

        # get contig information
        info_li = li[i].split()
        sample_name = info_li[0]
        ab = int(info_li[1])

        # if present
        if (ab == 1):
            li_1.append([sample_name])
    
    return li_1


file_li = sys.argv[1]
path = sys.argv[2]
sample = file_li.split("_a")[0]

# read csv files
with open(file_li,"r") as fi:
 	    li = fi.readlines()


res = lists(li)
# write csv file
with open('%s/%s_li.csv' %(path, sample), 'a', newline='') as csvfile:
    				spamwriter = csv.writer(csvfile, delimiter=' ',
    		                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
    				spamwriter.writerows(res)





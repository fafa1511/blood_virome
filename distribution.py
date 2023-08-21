# -*- coding: utf-8 -*-
"""
Distribution of contig length

    - input: all presence.csv files, length file
    - output: list of contig present with length

exemple: python distribution.py file_ab /work_projet/phages/hao/blood_virome/07-Presence

"""
import csv
import sys



def distribution(li_ab, li_len, li_dist):

    # for each contig in each presence.csv
    for i in range(len(li_len)):
        # get information
        info_ab = li_ab[i].split()
        info_len = li_len[i].split()

        # if same order of contig
        if (info_ab[0] == info_len[0]):
            sample_name = info_len[0]                    # get sample name
            ab = int(info_ab[1])                         # get presence or absence value
            length = int(info_len[1])                    # get length
            if (ab == 1):                                # if present
                li_dist.append([sample_name, length])    # add sample name with length

        else:
            print("error order")

    return li_dist




li_dist = []

file_ab = sys.argv[1]
file_len = sys.argv[2]
sample = file_ab.split("_abundance")[0]


# read csv files
with open(file_ab,"r") as fi_ab:
 	li_ab = fi_ab.readlines()

with open(file_len,"r") as fi_len:
 	li_len = fi_len.readlines()



li = distribution(li_ab, li_len, li_dist)
# write csv file
with open('%s_distribution.csv' % (sample), 'a', newline='') as csvfile:
				spamwriter = csv.writer(csvfile, delimiter=' ',
		                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
				spamwriter.writerows(li)
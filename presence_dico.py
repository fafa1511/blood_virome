# -*- coding: utf-8 -*-
"""
Determine the presence or absence of a contig

    - input: all covrage.txt files, threshold of covrage values, output path
    - output: all presence.csv files

example: python presence.py covrage file, 0.5 /work_projet/phages/hao/blood_virome/07-Presence
"""
import csv
import sys
import functions as fct


file_cov = sys.argv[1]
threshold_cov = float(sys.argv[2])
path_out = sys.argv[3]

sample = file_cov.split(".")[0]


with open(file_cov,"r") as fi_cov:
	coverage = fi_cov.readlines()


res = fct.create_dico_presence(fct.create_dico_cov(coverage), threshold_cov).items()

with open('%s/%s_dico_presence.csv' % (path_out, sample), 'a', newline='') as csvfile:
				spamwriter = csv.writer(csvfile, delimiter=' ',
		                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
				spamwriter.writerows(res)


                

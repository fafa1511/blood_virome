# -*- coding: utf-8 -*-
"""
Keep only viral contigs

    - input: presence.csv files after decontamination, viral list, output path
    - output: presence.csv files after decontamination with only viral contig

example: python univ_viral.py Univ_viral_contigs.txt presence.csv files output path

"""
import csv
import sys
import functions as fct

file_viral = sys.argv[1]
file_sample_decontamine = sys.argv[2]
path_output = sys.argv[3]

sample_name = file_sample_decontamine.split("_ConNeg23")[0]

# read csv files
with open(file_viral,"r") as fi_vi:
 	viral = fi_vi.readlines()

with open(file_sample_decontamine,"r") as fi_sam:
 	sample_decontamine = fi_sam.readlines()



dico_viral = fct.create_dico_viral(viral)
dico_sample_decontamine = fct.create_presence_dico(sample_decontamine)

res = fct.univ_viral_bis(dico_viral, dico_sample_decontamine).items()

with open('%s/%s_dico_univ.csv' % (path_output, sample_name), 'a', newline='') as csvfile:
				spamwriter = csv.writer(csvfile, delimiter=' ',
		                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
				spamwriter.writerows(res)
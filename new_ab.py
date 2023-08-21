# -*- coding: utf-8 -*-
"""
New table of abundance

    - input: file abundance, file presence, output path
    - output: new table of abundance

example: python new_ab.py file abundance file presence /work_projet/phages/hao/blood_virome/10-Abundance 

"""
import csv
import sys
import functions as fct



file_ab = sys.argv[1]
file_cov = sys.argv[2]

path_output = sys.argv[3]

sample_name = file_cov.split(".coverage")[0]

with open(file_ab,"r") as fi_ab:
	abundance = fi_ab.readlines()



with open(file_cov,"r") as fi_cov:
	coverage = fi_cov.readlines()




presence_dico = fct.create_dico_presence(fct.create_dico_cov(coverage), 0.5)
# precence_dico_diff_0 = 0
# for contig in presence_dico:
# 	if presence_dico[contig] != 0:
# 		precence_dico_diff_0 = precence_dico_diff_0+1
# print(precence_dico_diff_0)



dico_ab = fct.create_dico_ab(abundance)
# dico_ab_diff_0 = 0
# for contig in dico_ab:
# 	if dico_ab[contig] != 0:
# 		dico_ab_diff_0 = dico_ab_diff_0+1
# print(dico_ab_diff_0)

new_dico = fct.new_ab(dico_ab, presence_dico)
# res_dico_diff_0 = 0
# for contig in res_dico:
# 	if res_dico[contig] != 0:
# 		res_dico_diff_0 = res_dico_diff_0+1
# print(res_dico_diff_0)


res_dico = fct.recalcul(new_dico)

res = res_dico.items()
with open('%s/%s_dico_new_ab.csv' % (path_output, sample_name), 'a', newline='') as csvfile:
				spamwriter = csv.writer(csvfile, delimiter=' ',
		                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
				spamwriter.writerows(res)
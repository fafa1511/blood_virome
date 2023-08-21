# -*- coding: utf-8 -*-
"""
Decontamination

    - input: negatif_control file, negatif_control file, presence.csv files
    - output: new present contig decontamined list in csv file

example:decontamination.py /work_projet/phages/hao/blood_virome/07-Presence /work_projet/phages/hao/blood_virome/07-Presence/NegControl2_dico_presence.csv /work_projet/phages/hao/blood_virome/07-Presence/ControlNeg3_dico_presence.csv presence_files output_path
"""
import sys
import csv
import functions as fct

file_ControlNeg2 = sys.argv[1]
file_ControlNeg3 = sys.argv[2]
file_sample = sys.argv[3]
path_out = sys.argv[4]

with open(file_ControlNeg2,"r") as fi_con2:
	ControlNeg2 = fi_con2.readlines()

with open(file_ControlNeg3,"r") as fi_con3:
	ControlNeg3 = fi_con3.readlines()

with open(file_sample,"r") as fi_sam:
	Sample = fi_sam.readlines()

sample_name = file_sample.split("_dico")[0]

dico_ControlNeg23 = fct.ControlNeg23(fct.create_presence_dico(ControlNeg2), fct.create_presence_dico(ControlNeg3))
dico_sample = fct.create_presence_dico(Sample)

fct.decontamination(dico_ControlNeg23, dico_sample)



res = fct.decontamination(dico_ControlNeg23, dico_sample).items()

with open('%s/%s_ConNeg23_dico_presence.csv' % (path_out, sample_name), 'a', newline='') as csvfile:
				spamwriter = csv.writer(csvfile, delimiter=' ',
		                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
				spamwriter.writerows(res)


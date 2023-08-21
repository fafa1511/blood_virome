# -*- coding: utf-8 -*-
"""
Get summry information for all samples

    - input: all presence.csv files, output path
    - output: csv file with summary table

example: python prop.py /work_projet/phages/work_projet/phages/hao/blood_virome/07-Presence
"""
import glob
import csv
import sys

# count contig in different conditions
def prop(li_stool, li_plasma):
    all_1 = []
    s0_p1 = []
    s1_p0 = []
    all_0 = []
    prop_all_1 = 0
    prop_s0_p1 = 0
    prop_s1_p0 = 0
    prop_all_0 = 0
    tot = 0
    #for each pair for contig
    for i in range(len(li_stool)):
        # get contig information
        info_stool = li_stool[i].split()
        info_plasma = li_plasma[i].split()

        # if same order of contig
        if (info_stool[0] == info_plasma[0]):
            sample_name = info_stool[0]
            ab_stool = int(info_stool[1])
            ab_plasma = int(info_plasma[1])

            # all presence or all absence
            if (ab_stool == ab_plasma):
                if (ab_stool == 1):
                    prop_all_1 = prop_all_1+1
                    all_1.append(sample_name)
                    tot = tot+1
                    
                else:
                    prop_all_0 = prop_all_0+1
                    all_0.append(sample_name)
                    

            else:
                # present in stool, absent in plasma
                if (ab_stool > ab_plasma):
                    prop_s1_p0 = prop_s1_p0+1
                    s1_p0.append(sample_name)
                    tot = tot+1
                
                # present in plasma, absent in stool
                else:
                    prop_s0_p1 = prop_s0_p1+1
                    s0_p1.append(sample_name)
                    tot = tot+1
                    
        else:
            print("error of contig's order")

    # return all_1, s0_p1, s1_p0, all_0, [tot, prop_all_1, prop_s0_p1, prop_s1_p0, prop_all_0], [prop_all_1/tot, prop_s0_p1/tot, prop_s1_p0/tot, prop_all_0/tot]
    return str(tot), str(prop_all_1), str(prop_s0_p1), str(prop_s1_p0), str(prop_all_0)



# get path
path = sys.argv[1]

# find all file adundance (stool + plasma)
filex_stool = glob.glob("Stool*")
filex_plasma = glob.glob("Plasma*")
# sort the lists so that they are in the same order
filex_stool.sort()
filex_plasma.sort()

# browse all files two by two
for j in range(len(filex_stool)):
    sample = filex_stool[j].split("_a")[0]+"_"+filex_plasma[j].split("_a")[0]
    with open(filex_stool[j],"r") as fi_stool:
 	    li_stool = fi_stool.readlines()
    with open(filex_plasma[j],"r") as fi_plasma:
 	    li_plasma = fi_plasma.readlines()
         
    # checks if the files two by two have the same number of lines
    if (len(li_stool) != len(li_plasma)):                      # if diffrent number of lines

        if (len(li_stool) > len(li_plasma)):                   # if file stool has more lines
            unique_stool = list(dict.fromkeys(li_stool))       # remove duplicat
            
        else:                                                  # if plasma stool has more lines
            unique_plasma = list(dict.fromkeys(li_plasma))
            
    # if same number of lines we keep the initial list
    else:
        unique_stool = li_stool
        unique_plasma = li_plasma

        # tot = len(unique_stool)
        res = prop(unique_stool, unique_plasma)
        outpout = sample,res
    


    # write csv file
    with open('%s/res.csv' %(path), 'a', newline='') as csvfile:
    				spamwriter = csv.writer(csvfile, delimiter=' ',
    		                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
    				spamwriter.writerows(outpout)

    print("%s done" % (sample))




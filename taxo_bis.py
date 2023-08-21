"""
Create phyloseq objet: table of taxonomie

"""
import csv
import pandas as pd



# read files
with open("/work_projet/phages/hao/blood_virome/07-Counts/presence/test_50/li_presence/ControlNeg23.txt","r") as fi_NegControl:
	contig_contaminant = fi_NegControl.readlines()


with open("/work_projet/phages/Quentin/blood_virome/Univ_viral_contigs.txt","r") as fi_viral:
	contig_viral = fi_viral.readlines()


with open("/home/hzhang/work/work/dico/07-Presence/Plasma_A_R1R2_dico_presence.csv","r") as fi_contig:
	contig = fi_contig.readlines()









contig_commun = ["SPAdes_SampleA_NODE1177","SPAdes_SamplePM_NODE1617","SPAdes_SamplePM_NODE933","SPAdes_SampleST_NODE1913","SPAdes_SampleST_NODE3305","SPAdes_SampleST_NODE977","SPAdes_SampleL_NODE30","Cross_Assembly_MB_NODE3237","Cross_Assembly_MB_NODE8120","SPAdes_SampleL_NODE5485","Cross_Assembly_MB_NODE1712","SPAdes_SampleL_NODE30","SPAdes_SampleBO_NODE6328","Cross_Assembly_MB_NODE31620","SPAdes_SampleBM_NODE1147","SPAdes_SampleBM_NODE1842","SPAdes_SampleBM_NODE4186","SPAdes_SampleBM_NODE922","SPAdes_SampleBT_NODE2391","SPAdes_SampleC_NODE834","SPAdes_SampleMB_NODE709","SPAdes_SampleN_NODE1526","SPAdes_SampleNR_NODE227","SPAdes_SampleST_NODE702","SPAdes_SampleBP_NODE2059","Cross_Assembly_MB_NODE15411","SPAdes_SampleBP_NODE727","SPAdes_SampleMB_NODE234"]

def create_commun(liste_contig):
    dico = {}
    for i in range(len(liste_contig)):
        dico[liste_contig[i]] = 0
    return dico

# get contaminant/viral contig list
def create_dico(liste_contig):
    dico = {}
    for i in range(len(liste_contig)):
        dico[liste_contig[i].rstrip('\n')] = 0
    return dico

# get contig list
def create_presence_dico(input):
    dico = {}

    for i in range(len(input)):
        info = input[i].split()
        # contig name
        contig = info[0]
        # presence value
        presence = int(info[1])

        dico[contig] = presence

    return dico


# assignation contigs
def create_taxo_dico(dico_contig, dico_NegControl, dico_viral, dico_commun):
    dico_taxo = {}
    dico_taxo['contig_contaminant'] = []
    dico_taxo['contig_viral'] = []
    dico_taxo['contig_commun'] = []
    for contig in dico_contig:

        # if contig is found in contaminant list
        if contig in dico_NegControl.keys():
            dico_taxo['contig_contaminant'].append('oui')
        else:
            dico_taxo['contig_contaminant'].append('non')
        
        # if contig is found in viral list
        if contig in dico_viral.keys():
            dico_taxo['contig_viral'].append('oui')
        else:
            dico_taxo['contig_viral'].append('non')

        # if contig is found in commun list
        if contig in dico_commun.keys():
            dico_taxo['contig_commun'].append('oui')
        else:
            dico_taxo['contig_commun'].append('non')

    return dico_taxo


# get contaminant list
dico_NegControl = create_dico(contig_contaminant)

# get viral list
dico_viral = create_dico(contig_viral)

# get contig list
dico_contig = create_presence_dico(contig)

# get contig commun list
dico_commun = create_commun(contig_commun)

# assignation contigs
dico_taxo = create_taxo_dico(dico_contig, dico_NegControl, dico_viral, dico_commun)


# create index for dataframe (contig name)
ind = []
for contig in dico_contig:
    ind.append(contig)

# create dataframe
data = pd.DataFrame(data = dico_taxo, index = ind)

# get csv file for dataframe
data.to_csv('/work_home/hzhang/work/dico/res/taxo_bis.csv', sep = '\t')

"""
Functions
"""

# file coverage -> get dico_coverage
def create_dico_cov(coverage):
    dico_cov = {}
    # all_0 = 0
    # cov_0 = 0
    # prof_0 = 0
    for i in range(len(coverage)):
        info = coverage[i].split()
        # contig name
        contig = info[0]
        # fraction of sequence coverage value
        cov = float(info[1])
        # sequencing-coverage value
        prof = float(info[2])

        # if cov != 0 and prof != 0:
        #     all_0 = all_0+1

        # if cov == 0 and prof != 0:
        #     cov_0 = cov_0+1
        # if prof == 0 and cov != 0:
        #     prof_0 = prof_0+1

        
        dico_cov[contig] = {}
        dico_cov[contig]['cov'] = cov
        dico_cov[contig]['prof'] = prof

    # print("all 0:", all_0)
    # print("cov_0:", cov_0)
    # print("prof_0:", prof_0)

    return dico_cov



# with dico_coverage -> get dico_presence (0/1)
def create_dico_presence(dico_cov, threshold_cov):
    dico_presence = {}

    for contig in dico_cov:
        if (float(dico_cov[contig]['prof']) == 0 or float(dico_cov[contig]['prof']) < 1):
                dico_presence[contig] = 0

        else:
            if (float(dico_cov[contig]['cov']) > threshold_cov):  # > 0.5
                dico_presence[contig] = 1                         # present
                        
            else:                                                 # < 0.5  
                dico_presence[contig] = 0                         # absent

    return dico_presence



# file presence -> get dico_presence
def create_presence_dico(input):
    dico = {}

    for i in range(len(input)):
        info = input[i].split()
        # contig name
        contig = info[0]
        # fraction of sequence coverage value
        presence = int(info[1])

        dico[contig] = presence

    return dico



# get one contig list from two negatif controls
def ControlNeg23(dico_ControlNeg2, dico_ControlNeg3):
    new_dico_ControlNeg23 = {}

    for contig in dico_ControlNeg2:
        if dico_ControlNeg2[contig] == 1 or dico_ControlNeg3[contig] == 1:
            new_dico_ControlNeg23[contig] = 1
        else:
            new_dico_ControlNeg23[contig] = 0

    return new_dico_ControlNeg23




# remove contaminants
def decontamination(dico_ControlNeg23, dico_sample):
    for contig in dico_ControlNeg23:
        if dico_ControlNeg23[contig] == 1:
            dico_sample[contig] = 0

    return dico_sample



# from viral list -> get dico_viral
def create_dico_viral(list_viral):
    dico_viral = {}
    for i in range(len(list_viral)):
        dico_viral[list_viral[i].rstrip('\n')] = 0
    return dico_viral



# keep only viral contig's presence
def univ_viral(dico_viral, dico_sample_decontamine):
    for contig in dico_viral:
        if contig in dico_sample_decontamine:
            if dico_sample_decontamine[contig] == 1:
                dico_viral[contig] = 1
    return dico_viral


def univ_viral_bis(dico_viral, dico_sample_decontamine):
    for contig in dico_sample_decontamine:
        if contig in dico_viral:
            if dico_sample_decontamine[contig] == 1:
                dico_sample_decontamine[contig] = 1
        else:
            dico_sample_decontamine[contig] = 0
    return dico_sample_decontamine


        


def create_dico_ab(abundance):
    dico_ab = {}

    for i in range(9,len(abundance)):
        info_ab = abundance[i].split()
        # contig name
        contig = info_ab[0]
        # abundance relatif
        ab = float(info_ab[1])

        dico_ab[contig] = ab

    return dico_ab



def new_ab(dico_ab, dico_presence):
    # changement = 0
    for contig in dico_ab:
        if (dico_ab[contig] != 0) and (dico_presence[contig] == 0)  :
            dico_ab[contig] = 0
    #         changement = changement+1
    # print("changement:", changement)
    return dico_ab

def recalcul(dico_ab):
    # diff_0 = 0
    ab_tot = sum(dico_ab.values())
    for contig in dico_ab:
        if dico_ab[contig] != 0 :
            # diff_0 = diff_0+1
            dico_ab[contig] = dico_ab[contig]/ab_tot
    # print("nb diff 0: ", diff_0)
    # print("ab total:", ab_tot)
    return dico_ab


def create_new_ab_dico(input):
    dico = {}

    for i in range(len(input)):
        info = input[i].split()
        # contig name
        contig = info[0]
        # fraction of sequence coverage value
        ab = float(info[1])

        dico[contig] = ab

    return dico

# def summary(dico_presence, dico_decontamination, contig_viral):
#     dico_summary = {}
#     for contig in dico_presence:
#         dico_summary[contig] = {}
#         dico_summary[contig]['presence'] = dico_presence[contig]
#         dico_summary[contig]['decontamination'] = dico_decontamination[contig]
#         dico_summary[contig]['contig_viral'] = contig_viral[contig]
#     return dico_summary
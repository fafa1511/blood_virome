##################################### MAPPING #####################################
conda activate bwa-mem2-2.2.1

cd /work_projet/phages/Quentin/blood_virome/01-Non_Human_Samples/Maud_samples_non_human
for file in $(ls *.1.fastq.gz | sed 's/_NonHuman.1.fastq.gz//')
do
    bwa-mem2 mem -a -t 24 /work_projet/phages/Quentin/blood_virome/06-Mapping/Univ_dereplicated ${file}_NonHuman.1.fastq.gz ${file}_NonHuman.2.fastq.gz > /work_projet/phages/hao/blood_virome/06-Mapping/${file}_bwa.sam
    bwa-mem2 mem -a -t 24 /work_projet/phages/Quentin/blood_virome/06-Mapping/Univ_dereplicated ${file}_U_NonHuman.fastq.gz > /work_projet/phages/hao/blood_virome/06-Mapping/${file}_U_bwa.sam
    echo "$file"" done"
done


cd /work_projet/phages/Quentin/blood_virome/01-Non_Human_Samples/Ilias_samples_non_human
for file in $(ls *.1.fastq.gz | sed 's/_NonHuman.1.fastq.gz//')
do
    bwa-mem2 mem -a -t 24 /work_projet/phages/Quentin/blood_virome/06-Mapping/Univ_dereplicated ${file}_NonHuman.1.fastq.gz ${file}_NonHuman.2.fastq.gz > /work_projet/phages/hao/blood_virome/06-Mapping/${file}_bwa.sam
    bwa-mem2 mem -a -t 24 /work_projet/phages/Quentin/blood_virome/06-Mapping/Univ_dereplicated ${file}_U_NonHuman.fastq.gz > /work_projet/phages/hao/blood_virome/06-Mapping/${file}_U_bwa.sam
    echo "$file"" done"
done

conda deactivate

# -a = output all of the aligments for unpaired reads
# -t = number of threads to use

############################# CONVERSION: SAM -> BAM + SORT BY NAME #############################
conda activate samtools-1.9

cd /work_projet/phages/hao/blood_virome/06-Mapping
for file in $(ls *_U_bwa.sam | sed 's/_U_bwa.sam//')
do
    samtools view -bS ${file}_bwa.sam | samtools sort -n  > /work_projet/phages/hao/blood_virome/06-Mapping/bam/${file}_sortedByName_bwa.bam 
    samtools view -bS ${file}_U_bwa.sam | samtools sort -n  > /work_projet/phages/hao/blood_virome/06-Mapping/bam/${file}_U_sortedByName_bwa.bam 
    samtools merge /work_projet/phages/hao/blood_virome/06-Mapping/bam/${file}_sortedByName_bwa_all.bam /work_projet/phages/hao/blood_virome/06-Mapping/bam/${file}_sortedByName_bwa.bam /work_projet/phages/hao/blood_virome/06-Mapping/bam/${file}_U_sortedByName_bwa.bam
    echo "$file"" done"
done

conda deactivate

# view     SAM<->BAM<->CRAM conversion
# -b       output BAM
# -S       ignored (input format is auto-detected)
# sort     sort alignment file
# -n       sort by read name

#################################### MSAMTOOLS ######################################
conda activate msamtools-1.0.0

# filter
cd /work_projet/phages/hao/blood_virome/06-Mapping/bam
for file in $(ls *_bwa_all.bam | sed 's/_sortedByName_bwa_all.bam//')
do
    msamtools filter -b -l 80 -p 95 -z 80 --besthit ${file}_sortedByName_bwa_all.bam > /work_projet/phages/hao/blood_virome/07-Presence/Filter/${file}_bwa_all.filtered.bam
    echo "$file"" done"
done

conda deactivate

# -b                        output BAM (default: false)
# -l <int>                  min. length of alignment (default: 0)
# -p <int>                  min. sequence identity of alignment, in percentage, integer between 0 and 100; requires NM field to be present (default: 0)
# -z <int>                  min. percent of the query that must be aligned, between 0 and 100 (default: 0)
# --besthit                 keep all highest scoring hit(s) per read (default: false)

# profil
cd /work_projet/phages/hao/blood_virome/07-Presence/Filter
for file in $(ls *.bam | sed 's/_bwa_all.filtered.bam//')
do
    msamtools profile --multi=proportional --label=${file} --unit=rel -o /work_projet/phages/hao/blood_virome/07-Presence/Profil/${file}.profile.txt.gz ${file}_bwa_all.filtered.bam
    echo "$file"" done"
done

cd /work_projet/phages/hao/blood_virome/07-Presence/Profil
gzip -dr /work_projet/phages/hao/blood_virome/07-Presence/Profil


#   -o <file>                 name of output file (required)
#   --label=<string>          label to use for the profile; typically the sample id (required)
#   --unit=<string>           unit of abundance to report {ab | rel | fpkm | tpm} (default: rel)
#   --multi=<string>          how to deal with multi-mappers {all | equal | proportional} (default: proportional)


# coverage
cd /work_projet/phages/hao/blood_virome/07-Presence/Filter
for file in $(ls *.bam | sed 's/_bwa_all.filtered.bam//')
do
    msamtools coverage -z --summary -o /work_projet/phages/hao/blood_virome/07-Presence/Coverage/${file}.coverage.txt.gz ${file}_bwa_all.filtered.bam
    echo "$file"" done"
done

cd /work_projet/phages/hao/blood_virome/07-Presence/Coverage
gzip -dr /work_projet/phages/hao/blood_virome/07-Presence/Coverage

conda deactivate

#   -o <file>                 name of output file (required)
#   --summary                 do not report per-position coverage but report fraction of sequence covered (default: false)
#   -z, --gzip                compress output file using gzip (default: true)


###################################### PRESENCE ######################################
# Determine presence or absence of contigs
# - input: all covrage.txt file, values for threshold, output path
# - output: all presence.csv file


cd /work_projet/phages/hao/blood_virome/07-Presence/Coverage
for file in $(ls *.txt | sed 's/.coverage.txt//')
do
    python /work_projet/phages/hao/blood_virome/scripts/presence_dico.py ${file}.coverage.txt  0.5 /work_projet/phages/hao/blood_virome/07-Presence
    echo "$file"" done"
done


###################################### PROPORTION #####################################
# Get summry information for all samples (total present: all present/present only in plasma/present only in stool/all absent)
# - input: all presence.csv files, output path
# - output: csv file with summary table

cd /work_projet/phages/hao/blood_virome/07-Presence
python /work_projet/phages/hao/blood_virome/scripts/prop.py /work_projet/phages/hao/blood_virome/07-Presence

# verifie avec awk
cd /work_projet/phages/hao/blood_virome/07-Presence
for file in $(ls *.csv | sed 's/_dico_presence.csv//')
do
    awk -F " " '{CUMUL += $2} END {print CUMUL}' ${file}_dico_presence.csv
    echo "$file"" done"
done

###################################### DISTRIBUTION #####################################
# look the distribution of contig length
# - input: presence file, length file, output path
# - output: csv file with contig length

python /work_projet/phages/hao/blood_virome/scripts/distribution.py Plasma_A_R1R2_dico_presence.csv /work_projet/phages/Quentin/blood_virome/02-Contig_Clustering_Blat/Univ_dereplicated_lengths.tab /work_projet/phages/hao/blood_virome/07-Presence/Distribution
python /work_projet/phages/hao/blood_virome/scripts/distribution.py Stool_A_dico_presence.csv /work_projet/phages/Quentin/blood_virome/02-Contig_Clustering_Blat/Univ_dereplicated_lengths.tab /work_projet/phages/hao/blood_virome/07-Presence/Distribution

python /work_projet/phages/hao/blood_virome/scripts/distribution.py PlasmaBK_dico_presence.csv /work_projet/phages/Quentin/blood_virome/02-Contig_Clustering_Blat/Univ_dereplicated_lengths.tab /work_projet/phages/hao/blood_virome/07-Presence/Distribution
python /work_projet/phages/hao/blood_virome/scripts/distribution.py StoolBK_dico_presence.csv /work_projet/phages/Quentin/blood_virome/02-Contig_Clustering_Blat/Univ_dereplicated_lengths.tab /work_projet/phages/hao/blood_virome/07-Presence/Distribution

################################## DECONTAMINATION ####################################
# remove contaminant contigs
#  - input: negatif_control file, negatif_control file, all presence.csv files, output path
#  - output: new present contig decontamined list in csv file

cd /work_projet/phages/hao/blood_virome/07-Presence
for file in $(ls *.csv | sed 's/_dico_presence.csv//')
do
    python /work_projet/phages/hao/blood_virome/scripts/decontamination_dico.py /work_projet/phages/hao/blood_virome/07-Presence/NegControl2_dico_presence.csv /work_projet/phages/hao/blood_virome/07-Presence/ControlNeg3_dico_presence.csv ${file}_dico_presence.csv  /work_projet/phages/hao/blood_virome/08-Decontamination
    echo "$file"" done"
done


###################################### PROPORTION #####################################
cd /work_projet/phages/hao/blood_virome/08-Decontamination
python /work_projet/phages/hao/blood_virome/scripts/prop.py /work_projet/phages/hao/blood_virome/08-Decontamination


# verifie avec awk
for file in $(ls *.csv | sed 's/_ConNeg23_dico_presence.csv//')
do
    awk -F " " '{CUMUL += $2} END {print CUMUL}' ${file}_ConNeg23_dico_presence.csv
    echo "$file"" done"
done

############################### KEEP ONLY VIRAL CONTIG ##############################
# Keep only viral contigs with presence or absence status
# - input: presence.csv files after decontamination, viral list, output path
# - output: presence.csv files after decontamination with only viral contig
# liste des contig viraux :Univ_viral_contigs.txt 


cd /work_projet/phages/hao/blood_virome/08-Decontamination
for file in $(ls *.csv | sed 's/_ConNeg23_dico_presence.csv//')
do
    python /work_projet/phages/hao/blood_virome/scripts/univ_viral_dico.py /work_projet/phages/Quentin/blood_virome/Univ_viral_contigs.txt ${file}_ConNeg23_dico_presence.csv /work_projet/phages/hao/blood_virome/09-Viral
    echo "$file"" done"
done


###################################### PROPORTION #####################################

cd /work_projet/phages/hao/blood_virome/09-Viral
python /work_projet/phages/hao/blood_virome/scripts/prop.py /work_projet/phages/hao/blood_virome/09-Viral


# verifie avec awk
for file in $(ls *.csv | sed 's/_dico_univ.csv//')
do
    awk -F " " '{CUMUL += $2} END {print CUMUL}' ${file}_dico_univ.csv
    echo "$file"" done"
done

##################################### STATISTIQUE #####################################
# script: stat.R

##################################### ABUNDANCE #####################################

######### abundance brut
# get ew table of abundance
# - input: file abundance, file presence
# - output: new table of abundance

cd /work_projet/phages/hao/blood_virome/07-Presence/Coverage

for file in $(ls *.txt | sed 's/.coverage.txt//')
do
    python /work_projet/phages/hao/blood_virome/scripts/new_ab.py /work_projet/phages/hao/blood_virome/07-Presence/Profil/${file}.profile.txt ${file}.coverage.txt /work_projet/phages/hao/blood_virome/10-Abundance
    echo "$file"" done"
done


######### dataframe
# get dataframe for new tables of abundance
cd /work_projet/phages/hao/blood_virome/10-Abundance

conda activate pandas-1.5.3

python /work_projet/phages/hao/blood_virome/scripts/data_ab.py

conda deactivate


############################## LIST OF PRESENT CONTIG ##############################
# get lists of presence contig for 2 controlNg

cd /work_projet/phages/hao/blood_virome/07-Presence
for file in $(ls *Control* | sed 's/_dico_presence.csv//')
do
    python /work_projet/phages/hao/blood_virome/scripts/list.py ${file}_dico_presence.csv /work_projet/phages/hao/blood_virome/07-Presence/list_contig
    echo "$file"" done"
done


######### ONE LIST FOR NEGATIF CONTROL 
cd /work_projet/phages/hao/blood_virome/07-Presence/list_contig
cat NegControl2_dico_presence.csv_li.csv ControlNeg3_dico_presence.csv_li.csv | sort -u > ControlNeg23.txt

######### create phyloseq objet: tax_table() in Rstudio
conda activate pandas-1.5.3

python /work_projet/phages/hao/blood_virome/scripts/taxo_bis.py

conda deactivate

############################## ANALYZE CLASSIC: PHYLOSEQ ##############################
# R: phyloseq.R

######### get list of commun contigs
# Dirived from prop.py to et list of commun contigs for all samples

#- input: all presence.csv files, output path
#- output: csv file with list of commun contigs in plasma and in stool

cd /work_projet/phages/hao/blood_virome/09-Viral
python /work_projet/phages/hao/blood_virome/scripts/prop_contig_commun.py /work_projet/phages/hao/blood_virome/09-Viral


######## ABUNDANCE OF COMMUN CONTIGS
# R: ab_contigs_commun_all.R


############################## ANALYZE SAMPLE PLASMA ##############################
# R: Analyze_plasma.R

# dataframe of viral contig presence
cd /work_projet/phages/hao/blood_virome/09-Viral

conda activate pandas-1.5.3

python /work_projet/phages/hao/blood_virome/scripts/data_viral_presence.py

conda deactivate



########################## CHECK SEQUENCING-COVERAGE (profondeur)  ##########################
samtools sort Stool_M_bwa_all.filtered.bam > Stool_M_bwa_all.bam
samtools index Stool_M_bwa_all.bam
samtools view -c Stool_M_bwa_all.bam SPAdes_SampleA_NODE422



############################## GET SUBSAMPLE OF FASTA FILES ##############################
# fasta_reader_select_P3.py: this script is developed by a member of the PAGES team
python fasta_reader_select_P3.py Univ_dereplicated.fasta contig_selected genome_selected_new.fasta


############################## ANNOTATION DES CONTIGS ##############################
conda activate pharokka-1.3.2

pharokka.py -i genome_selected.fasta -o pharokka -d /db/outils/pharokka-1.3.2 -t 16 -p plasma -g prodigal


#   -h, --help            show this help message and exit
#   -i INFILE, --infile INFILE
#                         Input genome file in fasta format.
#   -o OUTDIR, --outdir OUTDIR
#                         Directory to write the output to.
#   -d DATABASE, --database DATABASE
#                         Database directory. If the databases have been installed in the default directory, this is not required. Otherwise specify the path.
#   -t THREADS, --threads THREADS
#                         Number of threads. Defaults to 1.
#   -p PREFIX, --prefix PREFIX
#                         Prefix for output files. This is not required.
#   -g GENE_PREDICTOR, --gene_predictor GENE_PREDICTOR
#                         User specified gene predictor. Use "-g phanotate" or "-g prodigal". Defaults to phanotate (not required unless prodigal is desired).


# get the subgroups
grep -wf partage_new.txt plasma.txt | awk ' { print $1, $9 $10 $11 $12 $13 } ' OFS=":" > annot_partage_new.txt


############################## HOST PREDICTION ##############################
# iPHoP v1.2.0: integrating Host Phage Predictions
# https://bitbucket.org/srouxjgi/iphop

# usage: iphop <task> [options]

# task:
#         predict         run full pipeline to generate host prediction for some input phage genome(s)
#         download        download and setup the latest host database
#         add_to_db       add some host genomes to the database of hosts (e.g. MAGs derived from the same metagenome(s))

# other:
#         split           small utility to split a fasta file into smaller batches (can be useful if you would like to split your input and run iPHoP separately for each batch)
#         clean           small utility to clean the output directory of iPHoP by compressing some of the larger files.


# the commands below are developed by a member of the PHAGES team
conda activate iphop-1.2.0

iphop predict --fa_file genome_selected.fasta \
--db_dir /db/outils/iphop-1.2.0/Sept_2021_pub/ \
--out_dir iphop/ \
--num_threads 24

conda deactivate
# Only keep the best hit for each contig
awk -F ',' 'NR==1 {print; next} {if ($5 > max[$1]) {max[$1]=$5; line[$1]=$0}} END {for (i in line) print line[i]}' Host_prediction_to_genome_m90.csv > Host_prediction_to_genome_m90_besthit.csv
awk -F ',' 'NR==1 {print; next} {if ($4 > max[$1]) {max[$1]=$4; line[$1]=$0}} END {for (i in line) print line[i]}' Host_prediction_to_genus_m90.csv > Host_prediction_to_genus_m90_besthit.csv

# -F ',' = use "," as separator to read the file
# NR==1 {print; next} = print the first line (header) and go to the next line
# {if ($5 > max[$1]) {max[$1]=$5; line[$1]=$0}} = store the line which has the highest value in the 5th column when multiple lines have the same value in the 1st column
# END {for (i in line) print line[i]} = print the stored lines

# Create similar tables for both "Genus" and "Genome" files
awk -F ',' '{OFS=","; print $1, $3, $4}' Host_prediction_to_genus_m90_besthit.csv > summary_genus_m90_besthit.csv
awk -F ',' '{OFS=","; print $1, $3, $5}' Host_prediction_to_genome_m90_besthit.csv > summary_genome_m90_besthit.csv

# -F ',' = use "," as a separator to read the file
# OFS="," = use "," as a separator to write the output
# print $1, $3, $4 / print $1, $3, $5 = only print the selected columns

# Add a column to indicate from which method they come from
awk -F ',' '{OFS=","; print $0, "Genus"}' summary_genus_m90_besthit.csv > tmp_file.txt && mv tmp_file.txt summary_genus_m90_besthit.csv
awk -F ',' ' {OFS=","; print $0, "Genome"}' summary_genome_m90_besthit.csv > tmp_file.txt && mv tmp_file.txt summary_genome_m90_besthit.csv

# -F ',' = use "," as a separator to read the file
# OFS="," = use "," as a separator to write the output
# print $0, "Genus" / print $0, "Genome" = write "Genus" or "Genome" at the end of the line

# Merge both "genus" and "genome" file
cat summary_genus_m90_besthit.csv summary_genome_m90_besthit.csv > merge_genus_genomes.csv

# In case of duplicates, prefer the "Genus" output
awk -F ',' '!seen[$1]++' merge_genus_genomes.csv > final_host_prediction.csv

# -F ',' = use "," as a separator to read the file
# !seen[$1]++ = only print the lines for which the value in the 1st column is encountered for the first time (this only works if you generated the file in the order "genus" - "genome", so that "genus" line are read first by awk going down the file)

# get the subgroups
grep -wf /work_projet/phages/hao/blood_virome/12-Annotation/partage_new.txt final_host_prediction.csv | awk ' { print $1 } ' > final_host_prediction_new.csv

############################## regarde qualite
cd /work_projet/phages/Quentin/blood_virome/04-Viral_Identification_VirVibCheckV/checkv


# verification nb de lignes avant sort
awk ' { print $1 } ' transfert.txt | wc -l

# enlever redondance
sort -u transfert.txt > transfert_u.txt

# avoir file quality
grep -wf abondant_u.txt /work_projet/phages/Quentin/blood_virome/04-Viral_Identification_VirVibCheckV/checkv/quality_summary.tsv | awk ' { print $8 } ' > abondant_quality.txt

# concatenation fichier contig et fichier quality
paste <(echo "$(awk ' {print $1}' abondant_u.txt)") <(echo "$(awk ' {print $1}' abondant_quality.txt)") > abondant_check.txt

for file in $(ls *.txt | sed 's/.txt//')
do
    awk ' { print $1 } ' ${file}.txt | wc -l
    sort -u ${file}.txt > ${file}_u.txt
    awk ' { print $1 } ' ${file}_u.txt | wc -l
    grep -wf ${file}_u.txt /work_projet/phages/Quentin/blood_virome/04-Viral_Identification_VirVibCheckV/checkv/quality_summary.tsv | awk ' { print $8 } ' > ${file}_quality.txt
    paste <(echo "$(awk ' {print $1}' ${file}_u.txt)") <(echo "$(awk ' {print $1}' ${file}_quality.txt)") > ${file}_check.txt
    echo "$file"" done"
done

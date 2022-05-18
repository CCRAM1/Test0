import pandas as pd
import scipy.stats as stats
import statsmodels.stats.api as sms
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#load into dataframes
bvirus = pd.read_csv('all.bigtable.virus.tsv.gz', compression='gzip', header=0, sep='\t')
amrf = pd.read_csv('all.amrfinder.tsv.gz', compression='gzip', header=0, sep='\t')
sfocus  = pd.read_csv('all.superfocus.tsv.gz', compression='gzip', header=0, sep='\t')
focus = pd.read_csv('all.focus.csv.gz', compression= 'gzip', header=None, sep=',')

#bvirus - contains only VIRUSES filtered from bigtable
#filter
bvirusfilt = bvirus[(bvirus.alnType=='aa') & (bvirus.evalue<1e-20)]
bvirusfiltgroup = bvirusfilt.groupby(by=['sampleID', 'family'], as_index=False).agg('sum','normCount')
#selecting columns
fbvirusfiltgroup = bvirusfiltgroup[['sampleID', 'family', 'CPM']]
len(fbvirusfiltgroup.sampleID.unique())


#focus - taxonomic ID and quantification of BACTERIAL sequences
focus.columns = ['sampleID', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain', 'forward', 'backward']
# forward % -> CPM (x a 10,000)
focus["CPM"] = focus["forward"] * 10000
#Need to filter for BACTERIAL sequences - take out archaea
focusfilt = focus[(focus.kingdom== 'Bacteria')]
#want sampleID, family, CPM
focusfiltgroup = focusfilt[['sampleID', 'family', 'genus', 'CPM']]

#sfocus - taxonomy via metabolic functions
#want sampleID, Subsystem Level 3, Function
sfocusgroup = sfocus[['sampleID', 'Subsystem Level 3', 'Function']]

#amrf - antimicrobial genes
#contig id -> sampleID
amrf. rename(columns = {'Contig id':'sampleID'}, inplace = True)
#want sampleID, Sequence name, Element type, Class, Subclass
amrfgroup = amrf[['sampleID', 'Sequence name', 'Element type','Class', 'Subclass']]

#ISSUE - for some reason python console wont execute..
# #Filter for AMR sequences
amrffilt = amrf[(amrf.Element type== 'AMR')]

concatenated = pandas.concat([fbvirusfiltgroup, amrfgroup])
concatenated = pandas.concat([fbvirusfiltgroup, amrfgroup], axis="columns")

fbvirusfiltgroup, focusfiltgroup, sfocusgroup, amrfgroup
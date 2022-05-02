import pandas as pd
import scipy.stats as stats
import statsmodels.stats.api as sms
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#load data into dataframes
bvirus = pd.read_csv('all.bigtable.virus.tsv.gz', compression='gzip', header=0, sep='\t')
amrf = pd.read_csv('all.amrfinder.tsv.gz', compression='gzip', header=0, sep='\t')
sfocus  = pd.read_csv('all.superfocus.tsv.gz', compression='gzip', header=5, sep='\t')
*need to convert forward (%) to CPM (x a million)
focus = pd.read_csv('all.focus.csv.gz', compression= 'gzip', header=0, sep=',')

#bvirus - contains only VIRUSES filtered from bigtable
#filter
bvirusfilt = bvirus[(bvirus.alnType=='aa') & (bvirus.evalue<1e-20)]
bvirusfiltgroup = bvirusfilt.groupby(by=['sampleID', 'family'], as_index=False).agg('sum','normCount')
#selecting columns
bvirusfiltgroup = virusesGroup[['sampleID', 'family', 'normCount']]
#rename nornmCount to CPM

#focus - taxonomic ID and quantification of BACTERIAL sequences
#Need to filter for BACTERIAL sequences - take out archaea
#PROBLEMS - sampleID (DDR form), no column titles,
# % -> CPM
#want sampleID, family, CPM

#sfocus - taxonomy via metabolic functions
#PROBLEMS column and data alignment (esp first 5 rows) - make header row 5 (skip first 5)
#want sampleID, Subsystem Level 3, Function

#amrf - antimicrobial genes
#contig id -> sampleID
#want sampleID, Gene symbol (or Sequence name?), Element type (AMR), Class (or Subclass?)
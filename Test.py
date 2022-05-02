import pandas as pd
import scipy.stats as stats
import statsmodels.stats.api as sms
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
meta = pd.read_csv('metadata.tsv.gz', compression='gzip', header=0, sep='\t')
data = pd.read_csv('bigtable.tsv.gz', compression='gzip', header=0, sep='\t')
taxonCounts = pd.read_csv('taxonLevelCounts.tsv.gz', compression='gzip', header=0, sep='\t')
dataMeta = pd.merge(data, meta, on=["sampleID"])
#viralHits
viruses = dataMeta[(dataMeta.kingdom == "Viruses")]
#filter
virusesGroup = viruses.groupby(by=['family', 'alnType', 'alnlen', 'pident'], as_index=False).count()
#styling
sizeScatter = 10 * virusesGroup['count']
sns.set_style("darkgrid")
sns.set_palette("colorblind")
sns.set(rc={'figure.figsize':(12, 8)})
#createFacetGrid
g = sns.FacetGrid(virusesGroup, col="family", hue="alnType", col_wrap=6)
#createPlot
g.map(sns.scatterplot, "alnlen", "pident", alpha=.1, sizes=(100,500), size=sizeScatter)
plt.legend(bbox_to_anchor=(5.0,1), loc=0, borderaxespad=2, ncol=6, shadow=True, labelspacing=1.5, borderpad=1.5)
plt.show()
#divideIntoQuadrants
g.map(sns.scatterplot, "alnlen", "pident", alpha=.1, sizes=(100,500), size=sizeScatter)
for ax in g.axes.flat:
    ax.tick_params(axis='both', labelleft=True, labelbottom=True)
    ax.axhline(y=75, c='red', linestyle='dashed')
    ax.axvline(x=150, c='red', linestyle='dashed')
plt.legend(bbox_to_anchor=(5.0,1), loc=0, borderaxespad=2, ncol=6, shadow=True, labelspacing=1.5, borderpad=1.5)
plt.show()

#Challenges
#bacterialHits
bacterial = dataMeta[(dataMeta.kingdom == "Bacteria")]
#filter
bacterialGroup = bacterial.groupby(by=['phylum', 'alnType', 'alnlen', 'pident'], as_index=False).count()
#styling
sizeScatter = 10 * bacterialGroup['count']
sns.set_style("darkgrid")
sns.set_palette("colorblind")
sns.set(rc={'figure.figsize':(12, 8)})
#createFacetGrid
g = sns.FacetGrid(bacterialGroup, col="phylum", hue="alnType", col_wrap=6)
g.map(sns.scatterplot, "alnlen", "pident", alpha=.1, sizes=(100,500), size=sizeScatter)
plt.legend(bbox_to_anchor=(5.0,1), loc=0, borderaxespad=2,ncol=6, shadow=True, labelspacing=1.5, borderpad=1.5)
plt.show()

#AdenoviridaeHits
Adenoviridae = viruses[(viruses.family=='Adenoviridae')]
#filter
AdenoviridaeGroup = Adenoviridae.groupby(by=['sampleID', 'alnType', 'alnlen', 'pident'], as_index=False).count()
#styling
sizeScatter = 10 * AdenoviridaeGroup['count']
sns.set_style("darkgrid")
sns.set_palette("colorblind")
sns.set(rc={'figure.figsize':(12, 8)})
#createFacetGrid
g = sns.FacetGrid(AdenoviridaeGroup, col="sampleID", hue="alnType", col_wrap=6)
#createPlot
g.map(sns.scatterplot, "alnlen", "pident", alpha=.1, sizes=(100,500), size=sizeScatter)
plt.legend(bbox_to_anchor=(5.0,1), loc=0, borderaxespad=2, ncol=6, shadow=True, labelspacing=1.5, borderpad=1.5)
plt.show()

#filteringStrategies
#filtering via evalue
# assign() will add or modify columns, np.where() will return a value base on a condition
viruses = viruses.assign(filter = np.where(viruses.evalue<1e-20,'pass','filter'))
#plot
virusesGroup = viruses.groupby(by=['family','alnType','alnlen','pident','filter'], as_index=False).count()
g = sns.FacetGrid(virusesGroup, col="family", hue="filter", col_wrap=6)
g.map(sns.scatterplot, "alnlen", "pident", alpha=.1)
plt.legend(bbox_to_anchor=(5.0,1), loc=0, borderaxespad=2,ncol=6, shadow=True, labelspacing=1.5, borderpad=1.5)
plt.show()
#filtering via alnlen and pident
# this will overwrite the flags with the new designations
viruses = viruses.assign(filter = np.where((viruses.alnlen>150) & (viruses.pident>75),'pass','filter'))
#plot Facet Grid with cut off lines

#Challenge - Filter raw viral hits to keep protein hits with an evalue of <1e-10
virusesFiltered = viruses[viruses.evalue<1e-10]
virusesFiltered = virusesFiltered.assign(filter = np.where((virusesFiltered.alnType 'aa'),'pass','filter'))

#visualising Annotations
# get adenoviridae counts
adenoCounts = taxonCounts[(taxonCounts.taxonLevel=='family') & (taxonCounts.taxonName=='Adenoviridae')]
# plot
sns.set_style("darkgrid")
sns.set_palette("colorblind")
sns.set(rc={'figure.figsize':(14,10)})
sns.barplot(x="count", y="sampleID", data=adenoCounts)
plt.subplots_adjust(left=0.2)
plt.grid(True)
plt.show()

# get all viral family counts
viralCounts = taxonCounts[(taxonCounts.taxonLevel=='family') & (taxonCounts.taxonPath.str.contains('k_Viruses'))]
# plot
sns.set_style("darkgrid")
sns.set(rc={'figure.figsize':(18,14)})
colors = plt.cm.nipy_spectral(np.linspace(0, 1, 50))
viralCountsChartPivot=pd.pivot_table(viralCounts, index=['sampleID'], columns=['taxonName'], values=['count'], aggfunc='sum')
plt.figure();
viralCountsChartPivot.plot.barh(stacked=True,color=colors);
plt.subplots_adjust(left=0.18)
# Shink current axis by 50% to allow legend to fit nicely
ax = plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.5, box.height])
plt.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), borderaxespad=1,ncol=2, shadow=True, labelspacing=1.5, borderpad=1.5)
plt.show()

# Heatmap
# get all viral family counts
viralFiltCounts = viralCounts.groupby(by=['sampleID','taxonName'], as_index=False)['count'].agg('sum')
#plot
sns.set_style("darkgrid")
sns.set(rc={'figure.figsize':(20,10)})
sns.heatmap(pd.crosstab([viralFiltCounts.sampleID], [viralFiltCounts.taxonName], values=viralFiltCounts['count'], aggfunc='sum', dropna=False).fillna(0),
            cmap="YlGnBu", annot=True, cbar=False, fmt=".0f")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

# Bubble plot
# get all viral family counts
viralFiltCounts = viralCounts.groupby(by=['sampleID','taxonName'], as_index=False)['count'].agg('sum')
#plot
sns.set_style("darkgrid")
sns.set(rc={'figure.figsize':(16,12)})
sns.scatterplot(x="taxonName", y="sampleID", data=viralFiltCounts, hue="count", s=viralFiltCounts['count'])
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

# Generating taxon counts
# Answer for "Challenge: Filter your raw viral hits to only keep protein hits with an evalue < 1e-10"
virusesFiltered = viruses[(viruses.alnType=='aa') & (viruses.evalue<1e-10)]
# Heatmap
sns.set_style("darkgrid")
sns.set(rc={'figure.figsize':(20,10)})
sns.heatmap(pd.crosstab([virusesFiltered.sampleID], [virusesFiltered.family], values=virusesFiltered.normCount, aggfunc='sum', dropna=False).fillna(0),
            cmap="YlGnBu", annot=True, cbar=False, fmt=".0f")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
# Bubble plot
sns.set_style("darkgrid")
sns.set(rc={'figure.figsize':(16,12)})
sns.scatterplot(x="family", y="sampleID", data=virusesFiltered, hue="normCount", s=virusesFiltered.normCount)
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

# Challenge - stacked bar chart, heatmap or bubble plot of the viral families for the Male and Female monkeys
# Heatmap
sns.set_style("darkgrid")
sns.set(rc={'figure.figsize':(20,10)})
sns.heatmap(pd.crosstab([virusesFiltered.sex], [virusesFiltered.family], values=virusesFiltered.normCount, aggfunc='sum', dropna=False).fillna(0),
            cmap="YlGnBu", annot=True, cbar=False, fmt=".0f")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
# Bubble plot
sns.set_style("darkgrid")
sns.set(rc={'figure.figsize':(16,5)})
sns.scatterplot(x="family", y="sex", data=virusesFiltered, hue="normCount", s=virusesFiltered.normCount)
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()
# stacked bar chart
print('hello')


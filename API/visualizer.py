import sys
sys.path.insert(1, r'C:\Users\aprak\PycharmProjects\TopProPack_v2_2\API')

import matplotlib as plt
import seaborn as sns
from collections import Counter
import API.centroid_protein as centroid_protein
import streamlit
import matplotlib.pyplot as plt
import pandas as pd



def draw_countplot(data, x, title):
	df = pd.DataFrame()
	df[x] = data
	ax = sns.countplot(x=x, data=df, order=list(range(max(data)+1)))
	plt.title("Clique Sizes Histogram")
	for p in ax.patches:
		ax.annotate('{}'.format(p.get_height()), (p.get_x() + 0.1, p.get_height()))
	plt.show()

def draw_histogram(data, title, normalized=False):
	ax = None
	if not normalized:
		ax = sns.distplot(data, kde=False, norm_hist=False)
	else:
		ax = sns.distplot(data)
	plt.show()

def draw_heatmap(name, heatmap_data, x_labels, y_labels, cmap, center=0):
	#maximum = 0
	#minimum = 0
	#for i in heatmap_data:
	#	for j in i:
	#		if j > maximum: maximum = j
	#		if j < minimum: minimum = j
	#if center is None: plot = sns.heatmap(heatmap_data, xticklabels=x_labels, yticklabels=y_labels, robust=True, cmap=cmap, vmax=maximum, vmin=minimum)
	plot = sns.heatmap(heatmap_data, xticklabels=x_labels, yticklabels=y_labels, center=0, vmin=-1, vmax=1, robust=True, cmap=cmap)
	plt.savefig("{}.png".format(name))
	plt.clf()
	#plt.show()

#draw_countplot([len(s) for s in centroid_protein.P.centroid_cliques], "clique_sizes", "Clique Size Frequencies")
#draw_histogram(centroid_protein.P.centroid_clique_distances, "Centroid Clique Distances", normalized=True)

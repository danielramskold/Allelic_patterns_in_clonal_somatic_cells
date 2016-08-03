from __future__ import division
import argparse, dr_tools, numpy, hcluster, random
import matplotlib.pyplot as pylab
import scipy.cluster.hierarchy as scipyhcluster

def Xi_activity_similarity(X, Y):
	num_different = sum(x != y for x,y in zip(X,Y))
	possibly_different = sum(X)+sum(Y)
	return num_different/possibly_different # corresponds to the binary distance i the R function dist

stateD = {'XI':1, 'bi':1, 'nd':0, 'xa':0}

if '__main__' == __name__:
	
	# load table
	linefeed = dr_tools.splitlines('chrX_clones_allelic_calls.txt')
	sample_labels = next(linefeed)[1:]
	character_matrixT = []
	for cells in linefeed:
		# values in cells are nd, XI, xa, bi, except first column which is gene symbol
		if any(c!='nd' for c in cells):
			character_matrixT.append([stateD[c] for c in cells[1:]])
	
	# make clusters
	character_matrix = numpy.array(character_matrixT).transpose()
	#hcdists = hcluster.pdist(character_matrix, metric='cityblock')
	hcdists = hcluster.pdist(character_matrix, metric=Xi_activity_similarity)
	hclinks = hcluster.linkage(hcdists, method='complete')
	draw_order = hcluster.leaves_list(hclinks)
	
	# draw tree
	scipyhcluster.dendrogram(hclinks, labels=sample_labels, leaf_rotation=90)
	pylab.subplots_adjust(bottom=0.3)
	pylab.savefig('tree_Xiexpr.pdf')

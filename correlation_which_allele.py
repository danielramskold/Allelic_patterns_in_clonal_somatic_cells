from __future__ import division
import argparse, dr_tools, numpy, hcluster, random
import matplotlib.pyplot as pylab
import scipy.cluster.hierarchy as scipyhcluster


def expandlinks(links, samples):
	# find which numbers are which groups (e.g. 6 for [1, 2] gives dictionary entry [6, 1, 2])
	dictionary = []
	highest_id = samples-1
	for link in links:
		highest_id += 1
		dictionary.append([highest_id] + link)
		
	# expand them (e.g. [6, 3] becomes [1, 2, 3] if 6 is for [1, 2])
	hasunexpanded = 1
	while hasunexpanded:
		hasunexpanded = 0
		for link in links:
			for element in range(len(link)):
				if link[element] >= samples:
					hasunexpanded = 1
					for entry in dictionary:
						if entry[0] == link[element]:		
							link.extend(entry[1:])	# add the values from dictionary
							break
					link.pop(element)	# remove old value
	for link in links:
		link.sort()

def gethclinks(exparray, method):
	hcdists = hcluster.pdist(exparray, method)
	hclinks = hcluster.linkage(hcdists)
	links = []
	for hclink in hclinks:
		links.append([int(hclink[0]), int(hclink[1])])
	return links
	
def generatebootstraparray(exparray):
	arraylength = len(exparray)
	arraywidth = len(exparray[0])
	newarray = []
	for row in range(arraylength):
		newarray.append([])
	for column in range(arraywidth):
		chosencolumn = random.randrange(0, arraywidth)
		for row in range(arraylength):
			newarray[row].append(exparray[row][chosencolumn])
	return newarray
	
class Clinkpattern:
	def __init__(self, linked):
		self.linked = linked	# e.g. [1, 2, 7]
		self.number = 0
	def matcheslinked(self, linked):
		if len(linked) != len(self.linked):
			return 0
		for ii in range(len(linked)):
			if linked[ii] != self.linked[ii]:
				return 0
		return 1
	def __str__(self):
		outstr = str(self.number) + " ["
		for l in self.linked:
			outstr += str(l) + " "
		outstr = outstr[:-1] + "]"
		return outstr
		
	def getstring(self, tissuenames):
		outstr = str(self.number) + " ["
		for l in self.linked:
			outstr += str(tissuenames[l]) + " "
		outstr = outstr[:-1] + "]"
		return outstr	
	
def bootstraphcluster(exparray, number, seed, method, originalhclink=None):
	# call this function to get bootstrap values
	# exparray is 2D array of values
	# number is how many selections should be generated
	# seed is an arbitrary number
	
	
	# generate the nodes from random selections from exparray
	random.seed(seed)
	all_links = []
	for ii in range(number):
		newlink = gethclinks(generatebootstraparray(exparray), method)
		expandlinks(newlink, len(exparray))
		all_links.append(newlink)
		
	# select the patterns in exparray
	patterns = []
	if originalhclink is None:
		originallink = gethclinks(exparray, method)
	else:
		hclinks = originalhclink
		links = []
		for hclink in hclinks:
			links.append([int(hclink[0]), int(hclink[1])])
		originallink = links
	expandlinks(originallink, len(exparray))
	for linkitem in originallink:
		patterns.append(Clinkpattern(linkitem))
		
	# find bootstrap values
	for links in all_links:
		for linked in links:
			for pattern in patterns:
				if pattern.matcheslinked(linked):
					pattern.number += 1
					break	
	return patterns

def state(maternal, paternal, o):
	if o.states == '3state':
		if (maternal>0, paternal>0) == (True,False): return 1
		elif (maternal>0, paternal>0) == (False,True): return -1
		else:
			return 0
			#return 0.67*(2*random.random()-1) # so that genes that express fewer genes shouldn't cluster better, or at least reduce the bias (weight factor uncertain)
	elif o.states == 'fraction':
		if maternal or paternal:
			return paternal/(paternal+maternal)
		else:
			return 0.5
	elif o.states == 'diff':
		return maternal-paternal
	elif o.states == 'monoallelic': #same as 3state
		return 0 if maternal>0 and paternal>0 else 1 if paternal>0 else -1 if maternal>0 else 0

def monoallelic_norm(X, Y):
	X = list(X)
	Y = list(Y)
	Xpos = X.count(1)
	Xneg = X.count(-1)
	Ypos = Y.count(1)
	Yneg = Y.count(-1)
	if Xpos+Xneg==0 or Ypos+Yneg==0: return 1
	Xposf = Xpos/(Xpos+Xneg)
	Yposf = Ypos/(Ypos+Yneg)
	Xnegf = Xneg/(Xpos+Xneg)
	Ynegf = Yneg/(Ypos+Yneg)
	scores = {(1,1): 0.5+0.125/Xposf/Yposf, (1,-1): 0.5-0.125/Xposf/Ynegf, (-1,1): 0.5-0.125/Xnegf/Yposf, (-1,-1): 0.5+0.125/Xnegf/Ynegf}
	return max(0,1-numpy.mean([scores[xi,yi] for xi,yi in zip(X,Y) if xi and yi]))

def monoallelic_norm2(X, Y):
	X = list(X)
	Y = list(Y)
	Xpos = X.count(1)
	Xneg = X.count(-1)
	Ypos = Y.count(1)
	Yneg = Y.count(-1)
	if Xpos+Xneg==0 or Ypos+Yneg==0: return 0
	raw_overlap = 1-monoallelic_raw(X, Y)
	Xposf = Xpos/(Xpos+Xneg)
	Yposf = Ypos/(Ypos+Yneg)
	Xnegf = Xneg/(Xpos+Xneg)
	Ynegf = Yneg/(Ypos+Yneg)
	# calc monoallelic_raw if there is no pattern
	expected_overlap = Xposf*Yposf + Xnegf*Ynegf
	if expected_overlap==0: return 0
	# rescale where expected_overlap=>1, 1=>0, 0=>2
	if raw_overlap >= expected_overlap:
		return 1 - (raw_overlap - expected_overlap)/(1 - expected_overlap)
	else:
		return 2 - raw_overlap/expected_overlap

def monoallelic_raw(X, Y):
	return 1-numpy.mean([int(xi==yi) for xi,yi in zip(X,Y) if xi and yi])

def c57overlap(X, Y):
	return 1-numpy.mean([int(xi==yi) for xi,yi in zip(X,Y) if (xi == 1 or yi == 1) and xi and yi])

def castoverlap(X, Y):
	return 1-numpy.mean([int(xi==yi) for xi,yi in zip(X,Y) if (xi == -1 or yi == -1) and xi and yi])

def c57overlap_assym(X, Y):
	return 1-numpy.mean([int(xi==yi) for xi,yi in zip(X,Y) if yi == 1 and xi])

def castoverlap_assym(X, Y):
	return 1-numpy.mean([int(xi==yi) for xi,yi in zip(X,Y) if xi == -1 and yi])

def monoallelic_numgenes(X, Y):
	return 1-0.002*sum(int(xi==yi) for xi,yi in zip(X,Y) if xi and yi)

def monoallelic_numgenes_100(X, Y):
	return 1-0.01*sum(int(xi==yi) for xi,yi in zip(X,Y) if xi and yi)

def c57overlap_numgenes(X, Y):
	return 1-0.002*sum(int(xi==1 and yi==1) for xi,yi in zip(X,Y))

def castoverlap_numgenes(X, Y):
	return 1-0.002*sum(int(xi==-1 and yi==-1) for xi,yi in zip(X,Y))

def monoallelic_numgenes_norm(X, Y):
	return 1-0.002*sum(1 if xi==yi else -1 for xi,yi in zip(X,Y) if xi and yi)

def monoallelic_numgenes_norm_100(X, Y):
	return 1-0.01*sum(1 if xi==yi else -1 for xi,yi in zip(X,Y) if xi and yi)

def c57overlap_numgenes_norm(X, Y):
	return 1-0.002*sum(1 if xi==1 and yi==1 else -0.5 if xi==1 or yi==1 else 0 for xi,yi in zip(X,Y) if xi and yi)

def castoverlap_numgenes_norm(X, Y):
	return 1-0.002*sum(1 if xi==-1 and yi==-1 else -0.5 if xi==-1 or yi==-1 else 0 for xi,yi in zip(X,Y) if xi and yi)

if '__main__' == __name__:
	opts = argparse.ArgumentParser()
	opts.add_argument('rpkmf_alleles')
	opts.add_argument('--filter', nargs='+')
	opts.add_argument('-M', '--method', default='monoallelic', choices=['m', 'monoallelic', 'monoallelic_norm', 'monoallelic_norm2', 'c57overlap', 'castoverlap', 'c57overlap_assym', 'castoverlap_assym', 'numsamemono', 'numsameC57', 'numsameCAST', 'numsamemono_norm', 'numsameC57_norm', 'numsameCAST_norm', 'spearman', 'pearson', 'braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule', 'numsamemono100', 'numsamemono100_norm'])
	opts.add_argument('-L', '--linkage', default='complete', choices=['single', 'average', 'complete', 'linkage', 'weighted', 'centroid', 'median', 'ward'])
	opts.add_argument('-s', '--bootstrap', type=int)
	opts.add_argument('-S', '--states', default='3state', choices=['3state', 'fraction', 'diff', 'monoallelic'])
	opts.add_argument('--fig', default='which_allele_tree.pdf')
	opts.add_argument('-r', '--rpkmf_total')
	opts.add_argument('-t', '--threshold_rpkm', help='requires --rpkmf_total', type=float)
	opts.add_argument('-R', '--randomize', action='store_true')
	o = opts.parse_args()
	
	# load expression data
	expr_alleles = dr_tools.loadexpr([o.rpkmf_alleles], counts=True)
	samples_alleles = sorted([e for e in expr_alleles if e not in ('IDs', 'symbols')])
	if o.rpkmf_total is not None:
		expr_total = dr_tools.loadexpr([o.rpkmf_total], counts=False)
		exprt_samples = set(expr_total.samples)
	
	character_matrix = [] # 2D, values from state()
	samplenames = []
	
	for s1, s2 in zip(samples_alleles[::2], samples_alleles[1::2]):
		if o.filter is not None and not any(part in s1.rsplit('_',1)[0] for part in o.filter): continue
		samplename = s1.rsplit('_',1)[0]
		# check that sample labels are consistent
		if samplename != s2.rsplit('_',1)[0] and samplename in expr_total:
			continue
		
		if o.threshold_rpkm is None:
			character_matrix.append([state(e1,e2,o) for e1,e2 in zip(expr_alleles[s1], expr_alleles[s2])])
		else:
			if s1.split('_c57')[0] not in exprt_samples: continue
			character_matrix.append([state(e1,e2,o)  if t>=o.threshold_rpkm else state(0,0,o) for e1,e2,t in zip(expr_alleles[s1], expr_alleles[s2], expr_total[s1.split('_c57')[0]])])
		samplenames.append(samplename)
		
		if o.randomize:
			random.shuffle(character_matrix[-1])
	
	# add correlation methods
	if o.method == 'pearson':
		from scipy.stats import *
		m = lambda x,y : 1.0 - pearsonr(x,y)[0]
	elif  o.method == 'spearman':
		from scipy.stats import *
		m = lambda x,y : 1.0 - spearmanr(x,y)[0]
	elif o.method == 'm':
		m = lambda x,y : sum(0 if xi*yi==1 else 1 if xi*yi==0 else 5 for xi,yi in zip(x,y))
	elif o.method == 'monoallelic':
		m = monoallelic_raw
	elif o.method == 'monoallelic_norm':
		m = monoallelic_norm
	elif o.method == 'monoallelic_norm2':
		m = monoallelic_norm2
	elif o.method == 'c57overlap':
		m = c57overlap
	elif o.method == 'castoverlap':
		m = castoverlap
	elif o.method == 'c57overlap_assym':
		m = c57overlap_assym
	elif o.method == 'castoverlap_assym':
		m = castoverlap_assym
	elif o.method == 'numsamemono':
		m = monoallelic_numgenes
	elif o.method == 'numsamemono100':
		m = monoallelic_numgenes_100
	elif o.method == 'numsameC57':
		m = c57overlap_numgenes
	elif o.method == 'numsameCAST':
		m = castoverlap_numgenes
	elif o.method == 'numsamemono_norm':
		m = monoallelic_numgenes_norm
	elif o.method == 'numsamemono100_norm':
		m = monoallelic_numgenes_norm_100
	elif o.method == 'numsameC57_norm':
		m = c57overlap_numgenes_norm
	elif o.method == 'numsameCAST_norm':
		m = castoverlap_numgenes_norm
	else:
		m = o.method
	
	# make clusters
	exparray = character_matrix
	hcdists = hcluster.pdist(exparray, metric=m)
	hclinks = hcluster.linkage(hcdists, method=o.linkage)
	draw_order = hcluster.leaves_list(hclinks)
	
	# draw tree
	scipyhcluster.dendrogram(hclinks, labels=samplenames, leaf_rotation=90)
	pylab.subplots_adjust(bottom=0.3)
	pylab.ylabel('%s (linkage=%s)'%(o.method, o.linkage))
	if o.method in ('numsamemono', 'numsameC57', 'numsameCAST', 'numsamemono_norm', 'numsameC57_norm', 'numsameCAST_norm'):
		pylab.yticks([1.0, 0.8, 0.6, 0.4, 0.2, 0.0], [0, 100, 200, 300, 400, 500])
	elif o.method in ('numsamemono100','numsamemono100_norm'):
		pylab.yticks([1.0, 0.8, 0.6, 0.4, 0.2, 0.0], [0, 20, 40, 60, 80, 100])
	pylab.savefig(o.fig)
	
	# bootstrap
	if o.bootstrap is not None:
		linkingpatterns = bootstraphcluster(exparray, o.bootstrap, None, m, hclinks)
		for pattern in linkingpatterns:
			print pattern.getstring(samplenames)

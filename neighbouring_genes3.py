from __future__ import division
import argparse, dr_tools, numpy, random, math, pylab

'''
altered from neighbouring_genes.py
calculates fraction of same parental origin for monoallelic neighbouring genes
'''

class Gene:
	def __init__(self, ID, TSS, strand):
		self.ID = ID
		self.TSS = TSS
		self.neighbours = []
		self.states = []
		self.strand = strand
		
	def set_state(self, state):
		self.states.append(state)
	
	def consistencyval(self, use_neighbours):
		if not self.neighbours: return []
		
		if use_neighbours:
			neighbours = self.neighbours
		else:
			neighbours = random.sample(self.non_chr_genes, len(self.neighbours))
		
		ret = []
		for neighbour in neighbours:
			ret.extend([state_overlap(state1, state2) for state1,state2 in zip(self.states, neighbour.states)])
		ret = [r for r in ret if r is not None]
		if not ret: return []
		return [numpy.mean(ret)]

def state_overlap(state1, state2):
	return int(state1 == state2) if state1 and state2 else None

def calc_state(maternal, paternal, o):
	if (maternal>0, paternal>0) == (True,False): return 1
	elif (maternal>0, paternal>0) == (False,True): return -1
	else:
		return 0

if '__main__' == __name__:
	opts = argparse.ArgumentParser()
	opts.add_argument('rpkmf_alleles')
	opts.add_argument('genePred', default='/mnt/crick/danielr/twocellstage/mouse/annotation/mm9_refGene_31Jul2011_norandom.txt', nargs='?')
	opts.add_argument('--numadjacent', type=int, default=1000)
	opts.add_argument('--filter', nargs='+')
	opts.add_argument('--exclude', nargs='+', default=[])
	opts.add_argument('--bindist', type=float, default=[0,1e20], nargs='+')
	opts.add_argument('--onlydiffstrand', action='store_true')
	opts.add_argument('-o', '--figure', default='neighbouring_genes3.pdf')
	args = opts.parse_args()
	
	# load expression data
	expr_alleles = dr_tools.loadexpr([args.rpkmf_alleles], counts=True)
	samples_alleles = sorted(e for e in expr_alleles if e not in ('IDs', 'symbols') and (args.filter is None or any(part in e for part in args.filter)) and not any(part in e for part in args.exclude))
	allowed_IDs = set(IDs.split('+')[0] for IDs in expr_alleles['IDs'])
	
	fractions_to_show = list()
	vals_real = list()
	vals_ctrl = list()
	labels = list()
	bootstrap_output = list()
	
	
	# sort the genes by posiotion
	# only include transcripts which are the first ID in the entry of the rpkm file
	for mindist,maxdist in zip(args.bindist[:-1], args.bindist[1:]):
		genes_per_chr = dict()
		ID_to_gene = dict()
		for p in dr_tools.splitlines(args.genePred):
			ID = p[1]
			if ID in allowed_IDs:
				chromosome = p[2]
				if chromosome in ('chrX', 'chrY'): continue
				if not chromosome in genes_per_chr: genes_per_chr[chromosome] = []
				genes_per_chr[chromosome].append(Gene(ID, int(p[4]) if p[3]=='+' else int(p[5]), p[3]))
		for chromosome in genes_per_chr:
			genes_per_chr[chromosome].sort(key=lambda gene: gene.TSS)
			non_chr_genes = [gene for chr_key in genes_per_chr if chr_key != chromosome for gene in genes_per_chr[chromosome]]
			for gene_i, gene in enumerate(genes_per_chr[chromosome]):
				ID_to_gene[gene.ID] = gene
				# only include neighbours in the 'forward' direction, to avoid dependence in stat tests later on
				gene.neighbours = genes_per_chr[chromosome][gene_i+1:args.numadjacent+gene_i+1]
				gene.neighbours = [other for other in gene.neighbours if mindist <= abs(other.TSS-gene.TSS) < maxdist]
				if args.onlydiffstrand:
					gene.neighbours = [other for other in gene.neighbours if other.strand != gene.strand]
			
				gene.non_chr_genes = non_chr_genes
	
		for s1, s2 in zip(samples_alleles[::2], samples_alleles[1::2]):
			samplename = s1.rsplit('_',1)[0]
		
			# check that sample labels are consistent
			if samplename != s2.rsplit('_',1)[0] and samplename in expr_total:
				continue
		
			# set values from gene expression
			for ID_field, e1, e2 in zip(expr_alleles['IDs'], expr_alleles[s1], expr_alleles[s2]):
				key = ID_field.split('+')[0]
				if not key in ID_to_gene: continue
				ID_to_gene[key].set_state(calc_state(e1,e2,args))
		
		# give list of 1 (same parental origin) and 0 (different) for each pair of monoallelic genes
		vals_real.append([agreement for chromosome in genes_per_chr for gene in genes_per_chr[chromosome] for agreement in gene.consistencyval(True)])
		vals_ctrl.extend([agreement for chromosome in genes_per_chr for gene in genes_per_chr[chromosome] for agreement in gene.consistencyval(False)])
		fractions_to_show.append(numpy.mean(vals_real[-1]))
		bootstrap_output.append(dr_tools.bootstrap(numpy.mean, [vals_real[-1]], nullval=0.5, processes=10))
		labels.append('%.1e-%.1e'%(mindist,maxdist))
	fractions_to_show.append(numpy.mean(vals_ctrl))
	labels.append('neg ctrl')
	bootstrap_output.append(dr_tools.bootstrap(numpy.mean, [vals_ctrl], nullval=0.5, processes=10))
	
	# plot
	Xpos = range(len(fractions_to_show))
	Xpos_mid = [x+0.4 for x in Xpos]
	pylab.bar(Xpos, fractions_to_show, facecolor='#cccccc', linewidth=0)
	pylab.errorbar(Xpos_mid, fractions_to_show, [[F-O[1] for F,O in zip(fractions_to_show,bootstrap_output)], [O[2]-F for F,O in zip(fractions_to_show,bootstrap_output)]], fmt=None, ecolor='k')
	
	for label_i, v_real in enumerate(vals_real):
		pval = dr_tools.permutationtest(numpy.mean, v_real, vals_ctrl, controls=10000)
		labels[label_i] += '\n%.1e'%pval
	
	pylab.ylabel('fraction same parental origin\n of pairs of monoallelic genes')
	pylab.subplots_adjust(bottom=0.3)
	pylab.xticks(Xpos_mid, labels, rotation=90)
	pylab.savefig(args.figure)

from __future__ import division
import argparse, dr_tools, numpy, random, pylab


class Gene:
	def __init__(self, ID, TSS):
		self.ID = ID
		self.TSS = TSS # genomic position on the chromosome
		self.maternalreads = []
		self.paternalreads = []
		
	def add_state(self, maternal, paternal):
		self.maternalreads.append(maternal)
		self.paternalreads.append(paternal)
		
def chrvalue(name):
	try:
		return int(name[3:])
	except:
		return {'X':101, 'Y':102, 'M':103}[name[3:]]

if '__main__' == __name__:
	opts = argparse.ArgumentParser()
	opts.add_argument('rpkmf_alleles')
	opts.add_argument('--genePred', default='/mnt/crick/danielr/twocellstage/mouse/annotation/mm9_refGene_31Jul2011_norandom.txt')
	opts.add_argument('--filter', nargs='+')
	opts.add_argument('-o', '--figurefile', default='monoallelic_by_chr.pdf')
	args = opts.parse_args()
	
	# load expression data
	expr_alleles = dr_tools.loadexpr([args.rpkmf_alleles], counts=True)
	samples_alleles = sorted(e for e in expr_alleles if e not in ('IDs', 'symbols') and (args.filter is None or any(part in e for part in args.filter)))
	
	# sort the genes by position
	# only include transcripts which are the first ID in the entry of the rpkm file
	allowed_IDs = set(IDs.split('+')[0] for IDs in expr_alleles['IDs'])
	genes_per_chr = dict()
	ID_to_gene = dict()
	for p in dr_tools.splitlines(args.genePred):
		ID = p[1]
		if ID in allowed_IDs:
			chromosome = p[2]
			if 'random' in chromosome: continue
			if not chromosome in genes_per_chr: genes_per_chr[chromosome] = []
			genes_per_chr[chromosome].append(Gene(ID, int(p[4]) if p[3]=='+' else int(p[5])))
	for chromosome in genes_per_chr:
		genes_per_chr[chromosome].sort(key=lambda gene: gene.TSS)
		for gene_i, gene in enumerate(genes_per_chr[chromosome]):
			ID_to_gene[gene.ID] = gene
	
	samples = []
	
	for s1, s2 in zip(samples_alleles[::2], samples_alleles[1::2]):
		samplename = s1.rsplit('_',1)[0]
		
		# check that sample labels are consistent
		if samplename != s2.rsplit('_',1)[0]: raise Exception
		if not 'c57only' in s1 or not 'castonly' in s2: raise Exception
		
		samples.append(samplename)
		
		# use expression values
		for ID, pat, mat in zip(expr_alleles['IDs'], expr_alleles[s1], expr_alleles[s2]):
			try:
				ID_to_gene[ID.split('+')[0]].add_state(mat, pat)
			except KeyError: pass
	
	# chromosome order, from chr1 to chrM
	chr_order = sorted(genes_per_chr.keys(), key=lambda n: chrvalue(n))
	labels = []
	
	maternal_fraction = []
	for chromosome in chr_order:
		m_reads = [sum(g.maternalreads[i] for g in genes_per_chr[chromosome] if g.paternalreads) for i,s in enumerate(samples)]
		p_reads = [sum(g.paternalreads[i] for g in genes_per_chr[chromosome] if g.paternalreads) for i,s in enumerate(samples)]
		if any(M+P==0 for M,P in zip(m_reads, p_reads)):
			continue
		maternal_fraction.append(numpy.mean([M/(M+P) for M,P in zip(m_reads, p_reads)]))
		labels.append(chromosome)
	
	# draw
	if args.filter: pylab.title(' '.join(args.filter))
	pylab.bar(range(len(labels)), maternal_fraction, facecolor='w')
	pylab.xticks([x+0.4 for x in range(len(labels))], labels, rotation=90, horizontalalignment='center')
	pylab.ylabel('maternal fraction of reads')
	pylab.ylim(0,1)
	pylab.subplots_adjust(bottom=0.3)
	pylab.savefig(args.figurefile)

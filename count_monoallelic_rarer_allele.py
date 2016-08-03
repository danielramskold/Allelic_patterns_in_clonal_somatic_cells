from __future__ import division
import argparse, dr_tools, numpy, random, pylab, math
from collections import defaultdict

def rarer_allele_only_fraction(s, expra, minor_alleles, minreads):
	only_rarer = 1
	bi_allelic = 2
	gene_status = defaultdict(int)
	for ai, allele_order in minor_alleles.items():
		minor_a, major_a = allele_order
		gene = expra['symbols'][ai]
		if expra[s+minor_a][ai]>=minreads:
			if expra[s+major_a][ai]>=minreads:
				gene_status[gene] = max(gene_status[gene], bi_allelic)
			elif expra[s+major_a][ai]==0:
				gene_status[gene] = max(gene_status[gene], only_rarer)
		
	mono_rarer = gene_status.values().count(only_rarer)
	bi = gene_status.values().count(bi_allelic)
	return mono_rarer, bi


def estimated_monoallelic_fraction(s, expra, minor_alleles, minreads):
	mono_rarer, bi = rarer_allele_only_fraction(s, expra, minor_alleles, minreads)
	if not mono_rarer and not bi: return float('NaN')
	return mono_rarer*2/(2*mono_rarer+bi)

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allelehits', required=True)
	parser.add_argument('-m', '--minreads', type=int, default=2)
	parser.add_argument('-gi', '--genelist_include')
	parser.add_argument('-ge', '--genelist_exclude', nargs='+')
	parser.add_argument('--figure')
	o = parser.parse_args()
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	
	a1 = '_c57only'
	a2 = '_castonly'
	all_samples = [s[:-len(a1)] for s in expra.samples[::2]]
	
	genelist_include = None if o.genelist_include is None else set(dr_tools.loadlist(o.genelist_include))
	genelist_exclude = None if o.genelist_exclude is None else set.union(*(set(dr_tools.loadlist(filename)) for filename in o.genelist_exclude))
	minor_alleles = dict()
	for ai, sym in enumerate(expra['symbols']):
		if genelist_include is not None and sym not in genelist_include: continue
		if genelist_exclude is not None and sym in genelist_exclude: continue
		countsum1 = sum(expra[s+a1][ai]>=o.minreads for s in all_samples)
		countsum2 = sum(expra[s+a2][ai]>=o.minreads for s in all_samples)
		minor_alleles[ai] = (a2,a1) if countsum1 >= countsum2 else (a1,a2)
	
	fractions = []
	for sample in all_samples:
		f = estimated_monoallelic_fraction(sample, expra, minor_alleles, o.minreads)
		print sample, f
		if not math.isnan(f): fractions.append(f)
	print 'average', numpy.mean(fractions)
	
	if o.figure:
		pylab.boxplot(fractions, whis=1000)
		pylab.ylim(0,1)
		pylab.savefig(o.figure)

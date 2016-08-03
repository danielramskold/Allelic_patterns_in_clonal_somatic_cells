from __future__ import division
import argparse, dr_tools, numpy, random
from scipy import stats
from collections import defaultdict

def give_sharing(clone_samples, minor_alleles, expra):
	num_shared = defaultdict(float)
	minreads = 1
	for ai, allele_order in minor_alleles.items():
		minor_a, major_a = allele_order
		cells_with_minor_allele = sum(expra[s+minor_a][ai]>=minreads and not expra[s+major_a][ai]>=minreads for s in clone_samples)
		expressing_cells = sum(expra[s+minor_a][ai]>=minreads or expra[s+major_a][ai]>=minreads for s in clone_samples)
		if cells_with_minor_allele >= max(2, expressing_cells*0.9):
			gene = expra['symbols'][ai]
			num_shared[gene] = 1
	return sum(num_shared.values())

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allelehits', required=True)
	parser.add_argument('-s', '--samplelist', required=True)
	parser.add_argument('-r', '--randomisations', type=int, default=100)
	parser.add_argument('-gi', '--genelist_include')
	parser.add_argument('-ge', '--genelist_exclude')
	o = parser.parse_args()
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	
	a1 = '_c57only'
	a2 = '_castonly'
	
	all_samples = [s[:-len(a1)] for s in expra.samples[::2]]
	clone_samples = dr_tools.loadlist(o.samplelist)
	genelist_include = None if o.genelist_include is None else set(dr_tools.loadlist(o.genelist_include))
	genelist_exclude = None if o.genelist_exclude is None else set(dr_tools.loadlist(o.genelist_exclude))
	
	minor_alleles = dict()
	for ai, sym in enumerate(expra['symbols']):
		if genelist_include is not None and sym not in genelist_include: continue
		if genelist_exclude is not None and sym in genelist_exclude: continue
		countsum1 = sum(expra[s+a1][ai]>0 for s in all_samples)
		countsum2 = sum(expra[s+a2][ai]>0 for s in all_samples)
		minor_alleles[ai] = (a2,a1) if countsum1 >= countsum2 else (a1,a2)
	
	clone_shared = give_sharing(clone_samples, minor_alleles, expra)
	random_shared = [give_sharing(random.sample(all_samples, len(clone_samples)), minor_alleles, expra) for i in range(o.randomisations)]
	p_crude =  max(1, sum(r>=clone_shared for r in random_shared))/len(random_shared)
	p_normal = stats.ttest_1samp(random_shared, clone_shared)[1]
	excess_genes = clone_shared - numpy.mean(random_shared)
	print p_crude, p_normal, excess_genes, random_shared[:10], clone_shared
	

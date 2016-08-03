from __future__ import division
import argparse, dr_tools, numpy, random
from scipy import stats
from collections import defaultdict

def give_sharing(clone_samples, minor_alleles, expra):
	num_shared = defaultdict(float)
	for ai, allele_order in minor_alleles.items():
		minor_a, major_a = allele_order
		cells_with_minor_allele = sum(expra[s+minor_a][ai]>0 and not expra[s+major_a][ai]>0 for s in clone_samples)
		if cells_with_minor_allele >= 2:
			expressing_cells = sum(expra[s+minor_a][ai]>0 or expra[s+major_a][ai]>0 for s in clone_samples)
			gene = expra['symbols'][ai]
			num_shared[gene] = max(cells_with_minor_allele/expressing_cells, num_shared[gene])
	return sum(num_shared.values())

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allelehits', required=True)
	parser.add_argument('-s', '--samplelist', required=True)
	parser.add_argument('-r', '--randomisations', type=int, default=100)
	parser.add_argument('-gi', '--genelist_include')
	parser.add_argument('-ge', '--genelist_exclude')
	parser.add_argument('--all_clonal')
	o = parser.parse_args()
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	
	a1 = '_c57only'
	a2 = '_castonly'
	
	all_samples = [s[:-len(a1)] for s in expra.samples[::2]]
	clone_samples = dr_tools.loadlist(o.samplelist)
	samples_to_randomly_pick_from = all_samples if o.all_clonal is None else set(all_samples)-set(dr_tools.loadlist(o.all_clonal))
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
	random_shared = [give_sharing(random.sample(samples_to_randomly_pick_from, len(clone_samples)), minor_alleles, expra) for i in range(o.randomisations)]
	p_crude =  max(1, sum(r>=clone_shared for r in random_shared))/len(random_shared)
	p_normal = stats.ttest_1samp(random_shared, clone_shared)[1]
	excess_genes = clone_shared - numpy.mean(random_shared)
	print p_crude, p_normal, excess_genes, random_shared[:10], clone_shared
	print clone_shared - dr_tools.pointinarray(random_shared, 0.025), clone_shared - dr_tools.pointinarray(random_shared, 0.975)

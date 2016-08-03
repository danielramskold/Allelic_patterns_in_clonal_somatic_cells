from __future__ import division
import argparse, dr_tools, numpy, random, pylab
from collections import defaultdict

def give_sharing_v1(clone_samples, minor_alleles, expra, minreads):
	num_shared = defaultdict(float)
	for ai, allele_order in minor_alleles.items():
		minor_a, major_a = allele_order
		cells_with_minor_allele = sum(expra[s+minor_a][ai]>=minreads and not expra[s+major_a][ai]>0 for s in clone_samples)
		if cells_with_minor_allele >= 2:
			expressing_cells = sum(expra[s+minor_a][ai]>=minreads or expra[s+major_a][ai]>=minreads for s in clone_samples)
			gene = expra['symbols'][ai]
			num_shared[gene] = max(cells_with_minor_allele/expressing_cells, num_shared[gene])
	return sum(num_shared.values())

def give_sharing_v3(clone_samples, minor_alleles, expra, minreads):
	num_shared = defaultdict(float)
	for ai, allele_order in minor_alleles.items():
		minor_a, major_a = allele_order
		cells_with_minor_allele = sum(expra[s+minor_a][ai]>=minreads and not expra[s+major_a][ai]>=minreads for s in clone_samples)
		expressing_cells = sum(expra[s+minor_a][ai]>=minreads or expra[s+major_a][ai]>=minreads for s in clone_samples)
		if cells_with_minor_allele >= max(2, expressing_cells*0.9):
			gene = expra['symbols'][ai]
			num_shared[gene] = 1
	return sum(num_shared.values())

def give_sharing_v0(clone_samples, minor_alleles, expra, minreads):
	num_shared = defaultdict(float)
	for ai, allele_order in minor_alleles.items():
		minor_a, major_a = allele_order
		cells_with_minor_allele = sum(expra[s+minor_a][ai]>=minreads and not expra[s+major_a][ai]>=minreads for s in clone_samples)
		if cells_with_minor_allele >=2:
			gene = expra['symbols'][ai]
			num_shared[gene] = 1
	return sum(num_shared.values())

def give_sharing_90pwMajor(clone_samples, minor_alleles, expra, minreads):
	num_shared = defaultdict(float)
	for ai, allele_order in minor_alleles.items():
		minor_a, major_a = allele_order
		cells_with_minor_allele = sum(expra[s+minor_a][ai]>=minreads and not expra[s+major_a][ai]>=minreads for s in clone_samples)
		cells_with_major_allele = sum(expra[s+major_a][ai]>=minreads and not expra[s+minor_a][ai]>=minreads for s in clone_samples)
		expressing_cells = sum(expra[s+minor_a][ai]>=minreads or expra[s+major_a][ai]>=minreads for s in clone_samples)
		if max(cells_with_minor_allele, cells_with_minor_allele) >= max(2, expressing_cells*0.9):
			gene = expra['symbols'][ai]
			num_shared[gene] = 1
	return sum(num_shared.values())

def give_sharing_min2wMajor(clone_samples, minor_alleles, expra, minreads):
	num_shared = defaultdict(float)
	for ai, allele_order in minor_alleles.items():
		minor_a, major_a = allele_order
		cells_with_minor_allele = sum(expra[s+minor_a][ai]>=minreads and not expra[s+major_a][ai]>=minreads for s in clone_samples)
		cells_with_major_allele = sum(expra[s+major_a][ai]>=minreads and not expra[s+minor_a][ai]>=minreads for s in clone_samples)
		expressing_cells = sum(expra[s+minor_a][ai]>=minreads or expra[s+major_a][ai]>=minreads for s in clone_samples)
		if max(cells_with_minor_allele, cells_with_minor_allele) >= 2:
			gene = expra['symbols'][ai]
			num_shared[gene] = 1
	return sum(num_shared.values())

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allelehits', required=True)
	parser.add_argument('-s', '--samplelist', required=True, nargs='+')
	parser.add_argument('-r', '--randomisations', type=int, default=100)
	parser.add_argument('-gi', '--genelist_include')
	parser.add_argument('-ge', '--genelist_exclude', nargs='+')
	parser.add_argument('-M', '--metric', default='v1', choices=['v0', 'v1', 'v3', '90pwM', '2wM'])
	parser.add_argument('-m', '--minreads', type=int, default=1)
	parser.add_argument('-o', '--figure', required=True)
	parser.add_argument('--all_clonal')
	parser.add_argument('-d', '--divide_by_total', action='store_true')
	o = parser.parse_args()
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	
	a1 = '_c57only'
	a2 = '_castonly'
	
	all_samples = [s[:-len(a1)] for s in expra.samples[::2]]
	samples_to_randomly_pick_from = all_samples if o.all_clonal is None else set(all_samples)-set(dr_tools.loadlist(o.all_clonal))
	genelist_include = None if o.genelist_include is None else set(dr_tools.loadlist(o.genelist_include))
	genelist_exclude = None if o.genelist_exclude is None else set.union(*(set(dr_tools.loadlist(filename)) for filename in o.genelist_exclude))
	
	minor_alleles = dict()
	included_genes = set()
	for ai, sym in enumerate(expra['symbols']):
		if genelist_include is not None and sym not in genelist_include: continue
		if genelist_exclude is not None and sym in genelist_exclude: continue
		countsum1 = sum(expra[s+a1][ai]>0 for s in all_samples)
		countsum2 = sum(expra[s+a2][ai]>0 for s in all_samples)
		if countsum1 and countsum2:
			minor_alleles[ai] = (a2,a1) if countsum1 >= countsum2 else (a1,a2)
			included_genes.add(sym)
	
	give_sharing = {'v0':give_sharing_v0, 'v1': give_sharing_v1, 'v3':give_sharing_v3, '90pwM':give_sharing_90pwMajor, '2wM':give_sharing_min2wMajor}[o.metric]
	
	total_genes = len(included_genes) if o.divide_by_total else 1
	
	mid = []
	error_below = []
	error_above = []
	for samplelist in o.samplelist:
		clone_samples = dr_tools.loadlist(samplelist)
		clone_shared = give_sharing(clone_samples, minor_alleles, expra, o.minreads)
		random_shared = [give_sharing(random.sample(samples_to_randomly_pick_from, len(clone_samples)), minor_alleles, expra, o.minreads) for i in range(o.randomisations)]
		
		low, median, high = dr_tools.pointinarray(random_shared, [0.025, 0.5, 0.975])
		
		mid.append((clone_shared - median)/total_genes)
		error_above.append((median - low)/total_genes)
		error_below.append((high - median)/total_genes)
		print clone_shared-low, clone_shared-median, clone_shared-high, error_above[-1], error_below[-1]
	
	pylab.bar(range(1, len(o.samplelist)+1), mid, yerr=[error_below, error_above], linewidth=0, ecolor='k', facecolor='#cccccc')
	pylab.savefig(o.figure)

from __future__ import print_function, division
import argparse, dr_tools, numpy, random
from collections import defaultdict
from scipy import stats

def parse_table(filename):
	column = ''
	header = []
	markers = defaultdict(list)
	marker_order = []
	convert_element = {'+':1, 'int':0, '':-1, '\x00':-1}
	with open(filename, 'rU') as infh:
		for line in infh:
			ele = line.rstrip()
			if ele == 'Subpop. #':
				column = 'header'
			elif column == 'header' and ele.isdigit():
				header.append(int(ele))
			elif (column == 'header' and not ele.isdigit()) or len(markers[column]) == len(header):
				column = ele
				marker_order.append(ele)
			else:
				markers[column].append(convert_element[ele])
	return header, markers, marker_order
			

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-r', '--rpkmfile')
	parser.add_argument('--tableS4', default='tableS4.txt')
	parser.add_argument('--to_cytof_markers', default='symbol_to_cytof_marker.txt')
	parser.add_argument('--shuffle_patterns', action='store_true')
	parser.add_argument('-o', '--sample_list_prefix')
	o = parser.parse_args()

	header, markers, marker_order = parse_table(o.tableS4)
	
	gene_to_marker = dict(dr_tools.splitlines(o.to_cytof_markers))
	marker_order = [m for m in marker_order if m in gene_to_marker.values()]
	if not o.shuffle_patterns:
		pop_cytof_pattern = dict((pop, [markers[m][popi] for m in marker_order]) for popi,pop in enumerate(header))
	else:
		pop_cytof_pattern = dict((pop, random.shuffle([markers[m][popi] for m in marker_order])) for popi,pop in enumerate(header))
	exprt = dr_tools.loadexpr(o.rpkmfile)
	random.seed(0)

	midexpr_symi_all_D = dict()
	for symi, sym in enumerate(exprt['symbols']):
		if sym not in gene_to_marker: raise Exception(dr_tools.join(sym, 'sym'))
		if gene_to_marker[sym] not in markers: raise Exception(dr_tools.join(gene_to_marker[sym], 'cytof'))
		midexpr_symi_all_D[gene_to_marker[sym]] = (numpy.mean([exprt[s][symi] for s in exprt.samples]), symi)
	midexpr_symi_all = [midexpr_symi_all_D[m] for m in marker_order]
	sym_order = [midexpr_symi_all_D[m][1] for m in marker_order]

	pop_counts = dict((pop,0) for pop in pop_cytof_pattern)
	pop_samples =defaultdict(list)

	for sample in exprt.samples:
		relexpr = [exprt[sample][symi]/midall for midall, symi in midexpr_symi_all]

		r_max = -2
		pop_max_L = []
		for pop, pattern in pop_cytof_pattern.items():
			try: r = stats.spearmanr(relexpr, pattern)[0]
			except ZeroDivisionError: r = 0
			if r>r_max:
				pop_max_L = [pop]
				r_max = r
			elif r==r_max:
				pop_max_L.append(pop)
		pop_max = random.choice(pop_max_L)
		
		pop_counts[pop_max] += 1
		pop_samples[pop_max].append(sample)
	print(pop_counts)

	if o.sample_list_prefix is not None:
		for pop, samples in pop_samples.items():
			dr_tools.printlist('%s%d.txt'%(o.sample_list_prefix, pop), samples)
		


		
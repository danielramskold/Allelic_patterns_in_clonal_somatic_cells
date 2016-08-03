import argparse
from collections import defaultdict

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('plot4bars_vocal_output')
	o = parser.parse_args()
	
	clonemono = defaultdict(set)
	with open(o.plot4bars_vocal_output, 'r') as infh:
		for line in infh:
			if ' shared genes, ' in line:
				clone = line.split(' shared genes, ')[0]
				genes = line.rstrip().split(':')[-1].split(', ')
				clonemono[clone].update(set(genes))
	
	num_shared = defaultdict(set)
	any_gene = set.union(*clonemono.values())
	for gene in any_gene:
		print gene

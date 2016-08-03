from __future__ import division
import argparse, dr_tools

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--genePred', required=True)
	parser.add_argument('-o', '--genelist_out', help="gives genes on -i chromosomes but not on -e", required=True)
	parser.add_argument('-i', '--chrom_include', nargs='+', help="'rest' matches all not in -i or -e, 'random' matches e.g. chr1_random", required=True)
	parser.add_argument('-e', '--chrom_exclude', nargs='+', help="'rest' matches all not in -i or -e, 'random' matches e.g. chr1_random", default=[])
	o = parser.parse_args()
	
	all_chromosomes = set(dr_tools.loadlist(o.genePred, 2))
	
	include = set(c for c in all_chromosomes if c in o.chrom_include)
	exclude = set(c for c in all_chromosomes if c in o.chrom_exclude)
	if 'random' in o.chrom_include:
		include.update(set(c for c in all_chromosomes if '_random' in c and c not in exclude))
	if 'random' in o.chrom_exclude:
		exclude.update(set(c for c in all_chromosomes if '_random' in c and c not in include))
	if 'rest' in o.chrom_include:
		include.update(all_chromosomes-exclude)
	if 'rest' in o.chrom_exclude:
		exclude.update(all_chromosomes-include)
	
	
	genes_incl = set()
	genes_excl = set()
	
	for p in dr_tools.splitlines(o.genePred):
		symbol = p[12]
		chromosome = p[2]
		if chromosome in include: genes_incl.add(symbol)
		elif chromosome in exclude: genes_excl.add(symbol)
	
	with open(o.genelist_out, 'w') as outfh:
		for gene in genes_incl-genes_excl:
			print >>outfh, gene

import dr_tools, argparse

def load_geneset(ID_to_symbol, filename):
	allowed_symbols = set(ID_to_symbol.values())
	allowed_IDs = set(ID_to_symbol.keys())
	geneset_symbols = set()
	entries = dict()
	for genes in dr_tools.splitlines(filename, sep=';'):
		genes_sym = set(gene for gene in genes if gene in allowed_symbols)
		genes_sym |= set(ID_to_symbol[gene] for gene in genes if gene in allowed_IDs and gene not in genes_sym)
		for sym in genes_sym:
			entries[sym] = frozenset(genes_sym)
	return entries

def overlap_of_2(entries1, entries2):
	set1 = set(entries1)
	set2 = set(entries2)
	set1_unique_c = len(set(entries1[sym] for sym in (set1-set2)))
	set2_unique_c = len(set(entries2[sym] for sym in (set2-set1)))
	common_c = len(set(entries1[sym] for sym in (set2&set1)))
	common_c2 = len(set(entries2[sym] for sym in (set2&set1)))
	if not common_c == common_c2: raise Exception
	saygenes = []
	for genes in set(entries2[sym] for sym in (set2&set1)):
		saygenes.append(';'.join(list(genes)))
	return set1_unique_c, common_c, set2_unique_c, ', '.join(saygenes)

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-A', '--annotationfile', default='/mnt/crick/danielr/Xandclones_BR/BR_fibroblasts/snp-call/more_formats/mm9_ensembl_refseq_norandom_11Apr2012_genesymbols.txt')
	parser.add_argument('-a', '--set1', required=True)
	parser.add_argument('-b', '--set2', required=True)
	parser.add_argument('-ge', '--disallowedgenes', nargs='+')
	o = parser.parse_args()
	
	if o.disallowedgenes:
		disallowedgenes = set()
		for filename in o.disallowedgenes:
			disallowedgenes.update(set(dr_tools.loadlist(filename)))
	else:
		disallowedgenes = None
	
	ID_to_symbol = dict((p[1], p[12]) for p in dr_tools.splitlines(o.annotationfile) if disallowedgenes is None or p[12] not in disallowedgenes)
	
	print dr_tools.join(overlap_of_2(load_geneset(ID_to_symbol, o.set1), load_geneset(ID_to_symbol, o.set2)))
	print len(set(ID_to_symbol.values()))

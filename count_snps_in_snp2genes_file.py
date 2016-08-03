import argparse, dr_tools

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('snp2genes', help='/mnt/kauffman/danielr/crick/Xandclones_BR/snp-validation/snp2genes/ensembl__nooverlap_ra5val_3percent.txt')
	parser.add_argument('-g', '--genelist')
	o = parser.parse_args()
	
	
	if o.genelist: allowedgenes = set(dr_tools.loadlist(o.genelist))
	
	num_snps = 0
	num_genes = 0
	for p in dr_tools.splitlines(o.snp2genes):
		if o.genelist and p[0] not in allowedgenes: continue
		num_snps += int(p[1])
		num_genes += 1
	print 'snps:', num_snps
	print 'genes:', num_genes

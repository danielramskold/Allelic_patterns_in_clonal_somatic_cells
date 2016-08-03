import argparse, dr_tools
from collections import defaultdict

def fromannotationline(line, l_annotationtype=0):
	inferred_strand = False
	if l_annotationtype == 0:
		# from refGene.txt
		p = line.rstrip('\r\n').split("\t")
		exonstarts = [int(f) for f in p[9].split(",")[:-1]]	# start positions for exons for the gene
		exonends = [int(f) for f in p[10].split(",")[:-1]]
		ID = p[1]
		chromosome = p[2]
		genename = p[12]
		strand = p[3]
		cdsstart = min(int(p[7]), int(p[6]))
		cdsend = max(int(p[7]), int(p[6]))
	return (chromosome, strand, cdsstart, cdsend, exonstarts, exonends, genename, ID, inferred_strand)

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('genepred')
	parser.add_argument('snplist', help='e.g. validated_cast_c57_snps.txt')
	parser.add_argument('outfile', help='same output format as validated_mm9_refseq_snp2genes.txt')
	parser.add_argument('--include_overlap', action='store_true')
	o = parser.parse_args()
	
	with open(o.genepred) as infh:
		for line in infh:
			chromosome, strand, cdsstart, cdsend, exonstarts, exonends, genename, ID, inferred_strand = fromannotationline(line)
			for start, end in zip(exonstarts, exonends):
				exon = dr_tools.Cregion(chromosome, start, end)
				exon.gene = genename
				exon.addtowindows()
	
	snps_per_gene = defaultdict(list)
	
	for p in dr_tools.splitlines(o.snplist):
		# e.g. chr11   117883408       C       A       0       1.00    -1.00   0.90    0.10    0.71    0.29
		
		
		if ',' in p[2] or ',' in p[3]: continue # added 18 Dec, since snp_stats2.py -S removes these SNPs anyway
		
		chromosome = p[0]
		position = int(p[1])-1
		genes = set(exon.gene for exon in dr_tools.Cregion.overlappingpoint(chromosome, position))
		if o.include_overlap:
			for gene in genes:
				snps_per_gene[gene].append('%s:%s'%(p[0], p[1]))
		else:
			if len(genes) == 1: # don't allow overlapping genes, exclude those SNPs
				gene = list(genes)[0]
				try: snps_per_gene[gene].append('%s:%s'%(p[0], p[1]))
				except:
					print p
					raise
	
	with open(o.outfile, 'w') as outfh:
		for gene, snps in snps_per_gene.items():
			print >>outfh, dr_tools.join(gene, len(snps), ';'.join(sorted(snps)))

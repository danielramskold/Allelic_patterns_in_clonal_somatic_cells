from __future__ import division
import argparse, dr_tools, numpy
from collections import defaultdict

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('sample_and_chromosome_list')
	parser.add_argument('output_file')
	parser.add_argument('-A', '--annotationfile')
	parser.add_argument('-a', '--allelehits')
	o = parser.parse_args()
	
	exprr = dr_tools.loadexpr(o.allelehits, False)
	expra = dr_tools.loadexpr(o.allelehits, True)
	
	chrom_to_IDs = defaultdict(set)
	for p in dr_tools.splitlines(o.annotationfile):
		chrom = p[2]
		sym = p[12]
		ID = p[1]
		chrom_to_IDs[chrom].add(ID)
	
	samples_set = set(expra.samples)
	
	with open(o.sample_and_chromosome_list) as infh:
		for line in infh:
			p = line.split()
			chrom = p[1]
			s_c57 = p[0]+'_c57only'
			s_cast = p[0]+'_castonly'
			if p[0] not in samples_set: continue
			for ai, ID in enumerate(expra['IDs']):
				if ID.split('+')[0] in chrom_to_IDs[chrom]:
					expra[s_c57][ai] = 0.0
					expra[s_cast][ai] = 0.0
	
	dr_tools.writeexpr(o.output_file, exprr, expra)

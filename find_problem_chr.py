from __future__ import division
import argparse, dr_tools, numpy
from collections import defaultdict

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-A', '--annotationfile')
	parser.add_argument('-a', '--allelehits')
	parser.add_argument('-t', '--threshold', type=float, default=0.5)
	parser.add_argument('-m', '--addition', nargs=2, action='append', default=[])
	parser.add_argument('-l', '--removal', nargs=2, action='append', default=[])
	o = parser.parse_args()
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	
	chrom_to_ai = defaultdict(list)
	for p in dr_tools.splitlines(o.annotationfile):
		chrom = p[2]
		sym = p[12]
		ID = p[1]
		chrom_to_ai[chrom].append(expra.ID_to_index[ID])
	
	
	for s_c57, s_cast in zip(expra.samples[::2], expra.samples[1::2]):
		sample = s_c57.rsplit('_',1)[0]
		for chrom in chrom_to_ai:
			for samplepart, chromosome in o.removal:
				if chrom==chromosome and samplepart in sample:
					continue
			for samplepart, chromosome in o.addition:
				if chrom==chromosome and samplepart in sample:
					print sample, chrom
					continue
			c57count, castcount, bicount = 0, 0, 0
			for ai in chrom_to_ai[chrom]:
				if expra[s_c57][ai] and not expra[s_cast][ai]:
					c57count += 1
				elif not expra[s_c57][ai] and expra[s_cast][ai]:
					castcount += 1
				elif expra[s_c57][ai] and expra[s_cast][ai]:
					bicount += 1
			totalcount = c57count + castcount + bicount
			
			if totalcount and abs(c57count-castcount)/totalcount > o.threshold:
				print sample, chrom, c57count/totalcount, castcount/totalcount

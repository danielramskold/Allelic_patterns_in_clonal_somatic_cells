from __future__ import division
import argparse, dr_tools, numpy

base_to_index = {'A':2, 'C':3, 'G':4, 'T':5}

class SNPinfo:
	def __init__(self, c57base, castbase):
		self.c57index = base_to_index[c57base]
		self.castindex = base_to_index[castbase]
		self.bases = '%s\t%s'%(c57base, castbase)
		self.c57count = 0
		self.castcount = 0

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-c', '--cellsums_files', metavar='cellsums', required=True, nargs='+')
	parser.add_argument('-v', '--snp_vcf', help='e.g. castsnp_alleles.txt')
	parser.add_argument('-o', '--outfile', default='/dev/null')
	parser.add_argument('-f', '--figure')
	parser.add_argument('-R', '--minratio', default=0, type=float)
	parser.add_argument('-rt', '--minreads_sum', default=1, type=int)
	parser.add_argument('-ra', '--minreads_allele', default=0, type=int)
	parser.add_argument('-V', '--snp_validatedbefore', help='e.g. validated_cast_c57_snps.txt')
	o = parser.parse_args()
	
	if o.snp_validatedbefore:
		allowed_coord = set()
		for p in dr_tools.splitlines(o.snp_validatedbefore, ignore='#'):
			chromosome = p[0]
			position = p[1]
			coord = '%s\t%s'%(chromosome, position)
			allowed_coord.add(coord)
	else: allowed_coord = None
	
	database_snps = dict()
	for p in dr_tools.splitlines(o.snp_vcf, ignore='#'):
		chromosome = 'chr'+p[0]
		position = p[1]
		c57base = p[3]
		castbase = p[4]
		if ',' in c57base or ',' in castbase: continue
		coord = '%s\t%s'%(chromosome, position)
		if allowed_coord is not None and coord not in allowed_coord: continue
		database_snps[coord] = SNPinfo(c57base, castbase)
	
	for filepath in o.cellsums_files:
		for p in dr_tools.splitlines(filepath):
			coord = '%s\t%s'%(p[0], p[1])
			if coord not in database_snps:
				# strange, it should be (in database_snps)
				# unless o.snp_validatedbefore was used
				continue
			snp = database_snps[coord]
			snp.c57count += int(p[snp.c57index])
			snp.castcount += int(p[snp.castindex])
	
	ratios = []
	with open(o.outfile, 'w') as outfh:
		for coord, snpinfo in database_snps.items():
			reads = snpinfo.c57count+snpinfo.castcount
			if reads == 0: ratio = 0
			else: ratio = snpinfo.c57count/reads
			if o.minratio <= ratio <= (1-o.minratio) and reads >= o.minreads_sum and min(snpinfo.c57count, snpinfo.castcount) >= o.minreads_allele:
				print >>outfh, dr_tools.join(coord, snpinfo.bases, '0', '1.00', '1.00', '1.00', '1.00', '1.00')
				ratios.append(ratio)
	
	if o.figure:
		import pylab
		step = 0.005
		xarr, yarr = dr_tools.bin(ratios, -step, 1+step, step, 1)
		#yarr = [y/len(ratios) for y in yarr]
		pylab.plot(xarr, yarr, 'k-')
		pylab.savefig(o.figure)

from __future__ import division
import argparse, pylab, dr_tools, math
from collections import defaultdict

def MAfraction(expra, sample, genelist):
	count_bi, count_mono = 0, 0
	for sym, c57, cast in zip(expra['symbols'], expra[sample+'_c57only'], expra[sample+'_castonly']):
		if sym in genelist:
			if c57 and cast: count_bi += 1
			elif c57 or cast: count_mono += 1
	return count_mono/(count_bi + count_mono)

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('rpkmfile')
	parser.add_argument('nondiatable')
	parser.add_argument('allelehits')
	parser.add_argument('genelist')
	o = parser.parse_args()
	
	genelist = set(dr_tools.loadlist(o.genelist))
	expra = dr_tools.loadexpr(o.allelehits, True)
	expr = dr_tools.loadexpr(o.rpkmfile, True)
	spikes_i = [i for i, ID in enumerate(expr['IDs']) if 'ERCC' in ID]
	
	xarr = defaultdict(list)
	yarr = defaultdict(list)
	
	for p in dr_tools.splitlines(o.nondiatable):
		if p[0] == '#sample':
			index_cellsource = p.index('cell.type')
		else:
			sample = p[0]
			try:
				ERCC_readsum = sum(expr[sample][spike] for spike in spikes_i)
				if ERCC_readsum < 100: continue
				sample_i = expr.samples.index(sample)
				mRNA_readsum = expr.normalizationreads[sample_i]
				cellsource = p[index_cellsource]
				
				yarr[cellsource].append(1-MAfraction(expra, sample, genelist))
				xarr[cellsource].append(math.log(mRNA_readsum/ERCC_readsum, 10))
				
				if xarr[cellsource][-1] >1.5 and yarr[cellsource][-1] < 0.6:
					print xarr[cellsource][-1], yarr[cellsource][-1], sample, ERCC_readsum, mRNA_readsum
			except KeyError:
				continue

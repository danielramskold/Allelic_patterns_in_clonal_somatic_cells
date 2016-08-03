from __future__ import division
import argparse, pylab, dr_tools
from scipy import stats

def MAfraction(expra, sample):
	count_bi, count_mono = 0, 0
	for sym, c57, cast in zip(expra['symbols'], expra[sample+'_c57only'], expra[sample+'_castonly']):
		if c57 and cast: count_bi += 1
		elif c57 or cast: count_mono += 1
	return count_mono/(count_bi + count_mono)

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('rpkmfile')
	parser.add_argument('diatable')
	parser.add_argument('allelehits')
	o = parser.parse_args()
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	expr = dr_tools.loadexpr(o.rpkmfile, True)
	spikes_i = [i for i, ID in enumerate(expr['IDs']) if 'ERCC' in ID]
	
	xarr = []
	yarr = []
	
	for p in dr_tools.splitlines(o.diatable):
		if p[0] == '#sample':
			index_dia = [p.index('cytoplasm.length'), p.index('cytoplasm.width')]
		else:
			sample = p[0]
			ERCC_readsum = sum(expr[sample][spike] for spike in spikes_i)
			sample_i = expr.samples.index(sample)
			mRNA_readsum = expr.normalizationreads[sample_i]
			try:
				width = float(p[index_dia[1]])
				length = float(p[index_dia[0]])
			except ValueError:
				continue
			
			xarr.append(1-MAfraction(expra, sample))
			yarr.append(mRNA_readsum/ERCC_readsum) # like D Edsgards' axis
	
	print stats.pearsonr(xarr, yarr)
	
	pylab.plot(xarr, yarr, '.')
	pylab.savefig('ERCCsum_vs_biallelic.pdf')

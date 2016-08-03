from __future__ import division
import argparse, pylab, dr_tools
from scipy import stats

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('rpkmfile')
	parser.add_argument('diatable')
	parser.add_argument('--dim', type=float, default=3)
	o = parser.parse_args()
	
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
			xarr.append((width*length)**(o.dim/2))
			yarr.append(mRNA_readsum/ERCC_readsum) # like D Edsgards' axis
	
	print stats.pearsonr(xarr, yarr)
	
	pylab.plot(xarr, yarr, '.')
	pylab.savefig('ERCCsum_vs_diameter.pdf')

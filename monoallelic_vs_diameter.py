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
	parser.add_argument('allelehits')
	parser.add_argument('diatable')
	parser.add_argument('--dim', type=float, default=3)
	o = parser.parse_args()
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	
	xarr = []
	yarr = []
	
	for p in dr_tools.splitlines(o.diatable):
		if p[0] == '#sample':
			index_dia = [p.index('cytoplasm.length'), p.index('cytoplasm.width')]
		else:
			sample = p[0]
			
			try:
				width = float(p[index_dia[1]])
				length = float(p[index_dia[0]])
			except ValueError:
				continue
			xarr.append((width*length)**(o.dim/2))
			yarr.append(MAfraction(expra, sample))
	
	print stats.pearsonr(xarr, yarr)
	
	pylab.plot(xarr, yarr, '.')
	pylab.savefig('monoallelic_vs_diameter.pdf')

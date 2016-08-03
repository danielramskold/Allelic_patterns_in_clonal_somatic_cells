from __future__ import division
import argparse, pylab, dr_tools, math, numpy
from collections import defaultdict

def MAfraction(expra, sample, genelist):
	count_bi, count_mono = 0, 0
	for sym, c57, cast in zip(expra['symbols'], expra[sample+'_c57only'], expra[sample+'_castonly']):
		if sym in genelist and expr[sample][expr.symbol_to_index[sym]] >= 20: # this is the same way to apply cutoff as in plot_monoallelic_by_cell_minrpkm.py
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
	expr = dr_tools.loadexpr(o.rpkmfile, False)
	spikes_i = [i for i, ID in enumerate(expr['IDs']) if 'ERCC' in ID]
	
	xarr = defaultdict(list)
	yarr = defaultdict(list)
	
	for p in dr_tools.splitlines(o.nondiatable):
		if p[0] == '#sample':
			index_cellsource = p.index('cell.type')
		else:
			sample = p[0]
			if sample == 'BQx46_indD_EmbryoMEF_BxC': continue # degraded sample
			try:
				ERCC_readsum = sum(expr[sample][spike] for spike in spikes_i)
				if ERCC_readsum < 100: continue
				sample_i = expr.samples.index(sample)
				mRNA_readsum = expr.normalizationreads[sample_i]
				cellsource = p[index_cellsource]
				if cellsource in ('mef', 'MAF'): cellsource='fibr'
			
				yarr[cellsource].append(1-MAfraction(expra, sample, genelist))
				xarr[cellsource].append(math.log(mRNA_readsum/ERCC_readsum, 2))
				print sample, cellsource, xarr[cellsource][-1], yarr[cellsource][-1]
			except KeyError:
				#print sample
				continue
	
	
	from scipy import stats
	from scipy import optimize
	def fitting_f(xvalues, *params):
		yvalues = []
		for x in xvalues:
			y = 1/(1+params[0]*2**(-params[1]*x))
			yvalues.append(y)
		return numpy.array(yvalues)
	
	X,Y = [], []
	for cellsource in xarr:
		X.extend(xarr[cellsource])
		Y.extend(yarr[cellsource])
		print cellsource, 'n=', len(xarr[cellsource])
	print 'spearman', stats.spearmanr(X, Y)
	coefficients_p, pcov = optimize.curve_fit(fitting_f, X, Y, [1,1])
	print 'line fit coefficients', coefficients_p
	x_line = numpy.arange(2, 11, 0.1)
	pylab.plot(x_line, fitting_f(x_line, *coefficients_p), 'k-')
	
	
	
	for cellsource in xarr:
		pylab.plot(xarr[cellsource], yarr[cellsource], 'o', label=cellsource, fillstyle='none', markeredgecolor=dr_tools.randomcolour(), markersize=5, markeredgewidth=2)
	pylab.legend(loc='lower right')
	pylab.ylabel('biallelic fraction (>=20 rpkm)')
	pylab.xlabel('mRNA amount, log2 (1(nonlog)=same as ERCC amount)')
	pylab.savefig('ERCCsum_vs_biallelic_nonDia_min20rpkm_v5.pdf')
	
	

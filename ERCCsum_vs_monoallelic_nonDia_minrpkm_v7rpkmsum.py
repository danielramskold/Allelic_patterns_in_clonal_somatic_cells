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

def calc_ERCC_moleculenumber(tablefile, before_dilution_vol_ul):
	Mix1_i = 3
	conc_attomolul = 0
	attomol = 602214.12927
	for i, p in enumerate(dr_tools.splitlines(tablefile)):
		if i == 0:
			if not 'attomoles/ul' in p[Mix1_i]: raise Exception
			#if not 'Mix 1' in p[Mix1_i]: raise Exception
		else:
			conc_attomolul += float(p[Mix1_i])
	return conc_attomolul * before_dilution_vol_ul * 602214.12927

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('rpkmfile')
	parser.add_argument('nondiatable')
	parser.add_argument('allelehits')
	parser.add_argument('genelist')
	o = parser.parse_args()
	
	ERCCvol_ul = 4e-7
	ERCC_moleculenumber = calc_ERCC_moleculenumber('ERCC.txt', ERCCvol_ul)
	
	genelist = set(dr_tools.loadlist(o.genelist))
	expra = dr_tools.loadexpr(o.allelehits, True)
	expr = dr_tools.loadexpr(o.rpkmfile, False)
	spikes_i = [i for i, ID in enumerate(expr['IDs']) if 'ERCC' in ID]
	genes_i = [i for i, ID in enumerate(expr['IDs']) if 'ERCC' not in ID]
	
	xarr = defaultdict(list)
	yarr = defaultdict(list)
	
	for p in dr_tools.splitlines(o.nondiatable):
		if p[0] == '#sample':
			index_cellsource = p.index('cell.type')
		else:
			sample = p[0]
			if sample == 'BQx46_indD_EmbryoMEF_BxC': continue # degraded sample
			try:
				ERCC_rpkmsum = sum(expr[sample][spike] for spike in spikes_i)
				if ERCC_rpkmsum < 100: continue
				sample_i = expr.samples.index(sample)
				mRNA_rpkmsum = sum(expr[sample][gene] for gene in genes_i)
				cellsource = p[index_cellsource]
				if cellsource in ('mef', 'MAF'): cellsource='fibroblast'
				if cellsource == 'brain': cellsource = 'neuron'
				if sample in ('bqb51_indE_3_liver_BxC', 'bqb52_indE_4_liver_BxC'): cellsource = 'liver_nonhep' # samples missing hepatocyte markers
				elif cellsource == 'liver': cellsource = 'hepatocyte'
				
				yarr[cellsource].append(1-MAfraction(expra, sample, genelist))
				xarr[cellsource].append(math.log(mRNA_rpkmsum/ERCC_rpkmsum * ERCC_moleculenumber, 2))
				print sample, cellsource, xarr[cellsource][-1], yarr[cellsource][-1]
			except KeyError:
				print 'missing', sample
				continue
	
	
	from scipy import stats
	from scipy import optimize
	def fitting_f(xvalues, *params):
		yvalues = []
		for x in xvalues:
			y = 1/(1+ 2.718281**(params[1]*x+params[0]))
			yvalues.append(y)
		return numpy.array(yvalues)
	
	X,Y = [], []
	for cellsource in xarr:
		X.extend(xarr[cellsource])
		Y.extend(yarr[cellsource])
		print cellsource, 'n=', len(xarr[cellsource])
	print 'spearman', stats.spearmanr(X, Y)
	coefficients_p, pcov = optimize.curve_fit(fitting_f, X, Y, [0,-1])
	print 'line fit coefficients', coefficients_p
	x_line = numpy.arange(math.floor(min(X)), math.ceil(max(X)), 0.1)
	pylab.plot(x_line, fitting_f(x_line, *coefficients_p), 'k-')
	
	
	import random; random.seed(10)
	for cellsource in xarr:
		pylab.plot(xarr[cellsource], yarr[cellsource], 'o', label=cellsource, fillstyle='none', markeredgecolor=dr_tools.randomcolour(), markersize=5, markeredgewidth=2)
	pylab.legend(loc='lower right')
	pylab.ylabel('biallelic fraction (>=20 rpkm)')
	pylab.xlabel('polyA+ RNA molecules, log2')
	pylab.savefig('ERCCsum_vs_biallelic_nonDia_min20rpkm_v7.pdf')
	
	

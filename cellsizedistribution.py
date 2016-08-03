from __future__ import division
import argparse, pylab, dr_tools, math, numpy
from collections import defaultdict

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
	o = parser.parse_args()
	
	ERCCvol_ul = 4e-7
	ERCC_moleculenumber = calc_ERCC_moleculenumber('ERCC.txt', ERCCvol_ul)
	
	expr = dr_tools.loadexpr(o.rpkmfile, False)
	spikes_i = [i for i, ID in enumerate(expr['IDs']) if 'ERCC' in ID]
	genes_i = [i for i, ID in enumerate(expr['IDs']) if 'ERCC' not in ID]
	
	xarr = defaultdict(list)
	
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
				
				xarr[cellsource].append(math.log(mRNA_rpkmsum/ERCC_rpkmsum * ERCC_moleculenumber, 2))
			except KeyError:
				print 'missing', sample
				continue
	
	
	X = xarr['fibroblast']
	#xarr_bin, yarr_bin = dr_tools.bin(X, math.floor(min(X)), math.ceil(max(X)), 0.5)
	pylab.hist(X, 20)
	pylab.xlabel('polyA+ RNA molecules, log2')
	pylab.ylabel('fibroblasts')
	pylab.savefig('cellsizedistribution_fibroblast.pdf')
	

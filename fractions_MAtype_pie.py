from __future__ import division
import argparse, pylab, dr_tools

def list_mono(samples, randomize_per_gene=False, use_ti=None):
	c57mono, castmono, bi = set(), set(), set()
	assert exprt['symbols'] == expra['symbols']
	ti_used = set()
	for ti, sym in enumerate(exprt['symbols']):
		if use_ti is not None:
			if ti not in use_ti: continue
		else:
			meanexpr = numpy.mean([exprt[s][ti] for s in samples])
			if meanexpr < o.minrpkm: continue
			if o.maxrpkm is not None and meanexpr >= o.maxrpkm: continue
			if allowedgenes and sym not in allowedgenes: continue
		ai = ti
		
		c57 = 0
		cast = 0
		for s in samples:
			if 'BxC' in s: swap = True
			elif 'CxB' in s: swap = False
			else: raise Exception
			
			if randomize_per_gene and random.random() < 0.5:
				swap = not swap
			
			if swap:
				c57 += bool(expra[s+'_castonly'][ai])
				cast += bool(expra[s+'_c57only'][ai])
			else:
				cast += bool(expra[s+'_castonly'][ai])
				c57 += bool(expra[s+'_c57only'][ai])
		if c57 and cast: bi.add(sym)
		elif c57: c57mono.add(sym)
		elif cast: castmono.add(sym)
		ti_used.add(ti)
	return c57mono, castmono, bi, ti_used

class MA_classes:
	def __init__(self):
		self.randomdynamic = 0
		self.randominherited = 0
		self.imprinted = 0
		self.chrX = 0

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('rpkmfile')
	parser.add_argument('allelehits')
	parser.add_argument('-X', '--chrX_genes', required=True)
	parser.add_argument('-I', '--imprinted_genes', required=True)
	parser.add_argument('-A', '--autosomal_genes_all', required=True)
	parser.add_argument('-S', '--autosomal_genes_selection')
	parser.add_argument('-m', '--minrpkm', type=float, default=20)
	parser.add_argument('-s', '--samplelist_clone', nargs='+', required=True)
	o = parser.parse_args()
	
	o.maxexpr = None
	if not o.autosomal_genes_selection:
		o.autosomal_genes_selection = o.autosomal_genes_all
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	exprt = dr_tools.loadexpr(o.rpkms, False)
	
	allowedgenes = set(dr_tools.loadlist(o.autosomal_genes_selection)) - set(dr_tools.loadlist(o.imprinted_genes))
	
	for samplelistfile in o.samplelist_clone:
		samples = dr_tools.loadlist(samplelistfile)
		c57mono, castmono, clonal_bi, ti_used = list_mono(samples)
		clonal_mono = c57mono | castmono
		
		num_random_mono = []
		for ri in range(10):
			c57mono, castmono, bi, ti_used = list_mono(samples, True, ti_used)
			num_random_mono.append(len(c57mono)+len(castmono))
		exp_random_mono = numpy.mean(num_random_mono)

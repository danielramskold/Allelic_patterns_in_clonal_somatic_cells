from __future__ import division
import argparse, dr_tools, numpy, pylab, random

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
	return c57mono, castmono, bi

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allelehits', required=True)
	parser.add_argument('-r', '--rpkms', required=True)
	parser.add_argument('-m', '--minrpkm', default=20, type=float)
	parser.add_argument('-M', '--maxrpkm', type=float)
	parser.add_argument('-gi', '--allowedgenes')
	parser.add_argument('--random_seed', type=int)
	parser.add_argument('-om', '--outputmode', choices=['stats', 'c57', 'cast', 'bi'], default='stats')
	o = parser.parse_args()
	
	if o.random_seed is not None: random.seed(o.random_seed)
	
	allowedgenes = set(dr_tools.loadlist(o.allowedgenes)) if o.allowedgenes else None
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	exprt = dr_tools.loadexpr(o.rpkms, False)
	
	c57mono, castmono, bi = list_mono(exprt.samples)
	if o.outputmode == 'stats':
		tot_genes = len(bi)+len(c57mono)+ len(castmono)
		print len(castmono), len(c57mono), (len(c57mono)+len(castmono))/tot_genes
		c57mono, castmono, bi = list_mono(exprt.samples, True)
		tot_genes = len(bi)+len(c57mono)+ len(castmono)
		print len(castmono), len(c57mono), (len(c57mono)+len(castmono))/tot_genes
	else:
		for sym in {'c57':c57mono, 'cast':castmono, 'bi':bi}[o.outputmode]:
			print sym

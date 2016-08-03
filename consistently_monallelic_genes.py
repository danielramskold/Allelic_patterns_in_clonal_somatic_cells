import argparse, dr_tools, random, numpy

def count_mono(samples, randomize_per_gene=False):
	c57mono_genes, castmono_genes = set(),set()
	assert exprt['symbols'] == expra['symbols']
	ti_used = set()
	for ti, sym in enumerate(exprt['symbols']):
		if numpy.mean([exprt[s][ti] for s in samples]) < o.minrpkm: continue
		if allowedgenes and sym not in allowedgenes: continue
		ai = ti
		
		c57 = 0
		cast = 0
		for s in samples:
			if randomize_per_gene and random.random() < 0.5:
				c57 += expra[s+'_castonly'][ai]
				cast += expra[s+'_c57only'][ai]
			else:
				cast += expra[s+'_castonly'][ai]
				c57 += expra[s+'_c57only'][ai]
		if c57 and cast: pass
		elif c57: c57mono_genes.add(sym)
		elif cast: castmono_genes.add(sym)
	return c57mono_genes, castmono_genes

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allelehits', required=True)
	parser.add_argument('-r', '--rpkms', required=True)
	parser.add_argument('-m', '--minrpkm', default=20, type=float)
	parser.add_argument('-gi', '--allowedgenes')
	parser.add_argument('-rg', '--randomize_per_gene', action='store_true')
	parser.add_argument('-G', '--genotype', choices=['cast', 'c57'], required=True)
	o = parser.parse_args()
	
	random.seed(20)
	
	allowedgenes = set(dr_tools.loadlist(o.allowedgenes)) if o.allowedgenes else None
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	exprt = dr_tools.loadexpr(o.rpkms, False)
	
	c57mono_genes, castmono_genes = count_mono(exprt.samples, o.randomize_per_gene)
	genes = c57mono_genes if o.genotype == 'c57' else castmono_genes
	for genesym in genes:
		print genesym

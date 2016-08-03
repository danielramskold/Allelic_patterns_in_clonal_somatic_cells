from __future__ import division
import argparse, dr_tools, numpy, pylab, random

def count_mono(samples, randomize_per_gene=False, use_ti=None):
	c57mono, castmono, bi = 0,0,0
	assert exprt['symbols'] == expra['symbols']
	ti_used = set()
	for ti, sym in enumerate(exprt['symbols']):
		if use_ti is not None:
			if ti not in use_ti: continue
		else:
			meanexpr = numpy.mean([exprt[s][ti] for s in samples])
			if meanexpr < o.minrpkm: continue
			if o.maxrpkm is not None and meanexpr >= o.maxrpkm: continue
			if disallowedgenes and sym in disallowedgenes: continue
			if allowedgenes and sym not in allowedgenes: continue
		ai = ti
		
		c57fraction = 0.5
		c57 = 0
		cast = 0
		both = 0
		for s in samples:
			a1 = bool(expra[s+'_c57only'][ai])
			a2 = bool(expra[s+'_castonly'][ai])
			if a1 and a2: both += 1
			elif a1 and not a2:
				if randomize_per_gene and random.random() < c57fraction:
					cast += 1
				else:
					c57 += 1
			elif a2 and not a1:
				if randomize_per_gene and random.random() > c57fraction:
					c57 += 1
				else:
					cast += 1
		
		tot_cells = max(1, c57+cast+both) # don't count unexpressing cells
		if c57/tot_cells >= o.min_cell_fraction:
			c57mono += 1
		elif cast/tot_cells >= o.min_cell_fraction:
			castmono += 1
		elif both: bi += 1
		
		ti_used.add(ti)
	tot_genes = c57mono + castmono + bi
	if tot_genes == 0 or o.genecount: tot_genes = 1
	return c57mono/tot_genes, castmono/tot_genes, bi/tot_genes, ti_used

class Bars:
	def __init__(self):
		self.xarr = []
		self.x_pos = 0
		self.mono = []
		self.cast = []
		self.tot = []
		self.labels = []
		self.labels_x = []
	
	def addbar(self, label, c57mono, castmono, bi=None, ti_used=None):
		self.cast.append(castmono)
		self.mono.append(castmono+c57mono)
		self.tot.append(c57mono+castmono+bi if o.genecount else 1)
		self.xarr.append(self.x_pos)
		if label:
			self.labels.append(label)
			self.labels_x.append(self.x_pos)
		print self.x_pos, label, c57mono, castmono
		self.x_pos += 1

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allelehits', required=True)
	parser.add_argument('-r', '--rpkms', required=True)
	parser.add_argument('-o', '--figure', required=True)
	parser.add_argument('-m', '--minrpkm', default=20, type=float)
	parser.add_argument('-M', '--maxrpkm', type=float)
	parser.add_argument('-c', '--clonal_groups', default='clonal_groups.txt')
	parser.add_argument('-gi', '--allowedgenes')
	parser.add_argument('-ge', '--disallowedgenes', nargs='+')
	parser.add_argument('--addgreen', action='store_true')
	parser.add_argument('--genecount', action='store_true')
	parser.add_argument('--random_seed', type=int)
	parser.add_argument('-R', '--random_bars', type=int, default=1)
	parser.add_argument('-f', '--min_cell_fraction', type=float, default=1.0)
	o = parser.parse_args()
	
	if o.random_seed is not None: random.seed(o.random_seed)
	
	allowedgenes = set(dr_tools.loadlist(o.allowedgenes)) if o.allowedgenes else None
	if o.disallowedgenes:
		disallowedgenes = set()
		for filename in o.disallowedgenes:
			disallowedgenes.update(set(dr_tools.loadlist(filename)))
	else:
		disallowedgenes = None
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	exprt = dr_tools.loadexpr(o.rpkms, False)
	
	bars = Bars()
	for clonal_group in dr_tools.loadlist(o.clonal_groups):
		samples = [s for s in exprt.samples if any(s.startswith(clonal_group_start) or s.startswith('pool.'+clonal_group_start) for clonal_group_start in clonal_group.split('\t'))]
		merged_output = count_mono(samples, False, None)
		use_ti = merged_output[-1]
		for sample in samples:
			bars.addbar('', *count_mono([sample], False, use_ti))
		bars.x_pos += 1
		bars.addbar(clonal_group, *merged_output)
		for ri in range(o.random_bars):
			bars.addbar('r', *count_mono(samples, True, use_ti))
		bars.x_pos += 2
	
	if o.addgreen: pylab.bar(bars.xarr, bars.tot, facecolor='#009933', linewidth=0, width=0.9)
	pylab.bar(bars.xarr, bars.mono, facecolor='#888800', linewidth=0, width=0.9)
	pylab.bar(bars.xarr, bars.cast, facecolor='#aa8800', linewidth=0, width=0.9)
	pylab.xticks(bars.labels_x, bars.labels, rotation=90)
	pylab.subplots_adjust(bottom=0.3)
	pylab.savefig(o.figure)

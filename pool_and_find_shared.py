from __future__ import division
import argparse, dr_tools, numpy, pylab, random
from collections import defaultdict

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

class Genecat:
	def __init__(self, label):
		self.label = label
		self.D = defaultdict(set)
	
	def addgenes(self, group, *genesets):
		for geneset in genesets:
			for gene in geneset:
				#print gene, geneset
				self.addgene(gene, group)
	
	def addgene(self, gene, group):
		self.D[gene].add(group)
	
	def print_counts(self):
		print self.label
		groupcounts = defaultdict(list)
		for gene in self.D:
			groupcounts[len(self.D[gene])].append(gene)
		for count in sorted(groupcounts):
			if len(groupcounts[count]) <= 5:
				print count, len(groupcounts[count]), ' '.join(groupcounts[count])
			else:
				print count, len(groupcounts[count])

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allelehits', required=True)
	parser.add_argument('-r', '--rpkms', required=True)
	parser.add_argument('-m', '--minrpkm', default=20, type=float)
	parser.add_argument('-M', '--maxrpkm', type=float)
	parser.add_argument('-c', '--clonal_groups', default='clonal_groups.txt')
	parser.add_argument('-gi', '--allowedgenes')
	parser.add_argument('-ge', '--disallowedgenes')
	parser.add_argument('--random_seed', type=int)
	o = parser.parse_args()
	
	if o.random_seed is not None: random.seed(o.random_seed)
	
	allowedgenes = set(dr_tools.loadlist(o.allowedgenes)) if o.allowedgenes else None
	disallowedgenes = set(dr_tools.loadlist(o.disallowedgenes)) if o.disallowedgenes else None
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	exprt = dr_tools.loadexpr(o.rpkms, False)
	
	randomC = Genecat('random')
	obsC = Genecat('observed')
	
	for clonal_group in dr_tools.loadlist(o.clonal_groups):
		samples = [s for s in exprt.samples if any(s.startswith(clonal_group_start) or s.startswith('pool.'+clonal_group_start) for clonal_group_start in clonal_group.split('\t'))]
		merged_output = list_mono(samples, False, None)
		obsC.addgenes(clonal_group, *merged_output[:2])
		use_ti = merged_output[-1]
		control_output = list_mono(samples, True, use_ti)
		randomC.addgenes(clonal_group, *control_output[:2])
	
	obsC.print_counts()
	randomC.print_counts()

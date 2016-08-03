from __future__ import division
import argparse, pylab, dr_tools, random, itertools

def pick_samples(samples_nom):
	samples = []
	for line in samples_nom:
		if ',' in line:
			while True:
				s = random.choice(line.split(','))
				if not s+'_c57only' in expra:
					print s
				else: break
				if s in multi_cell_samples:
					print 'multicellsample', s
			samples.append(s)
		else:
			samples.append(line)
	return samples

def monoallelic_fraction_pooled(samples, expra, allowed_genes, disallowed_genes, randomize_per_gene=False):
	monocount = 0
	bicount = 0
	for ai, sym in enumerate(expra['symbols']):
		if allowed_genes is not None and sym not in allowed_genes: continue
		if sym in disallowed_genes: continue
		if randomize_per_gene:
			if o.allelerand_skew:
				identities = ['c57' if expra[s+'_c57only'][ai]>0 and expra[s+'_castonly'][ai]==0 else 'cast' if expra[s+'_c57only'][ai]==0 and expra[s+'_castonly'][ai]>0 else '0' for s in samples]
				nonswapped = [True if I=='0' else random.random() < allelerand_skew[ai] if I=='c57' else 1-allelerand_skew[ai] for I in identities]
			else:
				nonswapped = [random.random() < 0.5 for s in samples]
			c57 = sum((expra[s+'_c57only'][ai]>0) if nonswap else (expra[s+'_castonly'][ai]>0) for s,nonswap in zip(samples, nonswapped))
			cast = sum((expra[s+'_castonly'][ai]>0) if nonswapped else (expra[s+'_c57only'][ai]>0) for s,nonswap in zip(samples, nonswapped))
			#print c57, cast, sum(expra[s+'_c57only'][ai]>0 for s in samples)
		else:
			c57 = sum(expra[s+'_c57only'][ai]>0 for s in samples)
			cast = sum(expra[s+'_castonly'][ai]>0 for s in samples)
		if c57 > 0 and cast > 0:
			bicount += 1
		elif c57 > 0 or cast > 0:
			monocount += 1
	return monocount/(monocount+bicount)

multi_cell_samples = {
"MAF_CxB_clone_A_10": 2, # two cells attached together
"MAF_CxB_clone_B_19": 2, # two cells attached together
"MAF_CxB_clone_B_31": 2, # two cells attached together
"MAF_CxB_clone_B_41": 2, # two small cells attached together
"MAF_CxB_clone_B_53": 2, # two small cells attached together
"MAF_CxB_clone_B_61": 2, # two small cells attached together
"MAF_CxB_clone_B_65": 2, # two small cells attached together
"MAF_CxB_clone_B_69": 2, # two small cells attached together. overblown
"MAF_CxB_clone_B_75": 2, # two small cells attached together
"MAF_CxB_clone_B_76": 2, # two cells attached together
"MAF_CxB_clone_B_87": 3, # big clump. 3? cells
"MAF_CxB_clone_B_88": 2, # two cells attached together
"MAF_CxB_clone_B_91": 2, # two small cells attached together
"MAF_CxB_clone_B_100": 2, # two small cells attached together
"MEF_E14_Clone_BxC_H_10": 2, # 2 cells picked
"MEF_E14_Clone_BxC_H_12": 2, # 2 cells picked
"MEF_E14_Clone_BxC_H_14": 2, # 2 cells picked
"pool.MAF_CxB_clone_B_58": 2, # two flat cells together
"pool.MAF_CxB_clone_B_77": 3, # clump of flat cells 2-4
"pool.MAF_CxB_clone_B_86": 3, # big dark clump. 3? cells
"pool.MAF_CxB_clone_B_95": 15, # 15 small cells in one go
"pool.MAF_CxB_clone_B_101": 3, # 3-4 cells in one go
"pool.MAF_CxB_RE_clone_B_102": 5, # 5 cells in one go. one of them is big
"pool.MAF_CxB_RE_clone_B_103": 3, # 3? big cells attached
"pool.MAF_CxB_RE_clone_B_104": 5, # 5 cells. one is big
"pool.MEF_E14_Clone_BxC_H_4": 2, # 2 cells picked
}

def ratio(expra, gi, samples):
	num_c57only = sum(expra[s+'_c57only'][gi] > 0 and expra[s+'_castonly'][gi] == 0 for s in samples)
	num_castonly = sum(expra[s+'_c57only'][gi] == 0 and expra[s+'_castonly'][gi] > 0 for s in samples)
	num_both = num_c57only+num_castonly
	if num_both == 0: return 0.5
	return num_c57only/num_both

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allelehits', required=True)
	parser.add_argument('-gi', '--allowedgenes')
	parser.add_argument('-ge', '--disallowedgenes', nargs='+', default=[])
	parser.add_argument('-R', '--random_dots', type=int, default=1)
	parser.add_argument('-s', '--samplelist', required=True, nargs='+')
	parser.add_argument('-n', default=4, type=int)
	parser.add_argument('-o', '--figure', default='poolN.pdf')
	parser.add_argument('-S', '--subtract_allelerand', action='store_true')
	parser.add_argument('-r', '--allelerand_skew', action='store_true')
	o = parser.parse_args()
	
	expra = dr_tools.loadexpr(o.allelehits, True)
	
	if o.allowedgenes:
		allowed_genes = set(dr_tools.loadlist(o.allowedgenes))
	else:
		allowed_genes = None
	disallowed_genes = set()
	for filename in o.disallowedgenes:
		disallowed_genes.update(set(dr_tools.loadlist(filename)))
	
	random.seed(0)
	
	samples_n = dict((samplelist, [random.sample(dr_tools.loadlist(samplelist, ignore='#'), o.n) for di in range(o.random_dots)]) for samplelist in o.samplelist)
	
	samples_all = [sa.split('_c57only')[0] for sa in expra.samples[::2]]
	allelerand_skew =  dict((gi, ratio(expra, gi, samples_all)) for gi in range(len(expra['symbols'])))
	
	n_cells = o.n
	
	# sim
	boxplot_y = []
	boxplot_x = []
	labels_x = []
	for x, samplelist in enumerate(o.samplelist):
		monofractions = []
		for di in range(o.random_dots):
			samples_l = pick_samples(samples_n[samplelist][di])
			if o.subtract_allelerand:
				monofractions.append(monoallelic_fraction_pooled(samples_l, expra, allowed_genes, disallowed_genes, False)-monoallelic_fraction_pooled(samples_l, expra, allowed_genes, disallowed_genes, True))
			else:
				monofractions.append(monoallelic_fraction_pooled(samples_l, expra, allowed_genes, disallowed_genes, False))
			
		boxplot_y.append(monofractions)
		boxplot_x.append(x)
		labels_x.append(samplelist.split('/')[-1].split('.txt')[0].split('.csv')[0])
	pylab.boxplot(boxplot_y, positions=boxplot_x, whis=10000)
	pylab.xticks(boxplot_x, labels_x, rotation=10)
	
	pylab.savefig(o.figure)
	
	from scipy import stats
	for i1, i2 in itertools.combinations(range(len(boxplot_y)), 2):
		print i1, i2, stats.ranksums(boxplot_y[i1], boxplot_y[i2])

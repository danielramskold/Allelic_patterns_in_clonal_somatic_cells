from __future__ import division
import argparse, dr_tools, numpy, pylab, random

def count_mono(samples, randomize_per_gene=False, use_ti=None, min_cell_fraction=1.0):
	c57mono, castmono, bi = 0,0,0
	ti_used = set()
	for ti, sym in enumerate(expra['symbols']):
		if use_ti is not None:
			if ti not in use_ti: continue
		else:
			if disallowedgenes and sym in disallowedgenes: continue
			if allowedgenes and sym not in allowedgenes: continue
		ai = ti
		
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
		
		#if (c57 and cast) or both: bi+=1
		#elif c57: c57mono += 1
		#elif cast: castmono += 1
		
		tot_cells = max(1, c57+cast+both) # don't count unexpressing cells
		if c57/tot_cells >= min_cell_fraction:
			c57mono += 1
		elif cast/tot_cells >= min_cell_fraction:
			castmono += 1
		elif both: bi += 1
		
		
		ti_used.add(ti)
	tot_genes = c57mono + castmono + bi
	if tot_genes == 0 or o.genecount: tot_genes = 1
	return c57mono/tot_genes, castmono/tot_genes, bi/tot_genes, ti_used

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

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--allelehits', required=True)
	parser.add_argument('-a2', '--allelehits_multicell')
	parser.add_argument('-c', '--clonal_groups', default='clonal_groups.txt')
	parser.add_argument('-gi', '--allowedgenes')
	parser.add_argument('-ge', '--disallowedgenes', nargs='+')
	parser.add_argument('--genecount', action='store_true')
	parser.add_argument('--random_seed', type=int)
	parser.add_argument('-R', '--random_dots', type=int, default=1)
	parser.add_argument('-C', '--clonal_group', default='clone_B')
	parser.add_argument('-n1', '--start_n', default=1, type=int)
	parser.add_argument('-n2', '--end_n', default=15, type=int)
	parser.add_argument('-o', '--figure', default='pool_n.pdf')
	parser.add_argument('--nonrandom_n1', action='store_true')
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
	
	c57fraction = 0.5
	
	for clonal_group in dr_tools.loadlist(o.clonal_groups):
		if not o.clonal_group in clonal_group: continue
		samples = [s.rsplit('_',1)[0] for s in expra.samples[::2]]
		samples = [s for s in samples if any(s.startswith(clonal_group_start) or s.startswith('pool.'+clonal_group_start) for clonal_group_start in clonal_group.split('\t'))]
		xarr_n = []
		yarr_mono = []
		xarr_ctrl_n = []
		yarr_ctrl = []
		xarr_n_line = []
		yarr_mono_line = []
		yarr_ctrl_line = []
		for n in range(o.start_n, o.end_n+1):
			y_last = []
			y_ctrl_last = []
			if n == 1 and o.nonrandom_n1:
				for s in samples:
					c57_f, cast_f, bi_f, ti_used = count_mono([s])
					yarr_mono.append(cast_f+c57_f)
					y_last.append(cast_f+c57_f)
					xarr_n.append(n)
			else:
				for i in range(o.random_dots):
					c57_f, cast_f, bi_f, ti_used = count_mono(random.sample(samples, n))
					yarr_mono.append(cast_f+c57_f)
					y_last.append(cast_f+c57_f)
					xarr_n.append(n)
			xarr_n_line.append(n)
			yarr_mono_line.append(numpy.mean(y_last))
		
			if n==1:
				yarr_ctrl_line.append(numpy.mean(y_last))
				print n, yarr_mono_line[-1]
			else:
				for i in range(o.random_dots):
					c57_f, cast_f, bi_f, ti_used = count_mono(random.sample(samples, n), True)
					yarr_ctrl.append(cast_f+c57_f)
					y_ctrl_last.append(cast_f+c57_f)
					xarr_ctrl_n.append(n)
				yarr_ctrl_line.append(numpy.mean(y_ctrl_last))
				print n, yarr_mono_line[-1], yarr_ctrl_line[-1]
		
		pylab.title(clonal_group)
		pylab.plot(xarr_ctrl_n, yarr_ctrl, 'go')
		pylab.plot(xarr_n_line, yarr_ctrl_line, 'g-')
		pylab.plot(xarr_n, yarr_mono, 'ko')
		pylab.plot(xarr_n_line, yarr_mono_line, 'k-')
		
		if o.allelehits_multicell:
			expra = dr_tools.loadexpr(o.allelehits_multicell, True)
			xarr_n = []
			yarr_mono = []
			for s, n in multi_cell_samples.items():
				if o.start_n <= n <= o.end_n:
					xarr_n.append(n)
					c57_f, cast_f, bi_f, ti_used = count_mono([s])
					yarr_mono.append(cast_f+c57_f)
		
			xarr_obs_line = sorted(set(xarr_n))
			yarr_obs_line = [numpy.mean([y for n,y in zip(xarr_n,yarr_mono) if n==x]) for x in xarr_obs_line]
			pylab.plot([xarr_n_line[0]]+xarr_obs_line, [yarr_mono_line[0]]+yarr_obs_line, 'r-')
		
		try:pylab.plot(xarr_n, yarr_mono, 'ro')
		except:
			print xarr_n, yarr_mono
		
		pylab.savefig(o.figure)

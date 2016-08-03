from __future__ import division
import argparse, dr_tools, numpy, pylab

def gene_i_by_listf(genelistf_arr, expr):
	allowedgenes = set()
	for genelistf in genelistf_arr:
		allowedgenes |= set(dr_tools.loadlist(genelistf))
	return set(i for i,sym in enumerate(expr['symbols']) if sym in allowedgenes)
	

def MAfraction(expra, exprt, sample, minrpkm, genelist=None):
	count_bi, count_mono = 0, 0
	for sym, c57, cast in zip(expra['symbols'], expra[sample+'_c57only'], expra[sample+'_castonly']):
		if (genelist is None or sym in genelist) and exprt[sample][exprt.symbol_to_index[sym]] >= minrpkm:
			if c57 and cast: count_bi += 1
			elif c57 or cast: count_mono += 1
	return count_mono/(count_bi + count_mono)

if '__main__' == __name__:
	opts = argparse.ArgumentParser()
	opts.add_argument('inf')
	opts.add_argument('rpkmf_total')
	opts.add_argument('min_rpkm', type=float)
	opts.add_argument('--filter', nargs='+')
	opts.add_argument('-f', '--figf', default='plot_monoallelic_by_cell_minrpkm.pdf')
	opts.add_argument('-gi', '--genelistf_include', nargs='+')
	opts.add_argument('-ge', '--genelistf_exclude', nargs='+')
	opts.add_argument('--castfather', action='store_true')
	opts.add_argument('--infercross', action='store_true')
	opts.add_argument('--alg2', action='store_true')
	o = opts.parse_args()

	expr = dr_tools.loadexpr([o.inf], counts=True)
	exprt = dr_tools.loadexpr([o.rpkmf_total], counts=False)
	
	allowed_gene_i = gene_i_by_listf(o.genelistf_include, expr) if o.genelistf_include else None
	excluded_gene_i = gene_i_by_listf(o.genelistf_exclude, expr) if o.genelistf_exclude else None
	
	def rpkm(Ai, sample):
		Ti = exprt.ID_to_index[expr['IDs'][Ai]]
		return exprt[sample][Ti]
	
	for p in dr_tools.splitlines(o.inf):
		if p[0] == '#samples': samples = p[1:]; break
	
	fractions = [] # maternal only + paternal only
	mfractions = [] # maternal only
	fractions_all3 = [] # maternal+parternal+biallelic
	labels = []
	
	for s1, s2 in zip(samples[::2], samples[1::2]):
		if o.filter is not None and not any(part in s1.rsplit('_',1)[0] for part in o.filter): continue
		if s1.rsplit('_',1)[0] != s2.rsplit('_',1)[0] or not 'c57' in s1 or 'c57' in s2:
			print 'Error in pair:', s1, s2
			continue # check for errors in input file format
		s = s1.rsplit('_',1)[0]
		if o.alg2:
			try: f = MAfraction(expr, exprt, s, o.min_rpkm, allowed_gene_i)
			except:
				print s
				continue
			mfractions.append(0)
			fractions.append(f)
			fractions_all3.append(1)
			labels.append(s)
			print s, f
		else:
			Z = zip(expr[s1], expr[s2])
			Z = [E for i,E in enumerate(Z) if (allowed_gene_i is None or i in allowed_gene_i) and (excluded_gene_i is None or i not in excluded_gene_i) and rpkm(i, s) >= o.min_rpkm] # only those includes in the gene list
			count_only_1 = sum((e1 > 0) and (e2 == 0) for e1,e2 in Z)
			count_only_2 = sum((e1 == 0) and (e2 > 0) for e1,e2 in Z)
			count_both = sum(e1 > 0 and e2 > 0 for e1,e2 in Z)
			try:
				if o.castfather or (o.infercross and 'BxC' in s1):
					mfractions.append((count_only_1)/(count_both+count_only_1+count_only_2))
				else:
					mfractions.append((count_only_2)/(count_both+count_only_1+count_only_2))
				fractions.append((count_only_1+count_only_2)/(count_both+count_only_1+count_only_2))
				fractions_all3.append((count_only_1+count_only_2+count_both)/(count_both+count_only_1+count_only_2))
				labels.append(s1.rsplit('_',1)[0])
				print s1.rsplit('_',1)[0], fractions[-1], mfractions[-1], fractions[-1]-mfractions[-1]
			except:
				print s1, s2, len(allowed_gene_i), sum(expr[s1]), sum(expr[s2])
				continue
	pylab.bar(range(len(fractions)), fractions_all3, facecolor='#009933', linewidth=0, width=0.95)
	pylab.bar(range(len(fractions)), fractions, facecolor='b', linewidth=0, width=0.95)
	pylab.bar(range(len(fractions)), mfractions, facecolor='#cc0000', linewidth=0, width=0.95)
	pylab.xticks([x+0.95/2 for x in range(len(fractions))], labels, rotation=90, horizontalalignment='center', fontsize=3)
	pylab.ylim(0,1)
	pylab.subplots_adjust(bottom=0.3)
	pylab.savefig(o.figf)
	#print fractions

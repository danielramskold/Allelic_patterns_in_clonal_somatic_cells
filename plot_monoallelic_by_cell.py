from __future__ import division
import argparse, dr_tools, numpy, pylab

def gene_i_by_listf(genelistf_arr, expr):
	allowedgenes = set()
	for genelistf in genelistf_arr:
		allowedgenes |= set(dr_tools.loadlist(genelistf))
	return set(i for i,sym in enumerate(expr['symbols']) if sym in allowedgenes)
	

if '__main__' == __name__:
	opts = argparse.ArgumentParser()
	opts.add_argument('inf')
	opts.add_argument('--filter', nargs='+')
	opts.add_argument('-f', '--figf', default='plot_monoallelic_by_cell.pdf')
	opts.add_argument('-gi', '--genelistf_include', nargs='+')
	opts.add_argument('-ge', '--genelistf_exclude', nargs='+')
	o = opts.parse_args()

	expr = dr_tools.loadexpr([o.inf], counts=True)
	#samples = sorted([e for e in expr if e not in ('IDs', 'symbols')])
	
	allowed_gene_i = gene_i_by_listf(o.genelistf_include, expr) if o.genelistf_include else None
	excluded_gene_i = gene_i_by_listf(o.genelistf_exclude, expr) if o.genelistf_exclude else None
	
	
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
		Z = zip(expr[s1], expr[s2])
		if o.genelistf_include or o.genelistf_exclude:
			Z = [E for i,E in enumerate(Z) if (allowed_gene_i is None or i in allowed_gene_i) and (excluded_gene_i is None or i not in excluded_gene_i)] # only those includes in the gene list
		count_only_1 = sum((e1 > 0) and (e2 == 0) for e1,e2 in Z)
		count_only_2 = sum((e1 == 0) and (e2 > 0) for e1,e2 in Z)
		count_both = sum(e1 > 0 and e2 > 0 for e1,e2 in Z)
		try:
			fractions.append((count_only_1+count_only_2)/(count_both+count_only_1+count_only_2))
			mfractions.append((count_only_2)/(count_both+count_only_1+count_only_2))
			fractions_all3.append((count_only_1+count_only_2+count_both)/(count_both+count_only_1+count_only_2))
			labels.append(s1.rsplit('_',1)[0])
		except:
			print s1, s2#, len(allowed_gene_i), sum(expr[s1]), sum(expr[s2])
			continue
	pylab.bar(range(len(fractions)), fractions_all3, facecolor='#009933', linewidth=0, width=0.95)
	pylab.bar(range(len(fractions)), fractions, facecolor='b', linewidth=0, width=0.95)
	pylab.bar(range(len(fractions)), mfractions, facecolor='#cc0000', linewidth=0, width=0.95)
	pylab.xticks([x+0.95/2 for x in range(len(fractions))], labels, rotation=90, horizontalalignment='center', fontsize=3)
	pylab.ylim(0,1)
	pylab.subplots_adjust(bottom=0.3)
	pylab.savefig(o.figf)

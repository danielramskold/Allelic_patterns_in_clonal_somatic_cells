import dr_tools, argparse, numpy, math, random
from itertools import chain

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('infile')
	parser.add_argument('outfile')
	parser.add_argument('-m', '--maxgenes', type=int)
	parser.add_argument('-S', '--maxgeneselection', choices=['max', 'mean', 'random'], default='max')
	parser.add_argument('-t', '--transform', choices=['none', 'log10+0.3'], default='none')
	parser.add_argument('-c', '--centering', choices=['none', 'mean'], default='none')
	parser.add_argument('-s', '--samplelist')
	parser.add_argument('-e', '--excludesample', nargs='+')
	o = parser.parse_args()

	# load input
	expr = dr_tools.loadexpr(o.infile)
	
	# select samples
	if o.samplelist is not None:
		samples = dr_tools.loadlist(o.samplelist)
	else:
		samples = expr.samples
	if o.excludesample:
		samples = [s for s in samples if s not in o.excludesample]

	# select genes
	genes_i = range(len(expr['symbols']))
	if o.maxgenes is not None:
		select_fn = {'max':max, 'mean':numpy.mean, 'random': (lambda v: random.random())}[o.maxgeneselection]
		sort_list = [(select_fn([expr[s][i] for s in samples]), i) for i in genes_i]
		sort_list.sort(reverse=True)
		genes_i = [i for sort_val, i in sort_list[:o.maxgenes]]

	# sort genes alphabetically
	sort_list = [(expr['symbols'][i], i) for i in genes_i]
	sort_list.sort()
	genes_i = [i for sort_val, i in sort_list]

	# tranform the data
	if o.transform != 'none' or o.centering != 'none':
		rpkm_transform = {'none':float, 'log10+0.3': (lambda v: math.log10(v+0.3))}[o.transform]
		center_fn = {'mean':numpy.mean, 'none':(lambda v: 0)}[o.centering]
		for i in genes_i:
			if o.transform != 'none':
				for s in samples:
					expr[s][i] = rpkm_transform(expr[s][i])
			if o.centering != 'none':
				
				mid = center_fn([expr[s][i] for s in samples])
				for s in samples:
					expr[s][i] -= mid

	# write output
	with open(o.outfile, 'w') as outfh:
		print >>outfh, ','.join([expr['symbols'][i] for i in genes_i])
		for s in samples:
			print >>outfh, ','.join([str(expr[s][i]) for i in genes_i])
		

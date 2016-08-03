from scipy import stats
import pylab, argparse, dr_tools, itertools, random, itertools

def maxpairs(samples, nmax):
	pairs = list(itertools.combinations(samples, 2))
	if len(pairs) > nmax:
		return random.sample(pairs, nmax)
	else:
		return pairs

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-r', '--rpkmfile', nargs='+', required=True)
	parser.add_argument('-s', '--samplelist', nargs='+', required=True, action='append')
	parser.add_argument('-o', '--figure', default='correlation.pdf')
	parser.add_argument('-n', '--names', action='append')
	parser.add_argument('-m', '--maxpergroup', type=int, default=300000000)
	o = parser.parse_args()
	
	expr = dr_tools.loadexpr(o.rpkmfile)
	boxplot_values = []
	labels = []
	for samplelistgroup, name in itertools.izip_longest(o.samplelist, o.names):
		if samplelistgroup is None: raise Exception
		if name is None: label = ''
		else: label = name + '\n'
		rho_values = []
		samples_used = set()
		possible_pairs = 0
		for samplelistfile in samplelistgroup:
			samples = set(dr_tools.loadlist(samplelistfile))
			rho_values.extend([stats.spearmanr(expr[s1], expr[s2])[0] for s1, s2 in maxpairs(samples, o.maxpergroup)])
			samples_used.update(samples)
			possible_pairs += len(samples) * (len(samples)-1) // 2
		boxplot_values.append(rho_values)
		print name, 'samples=%d correlations=%d genes=%d possible_pairs=%d'%(len(samples_used),len(rho_values), len(expr[s1]), possible_pairs)
		label += 'n=%d'%len(rho_values)
		labels.append(label)
	pylab.ylim(0,1)
	dr_tools.violin_plot(pylab.axes(), boxplot_values, range(len(boxplot_values)), bp=True)
	pylab.xticks(range(len(boxplot_values)), labels, rotation=90)
	pylab.subplots_adjust(bottom=0.3)
	pylab.savefig(o.figure)

# example command:
# python correlation_boxplot.py -r ../bjorn_reinius_complete_set/rpkmforgenes/star_merged_cast1_mm9/ensembl/rpkms_counts_rmnameoverlap.txt -s sample_lists/excltechrepl/* -s sample_lists/splitcell_big/* -s sample_lists/splitcell_small/* -s sample_lists/samples_big.txt -s sample_lists/samples_small.txt -s sample_lists/fibroblast_samples_ERCC.txt -s sample_lists/livercellsamples.txt -s sample_lists/brain_samples.txt -n tech_repl -n big_splits -n small_splits -n big_fibroblasts -n small_fibroblasts -n fibroblasts_other  -n liver -n brain

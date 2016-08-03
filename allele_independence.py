from __future__ import division
import argparse, dr_tools, numpy, pylab, random
from math import sqrt, isnan
from collections import defaultdict, namedtuple

def pairs(expr, stages, exclude):
	if stages is None:
		samples = expr.samples
	else:
		samples = [s for s in expr.samples if any(part in s for part in stages) and not any(part in s for part in exclude)]
	sample_pairs = zip(samples[::2], samples[1::2])
	if not all('c57only' in s_pat and 'castonly' in s_mat for s_pat, s_mat in sample_pairs):
		raise Exception
	return sample_pairs

'''
activation coefficient r, likelihood of activating a gene
biallelic genes b: r*r
paternal only genes p: r - r*r
maternal only genes m: r - r*r
non-expressed genes n: (1-r)**2
r=0.5-sqrt(0.25-(m+p)/2) or r=0.5+sqrt(0.25-(m+p)/2)
'''

class Bin:
#	def __init__(self, num_mono, num_samples):
#		self.frac_mono = num_mono/num_samples
#		self.N = num_samples
#		if self.frac_mono >= 0.5: # 0.5 = (0.5-0.5**2)*2
#			r = 0.5+sqrt(0.25-self.frac_mono/2)
#		else:
#			r = 0.5-sqrt(0.25-self.frac_mono/2)
#		self.exp_frac_bi = r**2 # expected
#		self.obs_frac_bi = [] # observed
	
	def __init__(self, num_silent, num_samples, simulate):
		self.frac_silent = num_silent/num_samples
		self.N = num_samples
		r = 1-sqrt(self.frac_silent)
		self.exp_frac_bi = [] if simulate else [r**2] # expected
		self.obs_frac_bi = [] # observed
	
	def add(self, num_bi):
		self.obs_frac_bi.append(num_bi/self.N)
	def addsim(self, num_bi):
		self.exp_frac_bi.append(num_bi/self.N)
	
	def hasNaN(self):
		return isnan(numpy.mean(self.exp_frac_bi)) or isnan(numpy.mean(self.obs_frac_bi))

def sem(values):
	if len(values) < o.minN: return float('nan')
	return numpy.std(values)/sqrt(len(values))

def violin_plot(ax,data,pos, bp=False, leftside=True, rightside=True, widthf=None, color='y'):
    '''
    create violin plots on an axis
    run with e.g violin_plot(pylab.axes(), [[3,4,5],[7,8]], [0, 1])
    '''
    # from http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
    from matplotlib.pyplot import figure, show
    from scipy.stats import gaussian_kde
    from numpy.random import normal
    from numpy import arange
    
    dist = max(pos)-min(pos)
    w = min(0.15*max(dist,1.0),0.5)
    for d,p in zip(data,pos):
        k = gaussian_kde(d) #calculates the kernel density
        m = k.dataset.min() #lower bound of violin
        M = k.dataset.max() #upper bound of violin
        x = arange(m,M,(M-m)/100.) # support for violin
        v = k.evaluate(x) #violin profile (density curve)
        v = v/v.max()*w #scaling the violin to the available space
        if widthf is not None: v = v*widthf
        if rightside: ax.fill_betweenx(x,p,v+p,facecolor=color,alpha=0.3)
        if leftside: ax.fill_betweenx(x,p,-v+p,facecolor=color,alpha=0.3)
    if bp:
        ax.boxplot(data,notch=1,positions=pos,vert=1)

'''
danielr@rna ~/casthybrid/one_chr_reads $ python allele_independence.py -i ~/casthybrid/snp_positions/allelecounts_from_pileup/v17S15_genomic_refseq_autosomes.txt --stages 8cell --exclude 8cell_8-3
to allele_independence_8cell.pdf
danielr@rna ~/casthybrid/one_chr_reads $ python allele_independence.py -i ~/casthybrid/snp_positions/allelecounts_from_pileup/v17S15_genomic_refseq_autosomes.txt --stages 4cell
to allele_independence_4cell.pdf
danielr@rna ~/casthybrid/one_chr_reads $ python allele_independence.py -i ~/casthybrid/snp_positions/allelecounts_from_pileup/v17S15_genomic_refseq_autosomes.txt --stages 16cell --exclude 16cell_4-2
to allele_independence_16cell.pdf
danielr@rna ~/casthybrid/one_chr_reads $ python allele_independence.py -i ~/casthybrid/snp_positions/allelecounts_from_pileup/v17S15_genomic_refseq_autosomes.txt --stages blast
to allele_independence_blastocyst.pdf
'''

if '__main__' == __name__:
	opts = argparse.ArgumentParser()
	opts.add_argument('-i', '--inf', nargs='+', required=True)
	opts.add_argument('--stages', nargs='+', help='when there is not a genomewide maternal bias')
	opts.add_argument('--exclude', nargs='+', help='remove from stages', default=[])
	opts.add_argument('--sim', action='store_true')
	opts.add_argument('-o', '--figure', default='allele_independence.pdf')
	opts.add_argument('--plotstyle', default=['mean_graph'], choices=['mean_graph', 'boxplot', 'mean_sem', 'violin', 'std', 'sayN', 'sayY'], nargs='+')
	opts.add_argument('--minN', default=0, type=int)
	o = opts.parse_args()
	
	expra = dr_tools.loadexpr(o.inf, counts=True)
	sample_pairs = pairs(expra, o.stages, o.exclude)
	
	bins = [Bin(num_cells, len(sample_pairs), o.sim) for num_cells in range(len(sample_pairs)+1)]
	
	for gene_i in range(len(expra['symbols'])):
		num_mono = sum((expra[s_pat][gene_i]>0)^(expra[s_mat][gene_i]>0) for s_pat,s_mat in sample_pairs)
		num_bi = sum((expra[s_pat][gene_i]>0)and(expra[s_mat][gene_i]>0) for s_pat,s_mat in sample_pairs)
		num_silent = sum((expra[s_pat][gene_i]==0)and(expra[s_mat][gene_i]==0) for s_pat,s_mat in sample_pairs)
		
		bins[num_silent].add(num_bi)
	
	if o.sim:
		while any(len(b.exp_frac_bi) < 10000 for b in bins):
			r = random.random()**2
			sim_states = [(random.random() < r, random.random() < r) for p in sample_pairs]
			num_bi = sum(p and m for p,m in sim_states)
			num_silent = sum(not (p or m) for p,m in sim_states)
			bins[num_silent].addsim(num_bi)
	
	if 'sayN' in o.plotstyle:
		Ns = [len(b.obs_frac_bi) for b in bins[:-1]]
		print 'min N =', min(Ns), 'max =', max(Ns)
	if 'violin' in o.plotstyle:
		violin_plot(pylab.axes(), [b.obs_frac_bi for b in bins[:-1]], [b.frac_silent for b in bins[:-1]], leftside=False, color='b', widthf=0.1)
		violin_plot(pylab.axes(), [random.sample(b.exp_frac_bi, len(o.obs_frac_bi)) for b in bins[:-1]], [b.frac_silent for b in bins[:-1]], rightside=False, color='r', widthf=0.1)
	if 'boxplot' in o.plotstyle:
		pylab.plot([b.frac_silent for b in bins], [numpy.mean(b.exp_frac_bi) for b in bins], color='r', label='expected')
		pylab.boxplot([b.obs_frac_bi for b in bins], positions=[b.frac_silent for b in bins], widths=0.5/len(sample_pairs))
	if 'mean_graph' in o.plotstyle:
		pylab.plot([b.frac_silent for b in bins], [numpy.mean(b.exp_frac_bi) for b in bins], color='r', label='expected')
		pylab.plot([b.frac_silent for b in bins], [numpy.mean(b.obs_frac_bi) for b in bins], color='b', label='observed')
	if 'mean_sem' in o.plotstyle:
		pylab.plot([b.frac_silent for b in bins], [numpy.mean(b.exp_frac_bi) for b in bins], color='r', label='expected')
		pylab.errorbar([b.frac_silent for b in bins], [numpy.mean(b.obs_frac_bi) for b in bins], [sem(b.obs_frac_bi) for b in bins], label='observed')
	if 'std' in o.plotstyle:
		pylab.plot([b.frac_silent for b in bins], [numpy.mean(b.exp_frac_bi) for b in bins], color='r', label='expected')
		pylab.errorbar([b.frac_silent for b in bins], [numpy.mean(b.obs_frac_bi) for b in bins], [numpy.std(b.obs_frac_bi) for b in bins], label='observed')
	if 'sayY' in o.plotstyle:
		print [numpy.mean(b.obs_frac_bi) for b in bins]
		print [numpy.mean(b.exp_frac_bi) for b in bins]
		from scipy import stats
		print 'paired t test', stats.ttest_rel([numpy.mean(b.obs_frac_bi) for b in bins if not b.hasNaN()], [numpy.mean(b.exp_frac_bi) for b in bins if not b.hasNaN()])
	pylab.title(' '.join(o.stages))
	pylab.xlabel('silent fraction of cells')
	pylab.ylabel('biallelic fraction of cells')
	pylab.legend()
	pylab.savefig(o.figure)

#http://nxn.se/post/46440196846/making-nicer-looking-pie-charts-with-matplotlib
import numpy as np
import matplotlib.pyplot as plt
import dr_tools, argparse

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('autosomal_mono', type=float)
	parser.add_argument('shared_mono_withimpr', type=float)
	parser.add_argument('shared_mono_exclimpr', type=float)
	parser.add_argument('shared_mono_allelerand', type=float)
	parser.add_argument('chrX_mono', type=float)
	parser.add_argument('chrX_mono_male', type=float)
	parser.add_argument('output_picture')
	o = parser.parse_args()
	
	Xpart = 1-0.9320448877805486
	
	imprinted = o.shared_mono_withimpr - o.shared_mono_exclimpr
	clonal = o.shared_mono_exclimpr - o.shared_mono_allelerand
	dynamic = o.autosomal_mono - o.shared_mono_withimpr
	autosomal_mono_factor = (1-Xpart)*o.autosomal_mono/(imprinted+clonal+dynamic)
	X_bi = (o.chrX_mono_male - o.chrX_mono)/o.chrX_mono_male
	
	fig = plt.figure(figsize=[10, 10])
	ax = fig.add_subplot(111)
	cmap = plt.cm.prism
	slices = [Xpart*X_bi, Xpart*(1-X_bi), autosomal_mono_factor*imprinted, autosomal_mono_factor*clonal, autosomal_mono_factor*dynamic, (1-Xpart)*(1-o.autosomal_mono)]
	colors = cmap(np.linspace(0., 1., len(slices)))
	colors = ['#cc55ff', '#22bb77', '#ee0033', '#00ee33', '#2233aa', '#6677aa']
	labels = ['X escapee', 'X mono', 'imprinted mono', 'clonal RME', 'dynamic RME', 'autosomal bi']
	ax.pie(slices, colors=colors, labels=labels, labeldistance=1.05)
	fig.savefig(o.output_picture)

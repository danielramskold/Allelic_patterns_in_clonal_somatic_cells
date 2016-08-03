import argparse, dr_tools
from collections import defaultdict
from itertools import chain

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('tcr_summary_file', metavar='tcr_summary.txt')
	parser.add_argument('output_prefix')
	parser.add_argument('--numbering_start', type=int, default=1)
	o = parser.parse_args()
	num = o.numbering_start

	clones_by_TCR_seq = defaultdict(lambda: defaultdict(list))
	for li, p in enumerate(dr_tools.splitlines(o.tcr_summary_file)):
		if li == 0:
			VDJ_i = []
			CDR3_i = []
			for column_header in ('CDR3 amino acid sequence',):
				CDR3_i.append(p.index(column_header)+1)
				CDR3_i.append(p.index(column_header)+1+len(p))
			for column_header in ('V alleles', 'J alleles'):
				VDJ_i.append(p.index(column_header)+1)
			for column_header in ('V alleles', 'D alleles', 'J alleles'):
				VDJ_i.append(p.index(column_header)+1+len(p))
		else:
			sample_name = p[0]
			source_person = sample_name.split('_')[0]
			try: TCR_features = tuple(p[i] for i in chain(CDR3_i, VDJ_i))
			except IndexError:
				# discard, it is missing info
				continue
			if '' in TCR_features: continue # discard if missing info
			clones_by_TCR_seq[source_person][TCR_features].append(sample_name)

	for source_person in clones_by_TCR_seq:
		for TCR_features, samples in clones_by_TCR_seq[source_person].items():
			if len(samples) >= 2:
				dr_tools.printlist('%s%d_%s.txt'%(o.output_prefix, num, source_person), samples)
			num += 1

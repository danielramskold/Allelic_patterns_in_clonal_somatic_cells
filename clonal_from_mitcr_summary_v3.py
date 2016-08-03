import argparse, dr_tools
from collections import defaultdict
from itertools import chain

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--tcr_summary_file', metavar='tcr_summary.txt', nargs='+', required=True)
	parser.add_argument('-o', '--output_prefix', required=True)
	parser.add_argument('-n', '--numbering_start', type=int, default=1)
	parser.add_argument('-p', '--oneperson', nargs='+', help='e.g. male_P1299_YFV2001 P1299 YFV2001', action='append', default=[])
	o = parser.parse_args()
	num = o.numbering_start

	source_person_mapping = dict()
	for names in o.oneperson:
		for name in names:
			source_person_mapping[name] = names[0]

	clones_by_TCR_seq = defaultdict(lambda: defaultdict(list))
	for tcr_summary_file in o.tcr_summary_file:
		for li, p in enumerate(dr_tools.splitlines(tcr_summary_file)):
			if li == 0:
				VDJ_i = []
				CDR3_i = []
				for column_header in ('CDR3 amino acid sequence',):
					CDR3_i.append(p.index(column_header)+1)
					CDR3_i.append(p.index(column_header)+1+len(p))
				for column_header in ('V segments', 'J segments'):
					VDJ_i.append(p.index(column_header)+1)
				for column_header in ('V segments', 'D segments', 'J segments', 'VD insertions', 'DJ insertions'):
					VDJ_i.append(p.index(column_header)+1+len(p))
			else:
				sample_name = p[0]
				if not p[1]: p = [p[0]]+p[2:]
				source_person = source_person_mapping.get(sample_name.split('_')[0], sample_name.split('_')[0])
				try: TCR_features = tuple(p[i] for i in chain(CDR3_i, VDJ_i))
				except IndexError:
					# discard, it is missing info
					continue
				if '' in TCR_features: continue # discard if missing info
				if any(cdr3.isdigit() for cdr3 in TCR_features[:2]):
					raise Exception, 'should be AA, but number in ' +repr(TCR_features[:2])
				clones_by_TCR_seq[source_person][TCR_features].append(sample_name)

	for source_person in clones_by_TCR_seq:
		for TCR_features, samples in clones_by_TCR_seq[source_person].items():
			if len(samples) >= 2:
				dr_tools.printlist('%s%d_%s.txt'%(o.output_prefix, num, source_person), samples)
				num += 1

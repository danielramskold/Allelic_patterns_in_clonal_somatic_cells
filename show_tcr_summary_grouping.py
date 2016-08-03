import argparse, dr_tools

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('tcr_summary')
	parser.add_argument('clone_list', nargs='+')
	o = parser.parse_args()

	summary_lines = dict()
	with open(o.tcr_summary, 'rU') as infh:
		for line in infh:
			sample = line.split('\t', 1)[0]
			summary_lines[sample] = line

	for clonelistfile in o.clone_list:
		samples = dr_tools.loadlist(clonelistfile)
		print '\n#', clonelistfile
		for sample in samples:
			print summary_lines[sample].rstrip('\r\n')
import argparse

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--files', nargs='+', required=True)
	parser.add_argument('-s', '--suffix', nargs='+', default=['_expression.txt', '.fastq.gz'])
	parser.add_argument('-l', '--samplelist', required=True)
	parser.add_argument('-L', '--outputlisting', choices=['matching', 'extra', 'missing', 'target', 'existing'], default='missing')
	parser.add_argument('-t', '--translation')
	o = parser.parse_args()
	
	
	with open(o.samplelist, 'rU') as infh:
		target_samples = [l.rstrip('\r\n') for l in infh.readlines()]
	
	if o.translation:
		translation = dict()
		with open(o.translation, 'rU') as infh:
			for l in infh:
				p = l.strip('\r\n').split('\t')
				translation[p[1]] = p[0]
		target_samples = [translation.get(s, s) for s in target_samples]
	
	target_files = set(sample+suffix for sample in target_samples for suffix in o.suffix)
	found_files = dict((path.split('/')[-1], path) for path in o.files)
	
	list_out = []
	if o.outputlisting == 'matching':
		m = set(found_files) & target_files
		list_out.extend([found_files[f] for f in m])
	elif o.outputlisting == 'extra':
		m = set(found_files) - target_files
		list_out.extend([found_files[f] for f in m])
	elif o.outputlisting == 'missing':
		m = target_files - set(found_files)
		list_out.extend(list(m))
	elif o.outputlisting == 'target':
		list_out.extend(list(target_files))
	elif o.outputlisting == 'existing':
		list_out.extend(list(found_files.keys()))
	
	try:
		for f in list_out: print f
	except IOError: pass

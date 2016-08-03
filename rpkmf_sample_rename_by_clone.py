import dr_tools, argparse

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--rpkmf_in', required=True)
	parser.add_argument('-o', '--rpkmf_out', required=True)
	parser.add_argument('-s', '--sample_lists', nargs='+', required=True)
	o = parser.parse_args()
	
	with open(o.rpkmf_out, 'w') as outfh:
		with open(o.rpkmf_in, 'r') as infh:
			for li, line in enumerate(infh):
				if li == 0:
					p = line.rstrip('\r\n').split('\t')
					sample_to_clone = dict((sample, filename) for filename in o.sample_lists for sample in dr_tools.loadlist(filename))
					for i, name in enumerate(p):
						if i==0: continue
						for suffix in ('', '_c57only', '_castonly'):
							if name.endswith(suffix) and name[:-len(suffix)] in sample_to_clone:
								clone_name = sample_to_clone[name[:-len(suffix)] ].split('/')[-1].split('.txt')[0]
								p[i] = clone_name + '-' + name
					print >>outfh, dr_tools.join(p)
				else:
					outfh.write(line)

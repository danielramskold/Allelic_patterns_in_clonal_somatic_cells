import argparse,gzip

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('fastqgz', nargs='+')
	o = parser.parse_args()
	
	for fqz in o.fastqgz:
		infh = gzip.open(fqz, 'r')
		try:
			line = next(infh)
			line = next(infh)
			print len(line.rstrip()), fqz
			l = len(line.rstrip())
			line = next(infh)
			line = next(infh)
			line = next(infh)
			line = next(infh)
			if not l == len(line.rstrip()):
				print 'Uneven length:', fqz
		except:
			print fqz
			raise
		infh.close()

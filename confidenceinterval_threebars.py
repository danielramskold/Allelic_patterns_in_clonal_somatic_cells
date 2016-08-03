import argparse, numpy, math

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('infile', metavar='file_from_pool_and_plot4bars')
	o = parser.parse_args()
	
	shared=[]
	randomshared=[]
	with open(o.infile, 'rU') as infh:
		for line in infh:
			if line.strip().endswith('X shared'):
				pass
			elif line.strip().endswith('random shared'):
				randomshared.append(float(line.split()[0]))
			elif line.strip().endswith('shared'):
				shared.append(float(line.split()[0]))
	
	X = [s-r for r,s in zip(randomshared, shared)]
	mid = numpy.mean(X)
	sem = numpy.std(X)/math.sqrt(len(X))
	print mid+sem*1.96, mid-sem*1.96

import argparse, dr_tools

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('infile')
	parser.add_argument('samplepart')
	o = parser.parse_args()
	
	samples = dr_tools.loadlist(o.infile)
	samples = [s for s in samples if o.samplepart in s]
	prefix = samples[0].rsplit(o.samplepart)[0]+o.samplepart
	dr_tools.printlist(prefix+'.txt', samples)

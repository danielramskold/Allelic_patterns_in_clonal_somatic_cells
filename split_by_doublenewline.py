import argparse, os

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('infile')
	parser.add_argument('output_prefix')
	parser.add_argument('numbering_start', type=int)
	o = parser.parse_args()
	num = o.numbering_start
	outfh = open(o.output_prefix+str(num)+'.txt', 'w')
	written = False
	with open(o.infile, 'rU') as infh:
		for line in infh:
			if line.strip():
				outfh.write(line)
				written = True
			elif written:
				outfh.close()
				num+=1
				outfh = open(o.output_prefix+str(num)+'.txt', 'w')
				written = False
			else: pass
	outfh.close()
	if not written:
		os.remove(o.output_prefix+str(num)+'.txt')

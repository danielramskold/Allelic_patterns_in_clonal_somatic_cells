import argparse

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('infile')
	parser.add_argument('--predicted_too', action='store_true')
	o = parser.parse_args()

	with open(o.infile, 'rU') as infh:
		for line in infh:
			p = line.split()
			if len(p) < 4: continue
			if not p[-1] in ('Maternal', 'Paternal'): continue
			if not (p[-2] == 'Imprinted' or (o.predicted_too and p[-2] == 'Predicted')): continue
			print p[0]

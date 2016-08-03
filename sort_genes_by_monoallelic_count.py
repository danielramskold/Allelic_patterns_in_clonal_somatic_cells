import argparse, dr_tools

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('rpkms', required=True)
	parser.add_argument('allelehits', required=True)
	parser.add_argument('--minrpkm', type=float, default=20)
	o = parser.parse_args()
	
	exprt = dr_tools.loadexpr(o.rpkms, False)
	expra = dr_tools.loadexpr(o.allelehits, True)
	
	samples = set(exprt.samples) & set(s.rsplit('_',1)[0] for s in expra.samples[::2])
	
	count_per_gene = dict()
	
	assert expra['symbols'] == exprt['symbols']
	
	for ti, sym in 
	# not done

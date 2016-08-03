import argparse, dr_tools, os

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('rpkmf_alleles', nargs='?')
	parser.add_argument('rpkmf_total')
	o = parser.parse_args()
	
	exprt = dr_tools.loadexpr(o.rpkmf_total, counts=False)
	counts = dr_tools.loadexpr(o.rpkmf_total, counts=True)
	
	if o.rpkmf_alleles:
		expra = dr_tools.loadexpr(o.rpkmf_alleles, counts=True)
	
	
		AiD = dict((ti, expra.ID_to_index[ID]) for ti, ID in enumerate(exprt['IDs']) if ID in expra.ID_to_index)
	
		for s in exprt.samples:
			if s+'_castonly' not in expra.samples: continue
			with open(s + '_expression.txt', 'w') as outfh:
				print >>outfh, dr_tools.join('#Gene_symbol', 'Refseq_IDs', 'RPKM', 'reads', 'CAST_hits', 'C57_hits')
				for ti in range(len(exprt['IDs'])):
					if ti in AiD:
						ai = AiD[ti]
						cast = int(expra[s+'_castonly'][ai])
						c57 = int(expra[s+'_c57only'][ai])
					else:
						cast = 0
						c57 = 0
					rpkm = exprt[s][ti]
					reads = int(round(counts[s][ti]))
					symbol = exprt['symbols'][ti].replace('+','|')
					ID = exprt['IDs'][ti].replace('+','|')
					print >>outfh, dr_tools.join(symbol, ID, rpkm, reads, cast, c57)
	else:
		for s in exprt.samples:
			with open(s + '_expression.txt', 'w') as outfh:
				print >>outfh, dr_tools.join('#Gene_symbol', 'Refseq_IDs', 'RPKM', 'reads')
				for ti in range(len(exprt['IDs'])):
					rpkm = exprt[s][ti]
					reads = int(round(counts[s][ti]))
					symbol = exprt['symbols'][ti].replace('+','|')
					ID = exprt['IDs'][ti].replace('+','|')
					print >>outfh, dr_tools.join(symbol, ID, rpkm, reads)

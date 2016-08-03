import argparse, sys, time, dr_tools

def rename(name):
	return name

def num(name):
	return 0

def swap_order(twovalues):
	V = twovalues.split()
	return V[1] + '\t' + V[0]

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('inf')
	parser.add_argument('outf')
	parser.add_argument('-r', '--rpkmf_getID')
	parser.add_argument('--noNA', action='store_true')
	parser.add_argument('-s', '--samplenames')
	parser.add_argument('--rpkmf_genes', action='store_true')
	parser.add_argument('--exclude_genes')
	o = parser.parse_args()
	
	exclude_genes = set(dr_tools.loadlist(o.exclude_genes)) if o.exclude_genes else set()
	
	# load columns
	column_to_name = dict()
	symbols = []
	IDs = []
	sample_values = {}
	with open(o.inf, 'r') as infh:
		for li, line in enumerate(infh):
			p = line.rstrip('\r\n').split('\t')
			if li == 0:
				for ci,name in enumerate(p):
					if name in ('name','transcript','nbSNPs','SNPlocations', 'c57:cast', 'chrom', 'pos', 'genexcells'):
						continue
					newname = rename(name)
					sample_values[newname] = []
					column_to_name[ci] = newname
			else:
				symbol = p[0]
				if symbol in exclude_genes:
					# pretend there wasn't any value in the input file
					continue
				symbols.append(symbol)
				IDs.append(p[1])
				for ci, name in column_to_name.items():
					try: sample_values[name].append(p[ci])
					except:
						print p, ci, len(p)
						raise
	
	# sort the columns
	sample_order = [name for num_out,name in sorted((num(name), name) for name in sample_values)]
	
	# add in removal of e.g. midblast_2-19,midblast_2-20,midblast_2-22
	if o.samplenames:
		with open(o.samplenames, 'r') as infh:
			requested_samples = set(line.split()[0] for line in infh)
		sample_order = [name for name in sample_order if name in requested_samples]
		if requested_samples - set(sample_order):
			print 'Missing:\n' + '\n'.join(list(requested_samples - set(sample_order)))
	
	# change ID column
	if o.rpkmf_getID:
		expr = dr_tools.loadexpr(o.rpkmf_getID)
		symbol_to_IDs = dict(zip(expr['symbols'],expr['IDs']))
		#IDs = [symbol_to_IDs.get(sym, prevID) for prevID, sym in zip(IDs, symbols)]
		IDs = [symbol_to_IDs.get(sym, 'NA') for prevID, sym in zip(IDs, symbols)]
		if o.rpkmf_genes:
			symbols_set = dict((s,i) for i,s in enumerate(symbols))
			new_sample_values = dict()
			for name in sample_order:
				new_sample_values[name] = []
				for i, symbol in enumerate(expr['symbols']):
					if symbol in symbols_set:
						new_sample_values[name].append(sample_values[name][symbols_set[symbol]])
					else:
						new_sample_values[name].append('0 0')
			sample_values = new_sample_values
			symbols = expr['symbols']
			IDs = expr['IDs']
	elif o.rpkmf_genes: raise Exception
	
	# write to file
	with open(o.outf, 'w') as outfh:
		print >>outfh, dr_tools.join('#samples', ['%s_c57only\t%s_castonly'%(s,s) for s in sample_order])
		print >>outfh,  dr_tools.join('#allmappedreads', ['0\t0' for s in sample_order])
		print >>outfh,  dr_tools.join('#normalizationreads', ['0\t0' for s in sample_order])
		print >>outfh,  dr_tools.join('#arguments', ' '.join(sys.argv), 'time: '+time.asctime())
		for i in range(len(symbols)):
			#if IDs[i] == '0 0':
			#	print symbols[i]
			#	IDs[i] = 'NA'
			if o.noNA and IDs[i] == 'NA': continue
			try: print >>outfh, dr_tools.join(symbols[i], IDs[i], ['0\t0' for name in sample_order], [swap_order(sample_values[name][i]) for name in sample_order])
			except:
				print symbols[i], sample_values[name][i]
				raise

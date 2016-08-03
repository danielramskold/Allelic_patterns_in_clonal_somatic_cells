import dr_tools, argparse, pylab, numpy, math

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-r', '--rpkmfile', required=True)
	parser.add_argument('-g', '--genelist', required=True, nargs='+')
	parser.add_argument('-o', '--figure', default='rpkm_distribution.pdf')
	o = parser.parse_args()
	
	expr = dr_tools.loadexpr(o.rpkmfile)
	
	genes = set.union(*(set(dr_tools.loadlist(filename)) for filename in o.genelist))
	
	rpkms_genelist = []
	for sym in genes:
		try: ti = expr.symbol_to_index[sym]
		except KeyError: continue
		meanrpkm = numpy.mean([expr[s][ti] for s in expr.samples])
		rpkms_genelist.append(math.log(max(2**-10, meanrpkm), 2))
	pylab.hist(rpkms_genelist, 10)
	pylab.xlabel('log2 mean rpkm (-10 = not expressed)')
	pylab.ylabel('# genes')
	pylab.savefig(o.figure)

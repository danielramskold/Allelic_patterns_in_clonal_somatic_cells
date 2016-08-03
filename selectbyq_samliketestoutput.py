import argparse, dr_tools

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('infile')
	parser.add_argument('genelist_out')
	parser.add_argument('-g', '--genelist_in')
	parser.add_argument('-q', '--maxq', default=0.05, type=float)
	parser.add_argument('--top', action='store_true')
	parser.add_argument('--bottom', action='store_true')
	parser.add_argument('--oneID', action='store_true')
	o = parser.parse_args()
	if not o.top and not o.bottom:
		o.top = True
		o.bottom = True

	at_list_top = True
	if o.genelist_in: allowedgenes = set(dr_tools.loadlist(o.genelist_in))
	last_q = 0
	genes_top = []
	genes_bottom = []
	for li, p in enumerate(dr_tools.splitlines(o.infile)):
		if li == 0:
			q_i = p.index('FDR')
			sym_i = 0
			if o.oneID: ID_i = p.index('IDs')
		else:
			qval = float(p[q_i])
			if at_list_top and qval < last_q: at_list_top = False
			sym = p[sym_i]
			sym_out = p[ID_i].split('+')[0] if o.oneID else sym
			last_q = qval
			if o.genelist_in and sym not in allowedgenes: continue
			if qval > o.maxq: continue
			if at_list_top:
				genes_top.append(sym_out)
			else:
				genes_bottom.append(sym_out)
	genes_out = []
	if o.top:
		genes_out.extend(genes_top)
	if o.bottom:
		if o.top:
			genes_out.extend(genes_bottom)
		else:
			genes_out.extend(reversed(genes_bottom))
	dr_tools.printlist(o.genelist_out, genes_out)

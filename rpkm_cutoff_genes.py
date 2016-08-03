from __future__ import division
import argparse, dr_tools, numpy, pylab, random

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('-r', '--rpkms', required=True)
	parser.add_argument('-m', '--minrpkm', default=20, type=float)
	parser.add_argument('-M', '--maxrpkm', type=float)
	parser.add_argument('-gi', '--allowedgenes')
	parser.add_argument('-ge', '--disallowedgenes', nargs='+')
	o = parser.parse_args()
	
	exprt = dr_tools.loadexpr(o.rpkms)
	allowedgenes = set(dr_tools.loadlist(o.allowedgenes)) if o.allowedgenes else None
	if o.disallowedgenes:
		disallowedgenes = set()
		for filename in o.disallowedgenes:
			disallowedgenes.update(set(dr_tools.loadlist(filename)))
	else:
		disallowedgenes = None
	
	samples = exprt.samples
	
	for ti, sym in enumerate(exprt['symbols']):
		meanexpr = numpy.mean([exprt[s][ti] for s in samples])
		if meanexpr < o.minrpkm: continue
		if o.maxrpkm is not None and meanexpr >= o.maxrpkm: continue
		if disallowedgenes and sym in disallowedgenes: continue
		if allowedgenes and sym not in allowedgenes: continue
		print sym

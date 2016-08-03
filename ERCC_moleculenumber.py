from __future__ import division
import argparse, pylab, dr_tools, math, numpy
from collections import defaultdict

def calc_ERCC_moleculenumber(tablefile, before_dilution_vol_ul):
	Mix1_i = 3
	conc_attomolul = 0
	attomol = 602214.12927
	for i, p in enumerate(dr_tools.splitlines(tablefile)):
		if i == 0:
			if not 'attomoles/ul' in p[Mix1_i]: raise Exception
			if not 'Mix 1' in p[Mix1_i]: raise Exception
		else:
			conc_attomolul += float(p[Mix1_i])
	return conc_attomolul * before_dilution_vol_ul * 602214.12927

if '__main__' == __name__:
	ERCCvol_ul = 0.1/40000
	ERCC_moleculenumber = calc_ERCC_moleculenumber('ERCC.txt', ERCCvol_ul)
	print ERCC_moleculenumber 

from __future__ import division
import argparse, dr_tools, numpy

if '__main__' == __name__:
	
	# load table
	linefeed = dr_tools.splitlines('chrX_clones_allelic_calls.txt')
	sample_labels = next(linefeed)[1:]
	character_matrixT = []
	for cells in linefeed:
		# values in cells are nd, XI, xa, bi, except first column which is gene symbol
		if any(c!='nd' for c in cells[1:]):
			Xi_expr = sum(c in ('XI','bi') for c in cells[1:])
			tot_expr = sum(c in ('XI','bi','xa') for c in cells[1:])
			if 1.0 > 1 - Xi_expr/tot_expr >= 0.9:
				character_matrixT.append(cells[1:])
	character_matrix = zip(*character_matrixT)
	for label, sample_values in zip(sample_labels, character_matrix):
		print label, sample_values.count('XI')+sample_values.count('bi')

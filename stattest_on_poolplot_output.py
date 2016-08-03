import argparse
from scipy import stats
from collections import defaultdict

if '__main__' == __name__:
	parser = argparse.ArgumentParser()
	parser.add_argument('input', type=argparse.FileType('r'))
	o = parser.parse_args()
	
	rand_c57 = defaultdict(list)
	rand_cast = defaultdict(list)
	real_c57 = dict()
	real_cast = dict()
	for line in o.input:
		p = line.rstrip().split(' ')
		if p[1] == '': pass
		elif p[1] == 'r':
			rand_c57[clonal_group].append(float(p[2]))
			rand_cast[clonal_group].append(float(p[3]))
		else:
			clonal_group = p[1]
			real_c57[clonal_group] = float(p[2])
			real_cast[clonal_group] = float(p[3])
	
	for clonal_group in real_c57:
		Z1, P1 = stats.ttest_1samp(rand_c57[clonal_group], real_c57[clonal_group])
		Z2, P2 = stats.ttest_1samp(rand_cast[clonal_group], real_cast[clonal_group])
		print clonal_group, P1, Z1<0, P2, Z2<0

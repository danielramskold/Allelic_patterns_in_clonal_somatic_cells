from __future__ import division
import argparse, dr_tools, numpy, math

if '__main__' == __name__:
	opts = argparse.ArgumentParser()
	opts.add_argument('inf')
	opts.add_argument('--filter', nargs='+')
	o = opts.parse_args()

	expr = dr_tools.loadexpr([o.inf], counts=True)
	samples = sorted([e for e in expr if e not in ('IDs', 'symbols')])
	
	fractions = []
	total_count_1 = 0
	total_count_2 = 0
	
	for s1, s2 in zip(samples[::2], samples[1::2]):
		if o.filter is not None and not any(part in s1.rsplit('_',1)[0] for part in o.filter): continue
		if s1.rsplit('_',1)[0] != s2.rsplit('_',1)[0]: continue # sanity check
		#if not 'c57only' in s1.rsplit('_',1)[0]: continue
		count_only_1 = sum((e1 > 0) and (e2 == 0) for e1,e2 in zip(expr[s1], expr[s2]))
		count_only_2 = sum((e1 == 0) and (e2 > 0) for e1,e2 in zip(expr[s1], expr[s2]))
		count_both = sum(e1 > 0 and e2 > 0 for e1,e2 in zip(expr[s1], expr[s2]))
		if count_only_1+count_only_2+count_both == 0: continue
		total_count_1 += count_only_1
		total_count_2 += count_only_2
		fractions.append((count_only_1+count_only_2)/(count_both+count_only_1+count_only_2))
		print s1.rsplit('_',1)[0], count_only_1, count_only_2, count_both, fractions[-1], count_only_1+count_only_2+count_both
	print 'average', numpy.mean(fractions), 'min', min(fractions), 'max', max(fractions), 'sem', numpy.std(fractions)/math.sqrt(len(fractions))
	print 'total counts', total_count_1, total_count_2

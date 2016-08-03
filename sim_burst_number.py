from __future__ import division
import random, math
print 'bursts', 'mono_fraction', 'log2(mono_fraction)'
for burst_n in range(2,10):
	mono = 0
	bi = 0
	for i in range(100000):
		a1 = 0
		a2 = 0
		for j in range(burst_n):
			if random.random() < 0.5: a1 += 1
			else: a2 += 1
		if a1 and a2: bi += 1
		elif a1 or a2: mono += 1
	print burst_n, mono/(mono+bi), math.log(mono/(mono+bi), 2)

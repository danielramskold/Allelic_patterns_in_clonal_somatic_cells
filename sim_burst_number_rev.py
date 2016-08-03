from __future__ import division
import random, math, pylab, dr_tools

def gives_mono(burst_n):
	a1 = 0
	a2 = 0
	burst_n_int = int(math.ceil(burst_n) if random.random() < burst_n%1 else math.floor(burst_n))
	for j in range(burst_n_int):
		if random.random() < 0.5: a1 += 1
		else: a2 += 1
	return not (a1 and a2)

mono_freq_init = [random.random()**2.5 for i in range(1000)]
mono_freq_out = []

for mono_freq in mono_freq_init:
	burst_n = 1-math.log(mono_freq, 2)
	mono = 0
	bi = 0
	for i in range(29):
		if gives_mono(burst_n): mono += 1
		else: bi += 1
	mono_freq_out.append(mono/(mono+bi))

xarr, yarr = dr_tools.bin(mono_freq_init, 0, 1, 0.05)
pylab.plot(xarr, yarr, 'g-')
xarr, yarr = dr_tools.bin(mono_freq_out, 0, 1, 0.05)
pylab.plot(xarr, yarr, 'b-')
pylab.savefig('plot_sim_burst_rev.pdf')
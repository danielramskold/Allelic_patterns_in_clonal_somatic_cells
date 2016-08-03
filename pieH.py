#http://nxn.se/post/46440196846/making-nicer-looking-pie-charts-with-matplotlib
import numpy as np
import matplotlib.pyplot as plt
import dr_tools
fig = plt.figure(figsize=[10, 10])
ax = fig.add_subplot(111)
cmap = plt.cm.prism
slices = [0.06848306332842416, 0.0007509148987481592, 0.0039718831428571435, 0.09162822808468336, 0.8351659105452872]
colors = cmap(np.linspace(0., 1., len(slices)))
colors = ['#22bb77', '#ee0033', '#00ee33', '#2233aa', '#6677aa']
labels = ['X and Y', 'imprinted mono', 'clonally inherited mono', 'random mono', 'autosomal bi']
ax.pie(slices, colors=colors, labels=['%s %.5f'%(l,f) for l,f in zip(labels, slices)], labeldistance=0.90)
fig.savefig('mono_bi_pie_H.pdf')
